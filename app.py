"""FastAPI backend for the local primer design web tool.

All file handling is in-memory: uploaded sequences are parsed from the
UploadFile stream directly, and outputs are returned in the JSON response.
"""

from __future__ import annotations

import base64
import tempfile
import uuid
from pathlib import Path
from typing import Any

from fastapi import FastAPI, File, HTTPException, UploadFile
from fastapi.responses import HTMLResponse
from fastapi.staticfiles import StaticFiles

from core import (
    COMMON_ENZYMES,
    PROTECTIVE_BASES,
    PROTECTIVE_SEQUENCES,
    design_gibson,
    design_restriction,
    enzyme_info,
    filter_enzymes_for_insert,
    get_enzyme_sites_in_insert,
    get_enzyme_sites_in_vector,
    parse_region,
    pcr_advice,
    read_record,
    revcomp,
    secondary_structure,
)
from genbank_out import write_genbank

BASE = Path(__file__).parent
app = FastAPI(title="Primer Design Agent")


def _read_upload(data: bytes, filename: str) -> dict[str, Any]:
    suffix = Path(filename).suffix.lower() or ".dna"
    with tempfile.NamedTemporaryFile(suffix=suffix, delete=False) as tmp:
        tmp.write(data)
        tmp_path = tmp.name
    try:
        return read_record(tmp_path)
    finally:
        Path(tmp_path).unlink(missing_ok=True)


@app.get("/", response_class=HTMLResponse)
async def index():
    return (BASE / "templates" / "index.html").read_text(encoding="utf-8")


@app.post("/api/upload")
async def upload(
    gene_file: UploadFile = File(...),
    vector_file: UploadFile = File(...),
):
    try:
        gene = _read_upload(await gene_file.read(), gene_file.filename)
        vector = _read_upload(await vector_file.read(), vector_file.filename)

        # Filter enzymes that cut within the insert (gene) sequence
        insert_seq = gene["sequence"]
        vector_seq = vector["sequence"]

        # Enzymes that do NOT cut the insert AND DO cut the vector
        not_in_insert = filter_enzymes_for_insert(insert_seq, COMMON_ENZYMES)
        in_vector = get_enzyme_sites_in_vector(vector_seq, COMMON_ENZYMES)
        available_enzymes = [e for e in not_in_insert if e in in_vector]

        # Enzymes that DO cut the insert (hidden from list)
        filtered_enzymes = get_enzyme_sites_in_insert(insert_seq, COMMON_ENZYMES)
        # Enzymes that do NOT cut the vector (also hidden, different reason)
        not_in_vector = {e: [] for e in COMMON_ENZYMES if e not in in_vector and e not in filtered_enzymes}

        enzyme_details = [
            {
                "name": e,
                **enzyme_info(e),
                "protective_bases": PROTECTIVE_BASES.get(e, 2),
                "protective_sequence": PROTECTIVE_SEQUENCES.get(e, "GCG"),
            }
            for e in available_enzymes
        ]

        return {
            "gene": {
                "name": gene["name"],
                "length": gene["length"],
                "topology": gene["topology"],
                "features": gene["features"],
                "primers_found": gene["primers_found"],
                "sequence": gene["sequence"],
            },
            "vector": {
                "name": vector["name"],
                "length": vector["length"],
                "topology": vector["topology"],
                "features": vector["features"],
                "sequence": vector["sequence"],
            },
            "enzymes": [e["name"] for e in enzyme_details],
            "enzyme_details": enzyme_details,
            "enzymes_filtered": {
                name: {"positions": positions, "count": len(positions)}
                for name, positions in filtered_enzymes.items()
            },
            "enzymes_not_in_vector": sorted(not_in_vector.keys()),
        }
    except Exception as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc


@app.post("/api/design")
async def design(request: dict):
    try:
        gene_seq = request["gene_seq"]
        vector_seq = request["vector_seq"]
        gene_features = request.get("gene_features", [])
        gene_length = request.get("gene_length", len(gene_seq))
        method = request["method"]
        insert_source = request.get("insert_source", "entire")
        insert_strand = request.get("insert_strand", "+")
        enzymes = request.get("enzymes", [])
        vector_insert_pos = request.get("vector_insert_pos")
        homology_len = int(request.get("homology_len", 20))
        target_tm = float(request.get("target_tm", 60.0))
        max_dtm = float(request.get("max_dtm", 2.0))
        min_len = int(request.get("min_len", 18))
        max_len = int(request.get("max_len", 36))

        if insert_source == "entire":
            insert_seq = gene_seq
            insert_start, insert_end = 0, len(gene_seq)
        elif insert_source.startswith("feature:"):
            feat_name = insert_source.split(":", 1)[1]
            feat = next((f for f in gene_features if f["name"] == feat_name), None)
            if not feat:
                raise ValueError(f"Feature {feat_name!r} not found in gene file.")
            insert_seq = gene_seq[feat["start"] : feat["end"]]
            if feat["strand"] == "-":
                insert_seq = revcomp(insert_seq)
            insert_start, insert_end = feat["start"], feat["end"]
        elif insert_source.startswith("region:"):
            insert_start, insert_end = parse_region(insert_source.split(":", 1)[1], gene_length)
            insert_seq = gene_seq[insert_start:insert_end]
        else:
            raise ValueError(f"Unknown insert_source: {insert_source}")

        if insert_strand == "-":
            insert_seq = revcomp(insert_seq)

        if method == "restriction":
            if len(enzymes) != 2:
                raise ValueError("Double digestion requires exactly 2 enzymes.")
            result = design_restriction(
                insert_seq,
                enzyme_f=enzymes[0],
                enzyme_r=enzymes[1],
                extra_5=int(request.get("extra_5", 0)),
                min_len=min_len,
                max_len=max_len,
                target_tm=target_tm,
                max_dtm=max_dtm,
            )
        elif method == "gibson":
            if vector_insert_pos is None:
                raise ValueError("Gibson assembly requires vector_insert_pos.")
            result = design_gibson(
                insert_seq,
                vector_seq=vector_seq,
                insert_pos=int(vector_insert_pos),
                homology_len=homology_len,
                min_len=min_len,
                max_len=max_len,
                target_tm=target_tm,
                max_dtm=max_dtm,
            )
        else:
            raise ValueError(f"Unknown method: {method}")

        if not result["pairs"]:
            raise ValueError(
                "No valid primer pair found. Try widening Tm range, length range, "
                "or max Tm difference."
            )

        best = result["pairs"][0]
        f_data, r_data = best["forward"], best["reverse"]
        f_full = f_data["tail"] + f_data["anneal"]
        r_full = r_data["tail"] + r_data["anneal"]
        amplicon = f_data["tail"] + insert_seq + revcomp(r_data["tail"])

        job_id = uuid.uuid4().hex[:8]
        method_label = result.get("method", "design")
        primer_entries = [
            {
                "name": f"{method_label}_F{f_data['length']}",
                "sequence": f_full,
                "start": 0,
                "end": f_data["length"] + len(f_data["tail"]),
                "strand": "+",
            },
            {
                "name": f"{method_label}_R{r_data['length']}",
                "sequence": r_full,
                "start": len(amplicon) - r_data["length"] - len(r_data["tail"]),
                "end": len(amplicon),
                "strand": "-",
            },
        ]
        gb_text = write_genbank(
            name=f"primer_design_{job_id}",
            sequence=amplicon,
            topology="linear",
            primers=primer_entries,
        )
        gb_b64 = base64.b64encode(gb_text.encode("utf-8")).decode("ascii")

        order_lines = ["Name\tSequence (5'->3')\tLength\tTm (C)\tGC (%)"]
        for label, data in [("F", f_data), ("R", r_data)]:
            seq = data["tail"] + data["anneal"]
            order_lines.append(
                f"{method_label}_{label}{data['length']}\t{seq}\t"
                f"{len(seq)}\t{data['tm']:.1f}\t{data['gc']}"
            )
        tsv_text = "\n".join(order_lines) + "\n"
        tsv_b64 = base64.b64encode(tsv_text.encode("utf-8")).decode("ascii")

        checks = secondary_structure(f_full, r_full)
        advice = pcr_advice(
            insert_seq,
            f_data["tm"],
            r_data["tm"],
            checks,
            method_label,
        )

        alts = []
        for alt in result["pairs"][1:6]:
            af, ar = alt["forward"], alt["reverse"]
            alts.append(
                {
                    "forward": af["tail"] + af["anneal"],
                    "reverse": ar["tail"] + ar["anneal"],
                    "tm_f": round(af["tm"], 1),
                    "tm_r": round(ar["tm"], 1),
                    "dtm": round(alt["dtm"], 1),
                    "len_f": af["length"],
                    "len_r": ar["length"],
                }
            )

        return {
            "ok": True,
            "job_id": job_id,
            "method": method_label,
            "insert_length": len(insert_seq),
            "insert_location": f"{insert_start + 1}..{insert_end}",
            "target_tm": target_tm,
            "tm_standard": "estimate:primer3-nn (50 mM Na+, no Mg2+)",
            "best": {
                "forward": {
                    "name": f"F{f_data['length']}",
                    "sequence": f_full,
                    "tail": f_data["tail"],
                    "anneal": f_data["anneal"],
                    "length": len(f_full),
                    "anneal_length": f_data["length"],
                    "tm": round(f_data["tm"], 1),
                    "gc": f_data["gc"],
                },
                "reverse": {
                    "name": f"R{r_data['length']}",
                    "sequence": r_full,
                    "tail": r_data["tail"],
                    "anneal": r_data["anneal"],
                    "length": len(r_full),
                    "anneal_length": r_data["length"],
                    "tm": round(r_data["tm"], 1),
                    "gc": r_data["gc"],
                },
                "tm_difference": round(abs(f_data["tm"] - r_data["tm"]), 1),
                "product_length": len(amplicon),
                "checks": checks,
            },
            "pcr_advice": advice,
            "alternatives": alts,
            "warnings": result.get("warnings", []),
            "output_gb_b64": gb_b64,
            "output_tsv_b64": tsv_b64,
            "output_gb_name": f"primer_design_{job_id}.gb",
            "output_tsv_name": f"primer_design_{job_id}_order.tsv",
            "amplicon_preview": amplicon[:60] + ("..." if len(amplicon) > 60 else ""),
        }
    except Exception as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc


@app.get("/api/enzymes")
async def list_enzymes():
    return [{"name": e, **enzyme_info(e)} for e in COMMON_ENZYMES]


app.mount("/static", StaticFiles(directory=str(BASE / "static")), name="static")
