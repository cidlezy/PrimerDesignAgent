"""Core primer design engine for restriction cloning and Gibson/HiFi assembly.

Uses Biopython for sequence I/O and restriction analysis, and Primer3-py for
melting temperature calculation. Outputs GenBank files compatible with
SnapGene's primer importer.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

from Bio import Restriction, SeqIO
from Bio.Seq import Seq
import primer3

COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def read_record(path: str | Path) -> dict[str, Any]:
    """Read one sequence record from .dna, GenBank, or FASTA."""

    src = Path(path)
    suffix = src.suffix.lower()
    if suffix == ".dna":
        fmt = "snapgene"
    elif suffix in (".gb", ".gbk", ".genbank"):
        fmt = "genbank"
    elif suffix in (".fa", ".fasta", ".fna"):
        fmt = "fasta"
    else:
        raise ValueError(f"Unsupported file format: {suffix}")

    record = SeqIO.read(str(src), fmt)
    features = []
    primers_found = []
    for feat in record.features:
        if feat.type == "source":
            continue
        label = (
            feat.qualifiers.get("label", [None])[0]
            or feat.qualifiers.get("gene", [None])[0]
            or feat.qualifiers.get("product", [None])[0]
            or feat.type
        )
        features.append(
            {
                "name": str(label),
                "type": feat.type,
                "location": str(feat.location),
                "start": int(feat.location.start),
                "end": int(feat.location.end),
                "strand": "+" if (feat.location.strand or 1) >= 0 else "-",
            }
        )
        if feat.type == "primer_bind":
            primers_found.append(
                {
                    "name": str(label),
                    "location": str(feat.location),
                    "start": int(feat.location.start),
                    "end": int(feat.location.end),
                    "strand": "+" if (feat.location.strand or 1) >= 0 else "-",
                    "sequence": "",
                }
            )

    # For .dna files, use sgffp to get full primer sequences and Tm
    if fmt == "snapgene":
        try:
            import sgffp as sgffp_mod
            sgff_reader = sgffp_mod.SgffReader.from_file(str(src))
            if sgff_reader.has_primers:
                primers_found = []
                for p in sgff_reader.primers:
                    best_site = max(
                        (s for s in p.binding_sites if not s.simplified),
                        key=lambda s: s.melting_temperature or 0,
                        default=None,
                    )
                    primers_found.append(
                        {
                            "name": p.name,
                            "location": f"{best_site.start + 1}..{best_site.end}" if best_site else "",
                            "start": best_site.start if best_site else 0,
                            "end": best_site.end if best_site else 0,
                            "strand": best_site.bound_strand if best_site else "+",
                            "sequence": p.sequence,
                            "tm": best_site.melting_temperature if best_site else None,
                            "anneal_length": len(best_site.annealed_bases) if best_site else 0,
                        }
                    )
        except ImportError:
            pass

    topology = "circular" if record.annotations.get("topology") == "circular" else "linear"
    return {
        "name": record.id or src.stem,
        "sequence": str(record.seq).upper(),
        "length": len(record.seq),
        "topology": topology,
        "features": features,
        "primers_found": primers_found,
        "format": fmt,
        "path": str(src),
    }


COMMON_ENZYMES = [
    "EcoRI", "BamHI", "XhoI", "SalI", "KpnI", "NdeI", "NotI", "HindIII",
    "EcoRV", "SacI", "SmaI", "XbaI", "SpeI", "PstI", "SphI", "NheI",
    "AgeI", "BglII", "ClaI", "ApaI", "NcoI", "NruI", "MluI", "AflII",
]

# NEB recommended extra bases 5' of the recognition site for efficient digestion.
# Actual 5' protective sequences based on NEB cleavage data and common lab practice.
# Each entry is the literal sequence to place 5' of the recognition site.
PROTECTIVE_SEQUENCES = {
    "EcoRI":   "GCG",
    "BamHI":   "GCG",
    "XhoI":    "CCG",
    "SalI":    "GCG",
    "KpnI":    "GCG",
    "NdeI":    "GCG",
    "NotI":    "GCGACC",
    "HindIII": "CCC",
    "EcoRV":   "GCG",
    "SacI":    "CGC",
    "SmaI":    "AAAC",
    "XbaI":    "GCG",
    "SpeI":    "GCG",
    "PstI":    "GCG",
    "SphI":    "GCG",
    "NheI":    "GCG",
    "AgeI":    "GCG",
    "BglII":   "GCG",
    "ClaI":    "GCG",
    "ApaI":    "GCG",
    "NcoI":    "GCG",
    "NruI":    "GCG",
    "MluI":    "GCG",
    "AflII":   "GCG",
}

PROTECTIVE_BASES = {k: len(v) for k, v in PROTECTIVE_SEQUENCES.items()}


def filter_enzymes_for_insert(sequence: str, enzymes: list[str]) -> dict[str, list[int]]:
    """Return only enzymes whose recognition sites are absent from the insert."""

    seq = sequence.upper()
    out = {}
    for name in enzymes:
        sites = find_sites(seq, [name]).get(name, [])
        if not sites:
            out[name] = []
    return out


def get_enzyme_sites_in_insert(sequence: str, enzymes: list[str]) -> dict[str, list[int]]:
    """Return enzymes that DO cut the insert and their positions."""

    seq = sequence.upper()
    return {name: sites for name, sites in find_sites(seq, enzymes).items() if sites}


def get_enzyme_sites_in_vector(sequence: str, enzymes: list[str]) -> dict[str, list[int]]:
    """Return enzymes that DO cut the vector and their positions."""

    return get_enzyme_sites_in_insert(sequence, enzymes)

def enzyme_info(name: str) -> dict[str, Any]:
    """Return recognition site and cut position for a named enzyme."""

    enzyme = getattr(Restriction, name, None)
    if enzyme is None:
        raise ValueError(f"Unknown enzyme: {name}")
    return {
        "name": enzyme.__name__,
        "site": enzyme.site,
        "size": enzyme.size,
        "fst5": enzyme.fst5,
        "fst3": enzyme.fst3,
        "ovhg": enzyme.ovhg,
    }


def find_sites(sequence: str, names: list[str]) -> dict[str, list[int]]:
    """Find 1-based positions of enzyme recognition sites in a sequence."""

    seq = Seq(sequence)
    batch = Restriction.RestrictionBatch([getattr(Restriction, n) for n in names])
    raw = batch.search(seq)
    out = {}
    for name in names:
        positions = raw.get(getattr(Restriction, name), [])
        out[name] = sorted(positions)
    return out


def gc_percent(seq: str) -> float:
    if not seq:
        return 0.0
    return 100.0 * sum(base in "GC" for base in seq) / len(seq)


def calc_tm(seq: str) -> float:
    """Nearest-neighbour Tm at 50 mM Na+, no Mg2+ (SantaLucia parameters)."""

    return primer3.calc_tm(seq, mv_conc=50.0, dv_conc=0.0, dntp_conc=0.0, dna_conc=250.0)


def secondary_structure(forward: str, reverse: str) -> dict[str, Any]:
    """Hairpin and dimer check using Primer3."""

    cond = dict(mv_conc=50.0, dv_conc=0.0, dntp_conc=0.0, dna_conc=250.0)
    return {
        "forward_hairpin_tm": round(primer3.calc_hairpin(forward, **cond).tm, 1),
        "reverse_hairpin_tm": round(primer3.calc_hairpin(reverse, **cond).tm, 1),
        "forward_self_dimer_dg": round(primer3.calc_homodimer(forward, **cond).dg / 1000, 1),
        "reverse_self_dimer_dg": round(primer3.calc_homodimer(reverse, **cond).dg / 1000, 1),
        "pair_dimer_dg": round(primer3.calc_heterodimer(forward, reverse, **cond).dg / 1000, 1),
    }


def pcr_advice(
    insert: str,
    tm_f: float,
    tm_r: float,
    checks: dict[str, Any],
    method: str,
) -> dict[str, Any]:
    """Generate amplification protocol suggestions based on sequence properties."""

    insert = insert.upper()
    gc = gc_percent(insert)
    length = len(insert)
    lower_tm = min(tm_f, tm_r)
    higher_tm = max(tm_f, tm_r)
    anneal = round(lower_tm - 3, 1)

    items: list[str] = []
    flags: list[str] = []

    # GC content
    if gc >= 70:
        flags.append("high_gc")
        items.append(
            f"插入片段 GC 含量 {gc:.0f}%，非常高。建议使用 Touchdown PCR（起始退火温度 "
            f"{round(lower_tm + 2, 0):.0f} °C，每 2 个循环降 1 °C，降至 {anneal:.0f} °C 后再扩 15-20 循环）。"
        )
        items.append("PCR 体系加 3-5% DMSO 或 1-1.5 M 甜菜碱（betaine）帮助打开高 GC 区域。")
        items.append("变性温度提高到 98 °C（如果聚合酶允许），变性时间延长到 10-15 s。")
        items.append("优先选择高 GC 模板专用聚合酶，如 Phusion HF、Q5 或 KOD FX。")
    elif gc >= 65:
        flags.append("high_gc")
        items.append(
            f"插入片段 GC 含量 {gc:.0f}%，偏高。建议退火温度从 {round(lower_tm + 2, 0):.0f} °C 起，"
            f"用 Touchdown 降到 {anneal:.0f} °C，或直接用 {anneal:.0f} °C 并加 3% DMSO。"
        )
    elif gc <= 35:
        flags.append("low_gc")
        items.append(
            f"插入片段 GC 含量仅 {gc:.0f}%，偏低。退火温度可以用 {max(50, round(anneal, 0)):.0f} °C 起试，"
            f"不必 Touchdown。如果非特异性条带多，再逐步升高退火温度。"
        )
    else:
        items.append(
            f"插入片段 GC 含量 {gc:.0f}%，在正常范围。退火温度建议 {anneal:.0f} °C。"
        )

    # Tm mismatch
    if higher_tm - lower_tm > 5:
        flags.append("tm_mismatch")
        items.append(
            f"两条引物 Tm 相差 {higher_tm - lower_tm:.1f} °C，比较大。"
            f"按较低的 Tm 设退火温度（{anneal:.0f} °C），但效率可能偏低，"
            f"可以考虑加 5-10% DMSO 或用 Touchdown 补偿。"
        )

    # Product length
    if length > 3000:
        flags.append("long_product")
        items.append(
            f"产物长度 {length} bp，较长。延伸时间按 30 s/kb 计算，约 {round(length / 1000 * 30, 0):.0f} s；"
            f"建议用高保真长片段聚合酶（如 Q5、Phusion、PrimeSTAR GXL）。"
        )
    elif length > 1500:
        flags.append("medium_product")
        items.append(
            f"产物长度 {length} bp。延伸时间按 30 s/kb，约 {round(length / 1000 * 30, 0):.0f} s。"
        )
    elif length < 200:
        flags.append("short_product")
        items.append(
            f"产物长度仅 {length} bp，延伸时间 15-20 s 足够。退火温度可稍高（{anneal + 2:.0f} °C）以减少非特异。"
        )

    # Secondary structure
    if checks.get("forward_hairpin_tm", 0) > 50 or checks.get("reverse_hairpin_tm", 0) > 50:
        flags.append("hairpin")
        items.append(
            f"引物有较强 hairpin（Tm {max(checks.get('forward_hairpin_tm', 0), checks.get('reverse_hairpin_tm', 0)):.1f} °C）。"
            f"建议加 5% DMSO，或在设计备选引物对中换用 hairpin 更弱的一条。"
        )
    if checks.get("pair_dimer_dg", 0) < -10:
        flags.append("dimer")
        items.append(
            f"两条引物之间有较强 dimer（ΔG {checks.get('pair_dimer_dg', 0):.1f} kcal/mol）。"
            f"建议优先看备选引物对，或加 5% DMSO 抑制。"
        )

    # Method-specific advice
    if method == "gibson":
        items.append("Gibson/HiFi 组装：同源臂区域 PCR 时，如产物较短，延伸时间可以适当缩短。")
    elif method == "restriction":
        items.append("双酶切克隆：PCR 产物建议先跑胶确认单一条带，再切胶回收做双酶切。")

    return {
        "insert_gc_percent": round(gc, 1),
        "insert_length": length,
        "recommended_anneal_c": anneal,
        "flags": flags,
        "advice": items,
    }


def parse_region(text: str, length: int) -> tuple[int, int]:
    """Parse 1-based inclusive region like '1626..2486' to 0-based half-open."""

    cleaned = text.strip().replace(",", "")
    for sep in ("..", "-", ":"):
        if sep in cleaned:
            left, right = cleaned.split(sep, 1)
            break
    else:
        raise ValueError(f"Cannot parse region: {text!r}")
    start, end = int(left) - 1, int(right)
    if not (0 <= start < end <= length):
        raise ValueError(f"Region {text!r} is outside 1..{length}.")
    return start, end


def design_candidates(
    insert: str,
    *,
    tail_f: str,
    tail_r: str,
    min_len: int = 18,
    max_len: int = 36,
    target_tm: float = 60.0,
    max_dtm: float = 2.0,
) -> dict[str, Any]:
    """Enumerate annealing lengths at both insert ends and pick the best pair.

    tail_f and tail_r are the 5' additions (enzyme sites, homology arms, etc.)
    on the forward and reverse primers respectively. tail_r should already be
    written as it appears on the reverse primer strand (i.e. if the vector
    has "AAAA..." downstream, tail_r = revcomp("AAAA...")).
    """

    insert = insert.upper()
    tail_f = tail_f.upper()
    tail_r = tail_r.upper()
    if len(insert) < min_len:
        raise ValueError(f"Insert is {len(insert)} bp, shorter than min primer length {min_len}.")

    forward: list[dict[str, Any]] = []
    reverse: list[dict[str, Any]] = []
    for length in range(min_len, min(max_len, len(insert)) + 1):
        f_anneal = insert[:length]
        r_anneal = revcomp(insert[-length:])
        forward.append(
            {
                "tail": tail_f,
                "anneal": f_anneal,
                "tm": calc_tm(f_anneal),
                "gc": round(gc_percent(f_anneal), 1),
                "length": length,
            }
        )
        reverse.append(
            {
                "tail": tail_r,
                "anneal": r_anneal,
                "tm": calc_tm(r_anneal),
                "gc": round(gc_percent(r_anneal), 1),
                "length": length,
            }
        )

    scored: list[dict[str, Any]] = []
    all_scored: list[dict[str, Any]] = []
    for f in forward:
        for r in reverse:
            dtm = abs(f["tm"] - r["tm"])
            mean_tm = (f["tm"] + r["tm"]) / 2
            tm_dev = abs(mean_tm - target_tm)
            cost = tm_dev * 2.0 + dtm * 0.5 + 0.01 * f["length"] + 0.01 * r["length"]
            if f["anneal"][-1] not in "GC":
                cost += 0.05
            if r["anneal"][-1] not in "GC":
                cost += 0.05
            entry = {
                "cost": cost,
                "forward": f,
                "reverse": r,
                "dtm": dtm,
                "tm_dev": tm_dev,
                "mean_tm": mean_tm,
            }
            all_scored.append(entry)
            if dtm <= max_dtm and tm_dev <= max_dtm * 2:
                scored.append(entry)

    if not scored:
        # Report the closest achievable pair so the user understands why.
        all_scored.sort(key=lambda x: x["cost"])
        closest = all_scored[0] if all_scored else None
        detail = ""
        if closest:
            cf, cr = closest["forward"], closest["reverse"]
            detail = (
                f"最接近的可实现引物对: Tm {cf['tm']:.1f} / {cr['tm']:.1f} °C "
                f"(均值 {closest['mean_tm']:.1f} °C，偏离目标 {closest['tm_dev']:.1f} °C)，"
                f"Forward GC {cf['gc']}% / Reverse GC {cr['gc']}%. "
                f"这条序列两端 GC 含量极不对称（5' 端 GC {forward[0]['gc']}%，3' 端 GC {reverse[0]['gc']}%），"
                f"在最大 Tm 差 {max_dtm} 的约束下找不到两端 Tm 都接近 {target_tm} °C 的引物。"
            )
        raise ValueError(
            "在给定的 Tm 约束下找不到引物对。"
            + detail
            + " 建议把'最大 Tm 差'放宽到 6-10 °C，或接受按较高 Tm 设退火温度。"
        )

    scored.sort(key=lambda x: x["cost"])
    return {
        "forward": forward,
        "reverse": reverse,
        "pairs": scored,
        "closest_achievable_tm_dev": min((e["tm_dev"] for e in all_scored), default=None),
    }


def design_restriction(
    insert: str,
    *,
    enzyme_f: str,
    enzyme_r: str,
    extra_5: int = 0,
    min_len: int = 18,
    max_len: int = 36,
    target_tm: float = 60.0,
    max_dtm: float = 2.0,
) -> dict[str, Any]:
    """Design primers with two enzyme sites as 5' tails."""

    info_f = enzyme_info(enzyme_f)
    info_r = enzyme_info(enzyme_r)
    site_f = info_f["site"]
    site_r = info_r["site"]

    # Actual protective sequences are selected from the local NEB-derived table.
    # `extra_5` is retained for API compatibility; a zero/positive value does not
    # replace these real bases because arbitrary N bases cannot be ordered.
    protect_f = PROTECTIVE_SEQUENCES.get(enzyme_f, "GCG")
    protect_r = PROTECTIVE_SEQUENCES.get(enzyme_r, "GCG")

    tail_f = protect_f + site_f
    tail_r = protect_r + revcomp(site_r)

    result = design_candidates(
        insert,
        tail_f=tail_f,
        tail_r=tail_r,
        min_len=min_len,
        max_len=max_len,
        target_tm=target_tm,
        max_dtm=max_dtm,
    )
    result["method"] = "restriction"
    result["enzyme_f"] = info_f
    result["enzyme_r"] = info_r
    result["protective_f"] = protect_f
    result["protective_r"] = protect_r
    result["tail_f"] = tail_f
    result["tail_r"] = tail_r
    result["warnings"] = []
    best = result["pairs"][0] if result["pairs"] else None
    if best and best["tm_dev"] > 3.0:
        bf, br = best["forward"], best["reverse"]
        result["warnings"].append(
            f"最优引物对 Tm {bf['tm']:.1f} / {br['tm']:.1f} °C，平均 {best['mean_tm']:.1f} °C，"
            f"偏离目标 {target_tm:.0f} °C 有 {best['tm_dev']:.1f} °C。"
            f"这通常是因为插入片段端部 GC 含量过高（Forward GC {bf['gc']}%），"
            f"即使最短引物（{bf['length']} nt）Tm 也已经超过目标。"
            f"建议：按较高 Tm 设退火温度，或使用 Touchdown PCR 从 {round(bf['tm'] + 2, 0):.0f} °C 开始。"
        )
    result["closest_achievable_tm_dev"] = result.get("closest_achievable_tm_dev")
    # Check that enzyme sites are not inside the insert
    sites_in_insert = find_sites(insert, [enzyme_f, enzyme_r])
    for enz, positions in sites_in_insert.items():
        if positions:
            result["warnings"].append(
                f"{enz} site found inside the insert at position(s) {positions}; "
                "digestion will cut the insert too."
            )
    return result


def design_gibson(
    insert: str,
    *,
    vector_seq: str,
    insert_pos: int,
    homology_len: int = 20,
    min_len: int = 18,
    max_len: int = 36,
    target_tm: float = 60.0,
    max_dtm: float = 2.0,
) -> dict[str, Any]:
    """Design primers with vector homology arms for Gibson/HiFi assembly.

    insert_pos is 0-based: the insert goes between vector[insert_pos-1] and
    vector[insert_pos] (i.e. after position insert_pos in 0-based terms).
    """

    vector_seq = vector_seq.upper()
    n = len(vector_seq)
    if not (0 <= insert_pos <= n):
        raise ValueError(f"Insertion position {insert_pos} is outside 0..{n}.")

    arm_f_start = max(0, insert_pos - homology_len)
    arm_f = vector_seq[arm_f_start:insert_pos]
    arm_r_end = min(n, insert_pos + homology_len)
    arm_r_top = vector_seq[insert_pos:arm_r_end]
    arm_r = revcomp(arm_r_top)

    result = design_candidates(
        insert,
        tail_f=arm_f,
        tail_r=arm_r,
        min_len=min_len,
        max_len=max_len,
        target_tm=target_tm,
        max_dtm=max_dtm,
    )
    result["method"] = "gibson"
    result["homology_len"] = homology_len
    result["insert_pos"] = insert_pos
    result["arm_f"] = arm_f
    result["arm_r"] = arm_r
    result["warnings"] = []
    if len(arm_f) < homology_len:
        result["warnings"].append(
            f"Forward homology arm shortened to {len(arm_f)} bp (near vector start)."
        )
    if len(arm_r) < homology_len:
        result["warnings"].append(
            f"Reverse homology arm shortened to {len(arm_r)} bp (near vector end)."
        )
    return result
