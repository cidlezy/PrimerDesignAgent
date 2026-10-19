"""Write SnapGene-compatible GenBank files with primer_bind features.

Based on snapgene-bridge's snapgene_genbank.py (Apache-2.0).
SnapGene imports primer_bind features as real primers when the JOURNAL line
mentions "from SnapGene".
"""

from __future__ import annotations

import datetime
import re

INDENT = " " * 21
WRAP = 58
LOCUS_RE = re.compile(r"[^A-Za-z0-9_.\-]")
FEATURE_KEY_RE = re.compile(r"^[A-Za-z0-9_'\-*]{1,15}$")


def _clean(text, limit=120):
    s = re.sub(r"[\r\n\t\"]", " ", str(text or "")).strip()
    return s[:limit] or "unnamed"


def _qual(key, value, *, quote=True):
    body = f'/{key}="{value}"' if quote else f"/{key}={value}"
    lines = []
    cur = ""
    for word in body.split(" "):
        cand = f"{cur} {word}" if cur else word
        if cur and len(cand) > WRAP:
            lines.append(cur)
            cur = word
        else:
            cur = cand
    lines.append(cur)
    return "".join(f"{INDENT}{line}\n" for line in lines)


def _feature_key(value):
    key = str(value or "misc_feature").strip().replace(" ", "_")
    return key if FEATURE_KEY_RE.match(key) else "misc_feature"


def write_genbank(
    *,
    name: str,
    sequence: str,
    topology: str = "linear",
    features=None,
    primers=None,
    snapgene_version: str = "8.0.0",
) -> str:
    """Return GenBank text that SnapGene imports with primer annotations."""

    sequence = sequence.upper()
    length = len(sequence)
    today = datetime.date.today()
    locus = (LOCUS_RE.sub("_", name) or "sequence")[:16]
    shape = "circular" if topology == "circular" else "linear  "

    out = [
        f"LOCUS       {locus:<16} {length:>11} bp    DNA     {shape} SYN "
        f"{today.strftime('%d-%b-%Y').upper()}\n",
        f"DEFINITION  {_clean(name, 60)}.\n",
        "ACCESSION   .\n",
        "VERSION     .\n",
        "KEYWORDS    .\n",
        "SOURCE      synthetic DNA construct\n",
        "  ORGANISM  synthetic DNA construct\n",
        f"REFERENCE   1  (bases 1 to {length})\n",
        "  AUTHORS   PrimerDesignAgent\n",
        "  TITLE     Direct Submission\n",
        f"  JOURNAL   Exported {today.strftime('%b')} {today.day}, {today.year} "
        f"from SnapGene {snapgene_version}\n",
        "            https://www.snapgene.com\n",
        "FEATURES             Location/Qualifiers\n",
        f"     source          1..{length}\n",
        _qual("mol_type", "other DNA"),
        _qual("organism", "synthetic DNA construct"),
    ]

    for feat in features or []:
        start = int(feat.get("start", 0))
        end = int(feat.get("end", 0))
        if end <= start:
            continue
        strand = feat.get("strand", "+")
        loc = f"{start + 1}..{end}"
        if strand == "-":
            loc = f"complement({loc})"
        out.append(f"     {_feature_key(feat.get('type')):<16}{loc}\n")
        out.append(_qual("label", _clean(feat.get("name", "feature")), quote=False))

    for primer in primers or []:
        start = int(primer.get("start", 0))
        end = int(primer.get("end", min(len(primer.get("sequence", "")), length)))
        strand = primer.get("strand", "+")
        loc = f"{start + 1}..{end}"
        if strand == "-":
            loc = f"complement({loc})"
        out.append(f"     primer_bind     {loc}\n")
        out.append(_qual("label", _clean(primer.get("name", "primer")), quote=False))
        out.append(_qual("note", f"color: black; sequence: {primer['sequence'].upper()}"))

    out.append("ORIGIN\n")
    for offset in range(0, length, 60):
        chunk = " ".join(sequence[i : i + 10] for i in range(offset, min(offset + 60, length), 10))
        out.append(f"{offset + 1:>9} {chunk}\n")
    out.append("//\n")
    return "".join(out)
