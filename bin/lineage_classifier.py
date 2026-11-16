#!/usr/bin/env python3

import argparse
import gzip
import json
import math
from datetime import datetime

PRESENT_FRACTION = 0.5
MINOR_FRACTION = 0.1
HIGH_FRACTION = 0.9
MEDIUM_FRACTION = 0.7
GENOME_WIDE_SNPS = 50
BASES = "ACGT"


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def load_barcode_bed(path):
    markers, families = [], {}
    with open(path) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            lineage = fields[3]
            markers.append({"chrom": fields[0], "pos": int(fields[2]), "lineage": lineage,
                            "allele": fields[4].upper()})
            families.setdefault(lineage, fields[5] if len(fields) > 5 and fields[5] else "Unknown")
    return markers, families


def parent_of(lineage, known):
    name = lineage
    while "." in name:
        name = name.rsplit(".", 1)[0]
        if name in known:
            return name
    return None


def ancestors_of(lineage, known):
    chain = []
    parent = parent_of(lineage, known)
    while parent:
        chain.append(parent)
        parent = parent_of(parent, known)
    return chain


def read_reference(path):
    sequences, name, chunks = {}, None, []
    with open_text(path) as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    sequences[name] = "".join(chunks).upper()
                name, chunks = line[1:].split()[0], []
            elif line:
                chunks.append(line)
    if name is not None:
        sequences[name] = "".join(chunks).upper()
    return sequences


def counts_to_fractions(counts):
    depth = sum(counts.values())
    if depth <= 0:
        return None, 0
    return {base: count / depth for base, count in counts.items()}, depth


def read_pileup(path, wanted):
    observed = {}
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5 or line.startswith("#"):
                continue
            chrom, pos, ref, alts, ad = fields[0], int(fields[1]), fields[2].upper(), fields[3], fields[4]
            if (chrom, pos) not in wanted or ad in ("", "."):
                continue
            alleles = [ref] + [a.upper() for a in alts.split(",")]
            values = [int(v) if v.isdigit() else 0 for v in ad.split(",")]
            counts = {}
            for allele, value in zip(alleles, values):
                if allele in BASES:
                    counts[allele] = counts.get(allele, 0) + value
            fractions, depth = counts_to_fractions(counts)
            if fractions:
                observed[(chrom, pos)] = {"fractions": fractions, "depth": depth, "inferred": False}
    return observed


def genotype_fractions(ref, alts, fmt, sample):
    keys = fmt.split(":")
    values = dict(zip(keys, sample.split(":")))
    alleles = [ref] + alts
    ad = values.get("AD", "")
    if ad and "." not in ad.split(","):
        counts = {}
        for allele, value in zip(alleles, ad.split(",")):
            counts[allele] = counts.get(allele, 0) + int(value)
        fractions, depth = counts_to_fractions(counts)
        if fractions:
            return fractions, depth
    gt = values.get("GT", "").replace("|", "/")
    calls = [int(g) for g in gt.split("/") if g.isdigit() and int(g) < len(alleles)]
    if not calls:
        return None, 0
    counts = {}
    for call in calls:
        counts[alleles[call]] = counts.get(alleles[call], 0) + 1
    fractions, _ = counts_to_fractions(counts)
    depth = values.get("DP", "")
    return fractions, int(depth) if depth.isdigit() else None


def read_vcf(path, wanted):
    observed, snps = {}, 0
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            chrom, pos, ref = fields[0], int(fields[1]), fields[3].upper()
            alts = [a.upper() for a in fields[4].split(",")]
            if len(ref) != 1 or not any(len(a) == 1 and a in BASES for a in alts):
                continue
            snps += 1
            if (chrom, pos) not in wanted:
                continue
            fractions, depth = (genotype_fractions(ref, alts, fields[8], fields[9]) if len(fields) > 9
                                else ({alts[0]: 1.0}, None))
            if fractions:
                observed[(chrom, pos)] = {"fractions": fractions, "depth": depth, "inferred": False}
    return observed, snps


def marker_evidence(markers, observed, reference, infer_reference):
    evidence = {}
    for marker in markers:
        key = (marker["chrom"], marker["pos"])
        if key in observed or key in evidence:
            continue
        if infer_reference:
            seq = reference.get(marker["chrom"], "")
            if 0 < marker["pos"] <= len(seq):
                evidence[key] = {"fractions": {seq[marker["pos"] - 1]: 1.0}, "depth": None, "inferred": True}
    evidence.update(observed)
    return evidence


def lineage_stats(markers, evidence):
    grouped = {}
    for marker in markers:
        grouped.setdefault(marker["lineage"], []).append(marker)
    stats = {}
    for lineage, group in grouped.items():
        fractions, details, inferred = [], [], 0
        for marker in group:
            item = evidence.get((marker["chrom"], marker["pos"]))
            if item is None:
                continue
            fraction = item["fractions"].get(marker["allele"], 0.0)
            fractions.append(fraction)
            if fraction >= PRESENT_FRACTION:
                details.append({"pos": marker["pos"], "allele": marker["allele"], "fraction": round(fraction, 3),
                                "depth": item["depth"], "inferred": item["inferred"]})
                inferred += item["inferred"]
        needed = max(1, math.ceil(len(group) / 2))
        mean = sum(fractions) / len(fractions) if fractions else None
        if len(fractions) < needed:
            state = "no_data"
        elif mean >= PRESENT_FRACTION:
            state = "present"
        else:
            state = "absent"
        stats[lineage] = {"state": state, "fraction": mean, "markers": len(group), "callable": len(fractions),
                          "supporting": len(details), "inferred": inferred, "details": details,
                          "minor_markers": sum(f >= MINOR_FRACTION for f in fractions)}
    return stats


def lowest_common_ancestor(lineages, known):
    paths = [[name] + ancestors_of(name, known) for name in lineages]
    shared = set(paths[0]).intersection(*paths[1:])
    for name in paths[0]:
        if name in shared:
            return name
    return None


def confidence_for(call, stats, known, ambiguous, mixed):
    path = [call] + ancestors_of(call, known)
    if ambiguous:
        return "low"
    fraction = stats[call]["fraction"] or 0.0
    inferred = any(stats[name]["inferred"] for name in path)
    strong = all(stats[name]["state"] == "present" and stats[name]["fraction"] >= HIGH_FRACTION for name in path)
    enough = stats[call]["supporting"] >= min(3, stats[call]["markers"])
    if strong and enough and not inferred and not mixed:
        return "high"
    if fraction >= MEDIUM_FRACTION:
        return "medium"
    return "low"


def unknown_result(reason, evidence_type, min_depth, stats=None):
    return {"lineage": "unknown", "confidence": "low", "score": 0, "matched_snps": 0, "total_snps": 0,
            "percentage": 0, "family": "Unknown", "lineage_path": [], "matched_details": [],
            "all_scores": scores_summary(stats or {}), "evidence": evidence_type, "min_depth": min_depth,
            "possible_mixed_infection": False, "minor_lineages": [], "ambiguous_sublineages": [],
            "reason": reason}


def scores_summary(stats):
    return {name: {"state": s["state"], "fraction": round(s["fraction"], 3) if s["fraction"] is not None else None,
                   "supporting": s["supporting"], "callable": s["callable"], "markers": s["markers"]}
            for name, s in sorted(stats.items()) if s["callable"] and s["supporting"]}


def classify_lineage(markers, families, evidence, evidence_type="pileup", min_depth=None):
    known = set(families)
    stats = lineage_stats(markers, evidence)
    if not any(s["callable"] for s in stats.values()):
        return unknown_result("no barcode site had usable data", evidence_type, min_depth, stats)

    def consistent(name):
        return all(stats[a]["state"] != "absent" for a in ancestors_of(name, known))

    candidates = sorted(name for name, s in stats.items() if s["state"] == "present" and consistent(name))
    if not candidates:
        return unknown_result("no lineage had most of its barcode markers", evidence_type, min_depth, stats)
    leaves = [name for name in candidates
              if not any(other != name and name in ancestors_of(other, known) for other in candidates)]
    ambiguous = []
    call = leaves[0]
    if len(leaves) > 1:
        ambiguous = leaves
        call = lowest_common_ancestor(leaves, known)
        while call is not None and stats.get(call, {}).get("state") != "present":
            call = parent_of(call, known)
        if call is None:
            result = unknown_result("barcode markers of unrelated lineages are all present: " + ", ".join(leaves),
                                    evidence_type, min_depth, stats)
            result["ambiguous_sublineages"] = leaves
            result["possible_mixed_infection"] = True
            return result

    path = [call] + ancestors_of(call, known)
    minors = sorted(
        ({"lineage": name, "fraction": round(s["fraction"], 3)} for name, s in stats.items()
         if name not in path and s["state"] == "absent" and s["fraction"] >= MINOR_FRACTION
         and s["minor_markers"] >= min(2, s["markers"]) and consistent_minor(name, stats, known)),
        key=lambda item: (-item["fraction"], item["lineage"]))
    mixed = bool(minors) or bool(ambiguous)
    best = stats[call]
    result = {
        "lineage": call,
        "confidence": confidence_for(call, stats, known, bool(ambiguous), mixed),
        "score": round(best["fraction"], 3),
        "matched_snps": best["supporting"],
        "total_snps": best["markers"],
        "percentage": round(best["fraction"] * 100, 2),
        "family": families.get(call, "Unknown"),
        "lineage_path": list(reversed(path)),
        "matched_details": best["details"],
        "all_scores": scores_summary(stats),
        "evidence": evidence_type,
        "min_depth": min_depth,
        "possible_mixed_infection": mixed,
        "minor_lineages": minors,
        "ambiguous_sublineages": ambiguous,
    }
    if ambiguous:
        result["reason"] = "sublineages " + " and ".join(ambiguous) + " are both supported; reporting their common parent"
    if any(stats[name]["inferred"] for name in path):
        result["reference_inferred"] = True
    return result


def consistent_minor(name, stats, known):
    return all(stats[a]["state"] != "absent" or (stats[a]["fraction"] or 0) >= MINOR_FRACTION
               for a in ancestors_of(name, known))


def main():
    parser = argparse.ArgumentParser(description="Classify the M. tuberculosis lineage from barcode markers")
    parser.add_argument("--vcf", required=True, help="Filtered VCF of the sample")
    parser.add_argument("--pileup", help="Allele counts at barcode sites (bcftools query of mpileup AD)")
    parser.add_argument("--barcode", required=True, help="Barcode BED file")
    parser.add_argument("--reference", required=True, help="Reference FASTA")
    parser.add_argument("--min_depth", type=int, default=5, help="Minimum read depth for a barcode site")
    parser.add_argument("--sample_id", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    markers, families = load_barcode_bed(args.barcode)
    if not markers:
        result = unknown_result("failed to load the lineage barcode", "none", args.min_depth)
    else:
        wanted = {(m["chrom"], m["pos"]) for m in markers}
        if args.pileup:
            observed = read_pileup(args.pileup, wanted)
            observed = {k: v for k, v in observed.items() if v["depth"] >= args.min_depth}
            evidence = marker_evidence(markers, observed, {}, False)
            result = classify_lineage(markers, families, evidence, "pileup", args.min_depth)
        else:
            observed, snps = read_vcf(args.vcf, wanted)
            infer = snps >= GENOME_WIDE_SNPS
            reference = read_reference(args.reference) if infer else {}
            evidence = marker_evidence(markers, observed, reference, infer)
            result = classify_lineage(markers, families, evidence, "vcf", args.min_depth)
            result["vcf_snps"] = snps
    result["sample_id"] = args.sample_id
    result["analysis_date"] = datetime.now().isoformat(timespec="seconds")
    with open(args.output, "w") as handle:
        json.dump(result, handle, indent=2)


if __name__ == "__main__":
    main()
