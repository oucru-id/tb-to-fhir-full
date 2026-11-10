
import csv
import glob
import os


def _strip(name, suffixes):
    for suffix in suffixes:
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return name


def from_samplesheet(path):
    found = {}
    with open(path, newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            sample = (row.get("sample") or "").strip()
            if sample:
                found[sample] = (row.get("platform") or "").strip() or "unknown"
    return found


def from_directories(reads_dir=None, vcf_dir=None):
    found = {}
    if reads_dir:
        for path in glob.glob(os.path.join(reads_dir, "*_illumina.fastq.gz")):
            name = os.path.basename(path)
            sample = _strip(name, ["_1_illumina.fastq.gz", "_2_illumina.fastq.gz"])
            if sample != name:
                found[sample] = "illumina"
        for path in glob.glob(os.path.join(reads_dir, "*_ont.fastq.gz")):
            name = os.path.basename(path)
            sample = _strip(name, ["_ont.fastq.gz"])
            if sample != name:
                found.setdefault(sample, "ont")
    if vcf_dir:
        for pattern in ("*.vcf", "*.vcf.gz"):
            for path in glob.glob(os.path.join(vcf_dir, pattern)):
                name = os.path.basename(path)
                sample = _strip(name, [".vcf.gz", ".vcf"])
                if sample != name:
                    found.setdefault(sample, "vcf")
    return found


def discover_sample_ids(samplesheet=None, reads_dir=None, vcf_dir=None):
    found = {}
    if samplesheet and os.path.exists(samplesheet):
        found.update(from_samplesheet(samplesheet))
    for sample, platform in from_directories(reads_dir, vcf_dir).items():
        found.setdefault(sample, platform)
    return found


def group_samples(discovered, existing_ids):
    discovered_ids = set(discovered)
    existing = set(existing_ids)
    return {
        "needs_metadata": sorted(discovered_ids - existing),
        "complete": sorted(discovered_ids & existing),
        "orphaned": sorted(existing - discovered_ids),
    }
