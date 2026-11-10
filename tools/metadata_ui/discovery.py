"""Discover sample IDs from the pipeline input directories.

Mirrors the channel definitions in main.nf so the list shown here matches what
Nextflow will actually process.
"""

import glob
import os


def _strip(name, suffixes):
    for suffix in suffixes:
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return name


def discover_sample_ids(base_dir):
    """Return {sample_id: platform} for every input the pipeline would pick up.

    Nanopore IDs are stored WITHOUT the _ont suffix. main.nf carries _ont through
    to the FHIR sample id, but find_matching_sample() in
    scripts/clinical_metadata_parser.py normalises it away, so the bare accession
    is the form that matches from either direction.
    """
    ngs_dir = os.path.join(base_dir, "data", "NGS")
    vcf_dir = os.path.join(base_dir, "data", "VCF")
    found = {}

    # Illumina: *_1_illumina.fastq.gz / *_2_illumina.fastq.gz -> one id per pair
    for path in glob.glob(os.path.join(ngs_dir, "*_illumina.fastq.gz")):
        name = os.path.basename(path)
        sample = _strip(name, ["_1_illumina.fastq.gz", "_2_illumina.fastq.gz"])
        if sample != name:
            found[sample] = "illumina"

    # Nanopore: *_ont.fastq.gz
    for path in glob.glob(os.path.join(ngs_dir, "*_ont.fastq.gz")):
        name = os.path.basename(path)
        sample = _strip(name, ["_ont.fastq.gz"])
        if sample != name:
            found.setdefault(sample, "nanopore")

    # Pre-called VCFs
    for pattern in ("*.vcf", "*.vcf.gz"):
        for path in glob.glob(os.path.join(vcf_dir, pattern)):
            name = os.path.basename(path)
            sample = _strip(name, [".vcf.gz", ".vcf"])
            if sample != name:
                found.setdefault(sample, "vcf")

    return found


def group_samples(discovered, existing_ids):
    """Split discovered vs stored ids into the three states the UI shows."""
    discovered_ids = set(discovered)
    existing = set(existing_ids)
    return {
        "needs_metadata": sorted(discovered_ids - existing),
        "complete": sorted(discovered_ids & existing),
        "orphaned": sorted(existing - discovered_ids),
    }
