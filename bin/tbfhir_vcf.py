
import gzip


class AnnotationError(Exception):
    pass


class VariantRecord:
    __slots__ = ("CHROM", "POS", "REF", "ALT", "QUAL", "INFO")

    def __init__(self, chrom, pos, ref, alt, qual, info):
        self.CHROM = chrom
        self.POS = int(pos)
        self.REF = ref
        self.ALT = [alt]
        self.QUAL = qual
        self.INFO = info


def _open_text(path):
    with open(path, "rb") as probe:
        magic = probe.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def read_vcf(path):
    variants = []
    with _open_text(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            chrom, pos, _id, ref, alt, qual_str, _filter, info_str = fields[:8]
            try:
                qual = float(qual_str) if qual_str not in (".", "") else None
            except ValueError:
                qual = None
            info = {}
            if info_str != ".":
                for item in info_str.split(";"):
                    if "=" in item:
                        key, value = item.split("=", 1)
                        info[key] = value
                    elif item:
                        info[item] = True
            variants.append(VariantRecord(chrom, pos, ref, alt, qual, info))
    return variants


def get_variant_value(record, field, default='unknown'):
    value = record.INFO.get(field)
    return value if value is not None else default


def get_variant_list(record, field):
    raw = get_variant_value(record, field, default=None)
    if raw is None or isinstance(raw, bool):
        return []
    values = [v.strip() for v in str(raw).split(',')]
    if not values or all(not v or v.lower() == 'unknown' for v in values):
        return []
    return values


WHO_GRADE_SEVERITY = {
    'assoc w r': 5,
    'assoc w r - interim': 4,
    'uncertain significance': 3,
    'not assoc w r - interim': 2,
    'not assoc w r': 1,
}

WHO_GRADE_GROUP = {
    'assoc w r': '1',
    'assoc w r - interim': '2',
    'uncertain significance': '3',
    'not assoc w r - interim': '4',
    'not assoc w r': '5',
}


def most_severe_grade(grades):
    if not grades:
        return 'unknown'
    return max(grades, key=lambda g: WHO_GRADE_SEVERITY.get(g.strip().lower(), 0))


def is_resistance_grade(grade):
    return grade.strip().lower().startswith('assoc w r')


def resolve_variant_drug_pairs(record, allow_legacy=False):

    drugs = get_variant_list(record, 'DRUG')
    grades = get_variant_list(record, 'WHO_CLASSIFICATION')
    variant_ids = get_variant_list(record, 'VARIANT_ID')

    if not drugs:
        return []

    if len(grades) == len(drugs):
        pairs = []
        for idx, drug in enumerate(drugs):
            variant_id = variant_ids[idx] if idx < len(variant_ids) else (
                variant_ids[0] if len(variant_ids) == 1 else None)
            pairs.append((drug, grades[idx], variant_id))
        return pairs

    location = f"{record.CHROM}:{record.POS} {record.REF}>{record.ALT[0]}"

    if not allow_legacy:
        raise AnnotationError(
            f"{location}: {len(drugs)} drugs but {len(grades)} WHO grades "
            f"({drugs} / {grades}).\n"
            f"  The annotation table has one row per (variant, drug) pair, so "
            f"bcftools annotate collapsed the grade and the per-drug association is "
            f"no longer recoverable.\n"
            f"  Fix: rebuild the table with tools/annotation/rebuild_annotation_table.sh, then "
            f"point params.annotation_table / params.annotation_header at the result.\n"
            f"  To proceed anyway with the most severe grade applied to every listed "
            f"drug (this over-calls resistance), rerun with --legacy-annotation-table."
        )

    grade = most_severe_grade(grades) if grades else 'unknown'
    variant_id = variant_ids[0] if variant_ids else None
    return [(drug, grade, variant_id) for drug in drugs]


def compute_allele_fraction(record):
    dp4 = get_variant_value(record, 'DP4', default=None)
    if dp4 and dp4 != 'unknown':
        try:
            counts = [int(x) for x in str(dp4).split(',')]
            if len(counts) == 4 and sum(counts) > 0:
                return (counts[2] + counts[3]) / sum(counts)
        except (ValueError, ZeroDivisionError):
            pass

    for field in ('AF', 'AF1'):
        value = get_variant_value(record, field, default=None)
        if value and value != 'unknown':
            try:
                fraction = float(str(value).split(',')[0])
                if 0.0 <= fraction <= 1.0:
                    return fraction
            except ValueError:
                continue
    return None


def genomic_hgvs(refseq, pos, ref, alt):

    if len(ref) == 1 and len(alt) == 1:
        return f"{refseq}:g.{pos}{ref}>{alt}"
    if len(alt) == 1 and len(ref) > 1 and ref[0] == alt[0]:
        start, end = pos + 1, pos + len(ref) - 1
        return f"{refseq}:g.{start}del" if start == end else f"{refseq}:g.{start}_{end}del"
    if len(ref) == 1 and len(alt) > 1 and ref[0] == alt[0]:
        return f"{refseq}:g.{pos}_{pos + 1}ins{alt[1:]}"
    end = pos + len(ref) - 1
    location = f"{pos}" if end == pos else f"{pos}_{end}"
    return f"{refseq}:g.{location}delins{alt}"
