
import os
import re
import uuid
from datetime import datetime, timezone

TBFHIR_NAMESPACE = uuid.uuid5(uuid.NAMESPACE_URL, "https://github.com/oucru-id/tb-to-fhir-full")

SOURCE_REPOSITORY = "https://github.com/oucru-id/tb-to-fhir-full"


def sanitize_id(id_string):
    sanitized = re.sub(r'[^A-Za-z0-9\-]', '-', str(id_string))
    sanitized = re.sub(r'-+', '-', sanitized)
    sanitized = sanitized.strip('-')
    return sanitized if sanitized else "unknown"


def fhir_id(*parts):
    value = "-".join(sanitize_id(p) for p in parts if p not in (None, ""))
    if len(value) > 64:
        raise ValueError(f"FHIR id longer than 64 characters: {value!r}")
    return value


def full_url(resource_type, resource_id):
    return f"urn:uuid:{uuid.uuid5(TBFHIR_NAMESPACE, f'{resource_type}/{resource_id}')}"


def walk_references(node, fn):
    if isinstance(node, dict):
        for key, value in list(node.items()):
            if key == "reference" and isinstance(value, str):
                node[key] = fn(value)
            else:
                walk_references(value, fn)
    elif isinstance(node, list):
        for item in node:
            walk_references(item, fn)


def point_references_at_entries(resources):
    index = {f"{r['resourceType']}/{r['id']}": full_url(r["resourceType"], r["id"]) for r in resources}
    for resource in resources:
        walk_references(resource, lambda value: index.get(value, value))
    return resources


def bundle_id(*parts):
    return str(uuid.uuid5(TBFHIR_NAMESPACE, "Bundle/" + "/".join(str(p) for p in parts)))


def _fixed_now():
    fixed = os.environ.get("TBFHIR_FIXED_TIMESTAMP")
    if fixed:
        value = datetime.fromisoformat(fixed.replace("Z", "+00:00"))
        return value if value.tzinfo else value.replace(tzinfo=timezone.utc)
    epoch = os.environ.get("SOURCE_DATE_EPOCH")
    if epoch:
        return datetime.fromtimestamp(int(epoch), tz=timezone.utc)
    return None


def now():
    return _fixed_now() or datetime.now(timezone.utc)


def now_iso():
    return now().isoformat()


def today_compact():
    fixed = _fixed_now()
    return (fixed or datetime.now()).strftime('%Y%m%d')


def xhtml_text(value):
    return str(value).replace("&", "&amp;").replace("<", "&lt;")


def normalise_sample_id(sample_id):
    value = str(sample_id)
    for suffix in ("_ont.fastq", "_ont", "_illumina"):
        if value.endswith(suffix):
            return value[: -len(suffix)]
    return value


def lookup_sample(records, sample_id):
    if not records:
        return None
    wanted = {str(sample_id), normalise_sample_id(sample_id)}
    for key in wanted:
        if key in records:
            return records[key]
    wanted_fid = {sanitize_id(w) for w in wanted}
    matches = [v for k, v in records.items()
               if sanitize_id(normalise_sample_id(k)) in wanted_fid]
    if len(matches) > 1:
        raise ValueError(f"sample {sample_id!r} matches more than one record after sanitising")
    return matches[0] if matches else None
