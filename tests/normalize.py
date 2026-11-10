"""Normalise FHIR bundles so runs can be compared byte for byte.

Volatile parts are replaced with fixed placeholders:

* ``Bundle.id``, ``Bundle.timestamp`` and every ISO-8601 date-time string
* entry ``fullUrl`` values become ``urn:tbfhir:<Type>/<id>``, and references
  that point at a fullUrl are rewritten to ``<Type>/<id>``. Bundles that
  reference by ``Type/id`` and bundles that reference by ``urn:uuid`` then
  normalise to the same thing.
* the date suffix of ``TB-GEN-<sample>-YYYYMMDD`` identifiers
* ``presentedForm.data`` is base64-decoded (and its timestamps normalised)
  so differences in the HTML report are readable in a diff.
"""
import base64
import copy
import json
import re

ISO_DATETIME = re.compile(
    r"\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}(?:\.\d+)?(?:Z|[+-]\d{2}:\d{2})?")
DATED_IDENTIFIER = re.compile(r"^(TB-GEN-.+)-\d{8}$")
TIMESTAMP = "<TIMESTAMP>"


def _replace_datetimes(text):
    return ISO_DATETIME.sub(TIMESTAMP, text)


def _walk(node, fullurl_map):
    if isinstance(node, dict):
        out = {}
        for key, value in node.items():
            if key == "reference" and isinstance(value, str) and value in fullurl_map:
                out[key] = fullurl_map[value]
            elif key == "data" and isinstance(value, str) and node.get("contentType", "").startswith("text/"):
                try:
                    decoded = base64.b64decode(value).decode("utf-8")
                    out[key] = "base64-decoded:" + _replace_datetimes(decoded)
                except (ValueError, UnicodeDecodeError):
                    out[key] = value
            elif key == "value" and isinstance(value, str) and DATED_IDENTIFIER.match(value):
                out[key] = DATED_IDENTIFIER.sub(r"\1-<DATE>", value)
            else:
                out[key] = _walk(value, fullurl_map)
        return out
    if isinstance(node, list):
        return [_walk(item, fullurl_map) for item in node]
    if isinstance(node, str):
        return _replace_datetimes(node)
    return node


def normalize_bundle(bundle, sort_entries=False):
    bundle = copy.deepcopy(bundle)
    fullurl_map = {}
    for entry in bundle.get("entry", []):
        resource = entry.get("resource", {})
        rtype, rid = resource.get("resourceType"), resource.get("id")
        if entry.get("fullUrl") and rtype and rid:
            fullurl_map[entry["fullUrl"]] = f"{rtype}/{rid}"
    for entry in bundle.get("entry", []):
        resource = entry.get("resource", {})
        rtype, rid = resource.get("resourceType"), resource.get("id")
        if "fullUrl" in entry:
            entry["fullUrl"] = f"urn:tbfhir:{rtype}/{rid}" if rtype and rid else "urn:tbfhir:anonymous"
    if "id" in bundle:
        bundle["id"] = "<BUNDLE-ID>"
    bundle = _walk(bundle, fullurl_map)
    if sort_entries:
        bundle["entry"] = sorted(
            bundle.get("entry", []),
            key=lambda e: (e.get("resource", {}).get("resourceType", ""),
                           e.get("resource", {}).get("id", "")))
    return bundle


def dumps(bundle):
    return json.dumps(bundle, indent=2, sort_keys=True, ensure_ascii=False) + "\n"


def load_normalized(path, sort_entries=False):
    with open(path, encoding="utf-8") as handle:
        return normalize_bundle(json.load(handle), sort_entries=sort_entries)
