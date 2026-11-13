"""The documented differences between v1.6 and 2.0 (``--fhir_profile satusehat``).

``expected_v2(v1_bundle)`` applies each intended change to a *normalised*
v1 merged bundle. The 2.0 linked bundle, normalised the same way and with
the new Encounter removed, must then be identical. Anything not listed here
is a regression.

1.  HGVS system is ``http://varnomen.hgvs.org`` (v1: https).
2.  g.HGVS is correct for indels and MNVs (v1 wrote ``g.{pos}{ref}>{alt}`` for all).
3.  The non-IG ``tbd-codes-cs`` codes (coverage-breadth, variant-quality,
    mapping-quality) use the project CodeSystem, whose display the coding
    must carry; the depth threshold stays in ``text``.
4.  Patient names are only in the Patient resource: no ``display`` on
    DiagnosticReport.subject / .specimen, Specimen.subject, ServiceRequest.subject,
    no name in the Specimen note or the HTML report.
5.  No placeholders: v1 wrote ``"Unknown"`` for missing names and birth dates
    and the run time for missing specimen dates; 2.0 omits them.
6.  DiagnosticReport.basedOn points at the ServiceRequest.
7.  Provenance.entity names the annotated VCF (v1: the intermediate bundle).
8.  An Encounter is created (v1 referenced one that did not exist); it is
    compared separately in tests/test_clinical.py.
9.  DiagnosticReport has a narrative (``text``) summarising the conclusion.
"""
import base64  # noqa: F401  (kept for readers decoding presentedForm by hand)
import copy
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "bin"))
from tbfhir_profiles import CS_OBSERVATION_CODE  # noqa: E402
from tbfhir_vcf import genomic_hgvs  # noqa: E402

TBD_CS = "http://hl7.org/fhir/uv/genomics-reporting/CodeSystem/tbd-codes-cs"
V1_GHGVS = re.compile(r"^(NC_000962\.3):g\.(\d+)([ACGTN]+)>([ACGTN]+)$")
TIMESTAMP = "<TIMESTAMP>"


def _replace_strings(node, old, new):
    if isinstance(node, dict):
        return {k: _replace_strings(v, old, new) for k, v in node.items()}
    if isinstance(node, list):
        return [_replace_strings(v, old, new) for v in node]
    if isinstance(node, str):
        return node.replace(old, new)
    return node


COVERAGE_BREADTH_DISPLAY = "Fraction of bases at or above the minimum depth"


def _fix_project_code_displays(observation):
    for comp in observation.get("component", []):
        for coding in comp.get("code", {}).get("coding", []):
            if coding.get("system") == CS_OBSERVATION_CODE and coding.get("code") == "coverage-breadth":
                coding["display"] = COVERAGE_BREADTH_DISPLAY


def _fix_ghgvs(observation):
    for comp in observation.get("component", []):
        if comp.get("code", {}).get("coding", [{}])[0].get("code") != "81290-9":
            continue
        coding = comp["valueCodeableConcept"]["coding"][0]
        match = V1_GHGVS.match(coding["code"])
        if match:
            refseq, pos, ref, alt = match.groups()
            coding["code"] = coding["display"] = genomic_hgvs(refseq, int(pos), ref, alt)


def expected_v2(v1_bundle, vcf_name):
    bundle = copy.deepcopy(v1_bundle)
    bundle = _replace_strings(bundle, "https://varnomen.hgvs.org", "http://varnomen.hgvs.org")
    bundle = _replace_strings(bundle, TBD_CS, CS_OBSERVATION_CODE)

    resources = [e["resource"] for e in bundle["entry"]]
    has = {r["resourceType"] for r in resources}
    for r in resources:
        rtype = r["resourceType"]
        if rtype == "Observation":
            _fix_ghgvs(r)
            _fix_project_code_displays(r)
        elif rtype == "DiagnosticReport":
            sid = r["id"][: -len("-genomic-report")]
            r["text"] = {"status": "generated",
                         "div": ('<div xmlns="http://www.w3.org/1999/xhtml">TB genomic analysis report '
                                 f'for sample {sid}. {r["conclusion"]}</div>')}
            r["subject"].pop("display", None)
            for spec in r.get("specimen", []):
                spec.pop("display", None)
            for form in r.get("presentedForm", []):
                form["data"] = re.sub(r"<p><strong>Patient:</strong> .*?</p>\n", "", form["data"])
            if "ServiceRequest" in has:
                sid = r["id"][: -len("-genomic-report")]
                r["basedOn"] = [{"reference": f"ServiceRequest/{sid}-service-request"}]
        elif rtype == "Specimen":
            sid = r["id"][: -len("-specimen")]
            r["subject"].pop("display", None)
            r["note"] = [{"text": f"Collected sputum sample ({sid}) for TB genetic testing purposes"}]
            if r.get("receivedTime") == TIMESTAMP:
                del r["receivedTime"]
            if r.get("collection", {}).get("collectedDateTime") == TIMESTAMP:
                del r["collection"]["collectedDateTime"]
        elif rtype == "ServiceRequest":
            r["subject"].pop("display", None)
        elif rtype == "Patient":
            names = []
            for name in r.get("name", []):
                name = {k: v for k, v in name.items() if v != "Unknown"}
                if "given" in name:
                    name["given"] = [g for g in name["given"] if g != "Unknown"] or None
                    if name["given"] is None:
                        del name["given"]
                if set(name) - {"use"}:
                    names.append(name)
            if names:
                r["name"] = names
            else:
                r.pop("name", None)
            if r.get("birthDate") == "Unknown":
                del r["birthDate"]
        elif rtype == "Provenance":
            for entity in r.get("entity", []):
                what = entity.get("what", {})
                if "display" in what:
                    what["display"] = re.sub(r"\(([^()]*)\.fhir\.json\)$", f"({vcf_name})", what["display"])
    return bundle


def without_encounters(bundle):
    bundle = copy.deepcopy(bundle)
    bundle["entry"] = [e for e in bundle["entry"] if e["resource"]["resourceType"] != "Encounter"]
    return bundle
