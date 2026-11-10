
from tbfhir_resistance import snomed_conclusion_coding
from tbfhir_vcf import WHO_GRADE_GROUP

LOINC = "http://loinc.org"
CANONICAL_BASE = "http://terminology.spheres.id"
CS_WHO_GRADE = f"{CANONICAL_BASE}/CodeSystem/who-catalogue-grade"
CS_OBSERVATION_CODE = f"{CANONICAL_BASE}/CodeSystem/tbfhir-observation-code"
CS_LINEAGE = f"{CANONICAL_BASE}/CodeSystem/mtb-lineage"
CS_LINEAGE_ATTRIBUTE = f"{CANONICAL_BASE}/CodeSystem/mtb-lineage-attribute"
CS_PIPELINE_PROPERTY = f"{CANONICAL_BASE}/CodeSystem/pipeline-property"

GENOMICS_SD = "http://hl7.org/fhir/uv/genomics-reporting/StructureDefinition"
FHIR_OBSERVATION_SD = "http://hl7.org/fhir/StructureDefinition/Observation"

KEMKES_SD = "https://fhir.kemkes.go.id/r4/StructureDefinition"
KEMKES_TAG = {"system": "http://terminology.kemkes.go.id/sp", "code": "genomics",
              "display": "Genomics"}
KEMKES_SYS_IDS = "http://sys-ids.kemkes.go.id"
KEMKES_NIK = "https://fhir.kemkes.go.id/id/nik"

WHO_GRADE_DISPLAY = {
    '1': 'Assoc w R', '2': 'Assoc w R - Interim', '3': 'Uncertain significance',
    '4': 'Not assoc w R - Interim', '5': 'Not assoc w R',
}

_SATUSEHAT_WHO_GRADE = {
    'Uncertain significance': {'system': LOINC, 'code': 'LA26333-7',
                               'display': 'Uncertain significance'},
    'Assoc w R': {'system': 'http://terminology.kemkes.go.id/sp', 'code': 'SP000478',
                  'display': 'Assoc w R'},
    'Assoc w R - Interim': {'system': 'http://terminology.kemkes.go.id/sp', 'code': 'SP000479',
                            'display': 'Assoc w R - Interim'},
    'Not assoc w R - Interim': {'system': 'http://terminology.kemkes.go.id/sp', 'code': 'SP000480',
                                'display': 'Not assoc w R - Interim'},
    'Not assoc w R': {'system': 'http://terminology.kemkes.go.id/sp', 'code': 'SP000481',
                      'display': 'Not assoc w R'},
}

_SATUSEHAT_CONCLUSION = {
    "No resistance detected": {
        "system": "https://terminology.kemkes.go.id/CodeSystem/episodeofcare-type",
        "code": "TB-SO", "display": "Tuberkulosis Sensitif Obat"},
    "Pre-XDR-TB": {
        "system": "http://terminology.kemkes.go.id/CodeSystem/clinical-term",
        "code": "OV000435", "display": "Pre-XDR"},
}

REPORT_CODE_V1 = {"system": LOINC, "code": "81247-9",
                  "display": "Master HL7 genetic variant reporting panel"}
REPORT_CODE_IG3 = {"system": LOINC, "code": "51969-4", "display": "Genetic analysis report"}


class Profile:

    name = "generic"
    include_nik = False
    lineage_legacy_coding = False
    kemkes_extensions = False

    def __init__(self, identifier_system_base=None):
        self.identifier_system_base = (identifier_system_base or "").rstrip("/") or None

    def meta(self, kind):
        profiles = self._profiles(kind)
        meta = {}
        if profiles:
            meta["profile"] = profiles
        tags = self._tags(kind)
        if tags:
            meta["tag"] = tags
        return meta or None

    def _profiles(self, kind):
        return {
            "variant": [f"{GENOMICS_SD}/variant"],
            "panel": [FHIR_OBSERVATION_SD],
            "lineage": [FHIR_OBSERVATION_SD],
            "report": [f"{GENOMICS_SD}/genomic-report"],
        }.get(kind, [])

    def _tags(self, kind):
        return []

    def bundle_meta(self, last_updated):
        return {"lastUpdated": last_updated}

    def who_grade_coding(self, grade):
        group = WHO_GRADE_GROUP.get(grade.strip().lower())
        if not group:
            return None
        return {"system": CS_WHO_GRADE, "code": group, "display": WHO_GRADE_DISPLAY[group]}

    def conclusion_coding(self, resistance_class):
        return snomed_conclusion_coding(resistance_class)

    def report_code(self):
        return {"coding": [dict(REPORT_CODE_IG3), dict(REPORT_CODE_V1)],
                "text": "TB Genomic Analysis Report"}

    def protein_change_component(self, clean_variant_id, variant_id):
        return {"code": {"coding": [{"system": CS_OBSERVATION_CODE, "code": "who-variant-id",
                                     "display": "WHO catalogue variant identifier"}],
                         "text": "WHO catalogue variant identifier"},
                "valueString": variant_id}

    def identifier(self, kind, value, org_id=None):
        """Business identifier for a resource kind, or None."""
        if not self.identifier_system_base or value in (None, ""):
            return None
        return {"system": f"{self.identifier_system_base}/{kind}", "value": str(value)}


class SatusehatProfile(Profile):
    name = "satusehat"
    include_nik = True
    lineage_legacy_coding = True
    kemkes_extensions = True

    _KEMKES = {
        "patient": "Patient", "specimen": "Specimen", "organization": "Organization",
        "practitioner": "Practitioner", "practitioner_role": "PractitionerRole",
        "service_request": "ServiceRequest", "encounter": "Encounter",
    }

    def _profiles(self, kind):
        if kind in self._KEMKES:
            return [f"{KEMKES_SD}/{self._KEMKES[kind]}"]
        if kind == "region":
            return [f"{GENOMICS_SD}/region-studied"]
        if kind == "report":
            return [f"{GENOMICS_SD}/genomics-report"]
        return super()._profiles(kind)

    def _tags(self, kind):
        if kind in ("variant", "region", "lineage", "report", "device", "provenance"):
            return [dict(KEMKES_TAG)]
        return []

    def bundle_meta(self, last_updated):
        return {"lastUpdated": last_updated, "profile": [f"{KEMKES_SD}/Bundle"]}

    def who_grade_coding(self, grade):
        coding = _SATUSEHAT_WHO_GRADE.get(grade.strip())
        return dict(coding) if coding else None

    def conclusion_coding(self, resistance_class):
        if resistance_class in _SATUSEHAT_CONCLUSION:
            return dict(_SATUSEHAT_CONCLUSION[resistance_class])
        return snomed_conclusion_coding(resistance_class)

    def report_code(self):
        return {"coding": [dict(REPORT_CODE_V1)], "text": "TB Genomic Analysis Report"}

    def protein_change_component(self, clean_variant_id, variant_id):
        return {"code": {"coding": [{"system": LOINC, "code": "48005-3",
                                     "display": "Amino acid change (pHGVS)"}]},
                "valueCodeableConcept": {
                    "coding": [{"system": "http://varnomen.hgvs.org",
                                "code": f"NC_000962.3:p.({clean_variant_id})",
                                "display": clean_variant_id}],
                    "text": variant_id}}

    def identifier(self, kind, value, org_id=None):
        if value in (None, ""):
            return None
        if kind == "organization":
            return {"use": "official", "system": f"{KEMKES_SYS_IDS}/organization", "value": str(value)}
        if not org_id:
            return None
        system_kind = {"service_request": "servicerequest"}.get(kind, kind)
        return {"system": f"{KEMKES_SYS_IDS}/{system_kind}/{org_id}", "value": str(value)}


PROFILES = {"generic": Profile, "satusehat": SatusehatProfile}


def get_profile(name, identifier_system_base=None):
    try:
        return PROFILES[name](identifier_system_base=identifier_system_base)
    except KeyError:
        raise SystemExit(f"unknown --fhir_profile {name!r}; choose one of {sorted(PROFILES)}")
