"""Table-driven tests for classify_drug_resistance (every classification branch)."""
import pytest

from tbfhir_resistance import classify_drug_resistance, extract_lineage_info
from tbfhir_terminology import get_drug_panel_config

R, S, I = "LA6676-6", "LA24225-7", "LA9663-1"
PANEL = {key: (code, display) for key, code, display in get_drug_panel_config()}


def panel(states, default=S):
    components = []
    for key, (code, display) in PANEL.items():
        value = states.get(key, default)
        components.append({"code": {"coding": [{"code": code, "display": display}]},
                           "valueCodeableConcept": {"coding": [{"code": value}]}})
    return {"code": {"coding": [{"code": "89486-5"}]}, "component": components}


def variant(gene, significance):
    return {"code": {"coding": [{"code": "69548-6"}]}, "component": [
        {"code": {"coding": [{"display": "Gene studied [ID]"}]}, "valueCodeableConcept": {"text": gene}},
        {"code": {"coding": [{"display": "Genetic variation clinical significance [Imp]"}]},
         "valueCodeableConcept": {"text": significance}}]}


CASES = [
    ("all susceptible", {}, S, "No resistance detected"),
    ("some indeterminate", {"pretomanid": I}, S, "No resistance detected - partial"),
    ("nothing assessable", {}, I, "Indeterminate"),
    ("rifampicin only", {"rifampicin": R}, S, "RR-TB"),
    ("rif + inh", {"rifampicin": R, "isoniazid": R}, S, "MDR-TB"),
    ("mdr + fq", {"rifampicin": R, "isoniazid": R, "levofloxacin": R}, S, "Pre-XDR-TB"),
    ("rr + fq", {"rifampicin": R, "moxifloxacin": R}, S, "Pre-XDR-TB"),
    ("mdr + fq + bdq", {"rifampicin": R, "isoniazid": R, "levofloxacin": R, "bedaquiline": R}, S, "XDR-TB"),
    ("mdr + fq + lzd", {"rifampicin": R, "isoniazid": R, "moxifloxacin": R, "linezolid": R}, S, "XDR-TB"),
    ("inh, rif susceptible", {"isoniazid": R}, S, "HR-TB"),
    ("inh, rif indeterminate", {"isoniazid": R, "rifampicin": I}, S,
     "Drug-resistant - rifampicin not assessable"),
    ("streptomycin mono", {"streptomycin": R}, S, "Streptomycin-resistant TB"),
    ("ethionamide mono", {"ethionamide": R}, S, "Ethionamide-resistant TB"),
    ("pyrazinamide mono", {"pyrazinamide": R}, S, "Pyrazinamide-resistant TB"),
    ("ethambutol mono", {"ethambutol": R}, S, "Ethambutol-resistant TB"),
    ("fq mono", {"levofloxacin": R, "moxifloxacin": R}, S, "Drug-resistant"),
    ("two non-rif groups", {"ethambutol": R, "pyrazinamide": R}, S, "Drug-resistant"),
]


@pytest.mark.parametrize("label,states,default,expected", CASES, ids=[c[0] for c in CASES])
def test_classification(label, states, default, expected):
    classification, *_ = classify_drug_resistance([panel(states, default)])
    assert classification == expected


def test_rr_with_isoniazid_indeterminate_says_mdr_not_excluded():
    _cls, description, *_ = classify_drug_resistance([panel({"rifampicin": R, "isoniazid": I})])
    assert "MDR-TB cannot be excluded" in description


def test_resistance_genes_need_an_assoc_w_r_grade():
    obs = [panel({"rifampicin": R}), variant("rpoB", "Assoc w R (Rifampicin)"),
           variant("embB", "Uncertain significance (Ethambutol)")]
    _cls, _desc, genes, drugs, _ind = classify_drug_resistance(obs)
    assert genes == ["rpoB"]
    assert drugs == ["rifampin"]


def test_indeterminate_drugs_are_listed():
    *_rest, indeterminate = classify_drug_resistance([panel({"pretomanid": I, "cycloserine": I})])
    assert indeterminate == ["cycloserine", "pretomanid"]


def test_lineage_text_is_extracted():
    lineage = {"code": {"coding": [{"code": "614-8"}]}, "valueCodeableConcept": {"text": "Lineage lineage4"}}
    assert extract_lineage_info([panel({}), lineage]) == "Lineage lineage4"
    assert extract_lineage_info([panel({})]) is None
