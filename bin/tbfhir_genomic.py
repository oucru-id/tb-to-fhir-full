#!/usr/bin/env python3

import argparse
import base64
import json
import os
import sys

from tbfhir_common import (SOURCE_REPOSITORY, bundle_id, fhir_id, full_url, lookup_sample,
                           now_iso, point_references_at_entries, sanitize_id, today_compact,
                           xhtml_text)
from tbfhir_profiles import (CS_LINEAGE, CS_LINEAGE_ATTRIBUTE, CS_OBSERVATION_CODE,
                             CS_PIPELINE_PROPERTY, get_profile)
from tbfhir_resistance import classify_drug_resistance, extract_lineage_info
from tbfhir_terminology import (HGVS, LOINC, REFERENCE_SEQUENCE, REFSEQ, SEQUENCE_ONTOLOGY,
                                SNOMED, GeneMap, get_drug_panel_config,
                                get_drug_snomed_mapping, get_effect_so_mapping)
from tbfhir_vcf import (AnnotationError, compute_allele_fraction, genomic_hgvs,
                        get_variant_list, get_variant_value, is_resistance_grade, read_vcf,
                        resolve_variant_drug_pairs)

CATEGORY_LAB = {"coding": [{"system": "http://terminology.hl7.org/CodeSystem/observation-category",
                            "code": "laboratory", "display": "Laboratory"}]}
CATEGORY_GE = {"coding": [{"system": "http://terminology.hl7.org/CodeSystem/v2-0074",
                           "code": "GE", "display": "Genetics"}]}
REFSEQ_COMPONENT = {
    "code": {"coding": [{"system": LOINC, "code": "48013-7",
                         "display": "Genomic reference sequence ID"}]},
    "valueCodeableConcept": {"coding": [{"system": REFSEQ, "code": REFERENCE_SEQUENCE,
                                         "display": REFERENCE_SEQUENCE}],
                             "text": REFERENCE_SEQUENCE},
}


def div(text):
    return f"<div xmlns=\"http://www.w3.org/1999/xhtml\">{xhtml_text(text)}</div>"


def with_meta(resource, meta):
    if not meta:
        return resource
    out = {}
    for key, value in resource.items():
        out[key] = value
        if key == "id":
            out["meta"] = meta
    return out


class SampleContext:
    def __init__(self, sample_id, fid, profile, performer, timestamp):
        self.sample_id = sample_id
        self.fid = fid
        self.profile = profile
        self.performer = performer            
        self.timestamp = timestamp

    def subject_refs(self):
        refs = {"subject": {"reference": f"Patient/{self.fid}-patient"},
                "specimen": {"reference": f"Specimen/{self.fid}-specimen"},
                "effectiveDateTime": self.timestamp}
        if self.performer:
            refs["performer"] = [{"reference": self.performer["reference"]}]
        return refs

TBFHIR_CODE_DISPLAY = {
    "coverage-breadth": "Fraction of bases at or above the minimum depth",
    "variant-quality": "Variant call quality (QUAL)",
    "mapping-quality": "Mapping quality (MQ)",
}


def tbfhir_code(code, text=None):
    display = TBFHIR_CODE_DISPLAY[code]
    return {"coding": [{"system": CS_OBSERVATION_CODE, "code": code, "display": display}],
            "text": text or display}


def build_variant_observation(ctx, idx, record, drug_pairs, gene_map):
    pos = str(record.POS)
    ref = record.REF
    alt = str(record.ALT[0])
    gene = get_variant_value(record, 'GENE')
    effect = get_variant_value(record, 'EFFECT')
    depth_val = get_variant_value(record, 'DP', None)
    variant_ids = [v for _d, _g, v in drug_pairs if v]
    hgvs_g = genomic_hgvs(REFERENCE_SEQUENCE, record.POS, ref, alt)

    components = [
        {"code": {"coding": [{"system": LOINC, "code": "81290-9",
                              "display": "Genomic DNA change (gHGVS)"}]},
         "valueCodeableConcept": {"coding": [{"system": HGVS, "code": hgvs_g, "display": hgvs_g}]}},
        REFSEQ_COMPONENT,
    ]

    if depth_val is not None and depth_val != 'unknown':
        try:
            components.append({
                "code": {"coding": [{"system": LOINC, "code": "82121-5", "display": "Allelic read depth"}]},
                "valueQuantity": {"value": int(depth_val), "unit": "reads per base pair",
                                  "system": "http://unitsofmeasure.org", "code": "[1]"}})
        except ValueError:
            pass

    allele_fraction = compute_allele_fraction(record)
    if allele_fraction is not None:
        components.append({
            "code": {"coding": [{"system": LOINC, "code": "81258-6",
                                 "display": "Sample variant allelic frequency [NFr]"}]},
            "valueQuantity": {"value": round(allele_fraction, 4),
                              "system": "http://unitsofmeasure.org", "code": "1"}})

    if record.QUAL is not None:
        components.append({
            "code": tbfhir_code("variant-quality", "Variant call quality (QUAL)"),
            "valueQuantity": {"value": round(record.QUAL, 2),
                              "system": "http://unitsofmeasure.org", "code": "1"}})

    mapping_quality = get_variant_value(record, 'MQ', default=None)
    if mapping_quality is not None and mapping_quality != 'unknown':
        try:
            components.append({
                "code": tbfhir_code("mapping-quality", "Mapping quality (MQ)"),
                "valueQuantity": {"value": round(float(mapping_quality), 2),
                                  "system": "http://unitsofmeasure.org", "code": "1"}})
        except ValueError:
            pass

    if gene != "unknown":
        components.append({
            "code": {"coding": [{"system": LOINC, "code": "48018-6", "display": "Gene studied [ID]"}]},
            "valueCodeableConcept": gene_map.concept(gene)})

    if effect != "unknown":
        so = get_effect_so_mapping(effect)
        value = {"text": effect}
        if so:
            value = {"coding": [{"system": SEQUENCE_ONTOLOGY, "code": so['code'],
                                 "display": so['display']}], "text": effect}
        components.append({
            "code": {"coding": [{"system": LOINC, "code": "48019-4", "display": "DNA change type"}]},
            "valueCodeableConcept": value})

    for drug, grade, _variant_id in drug_pairs:
        snomed = get_drug_snomed_mapping(drug)
        drug_value = ({"coding": [{"system": SNOMED, "code": snomed['code'],
                                   "display": snomed['display']}], "text": drug}
                      if snomed else {"text": drug})
        components.append({
            "code": {"coding": [{"system": LOINC, "code": "51963-7",
                                 "display": "Medication assessed [Identifier]"}]},
            "valueCodeableConcept": drug_value})

        who_coding = ctx.profile.who_grade_coding(grade)
        significance = ({"coding": [who_coding], "text": f"{grade} ({drug})"}
                        if who_coding else {"text": f"{grade} ({drug})"})
        components.append({
            "code": {"coding": [{"system": LOINC, "code": "53037-8",
                                 "display": "Genetic variation clinical significance [Imp]"}]},
            "valueCodeableConcept": significance})

    for variant_id in dict.fromkeys(variant_ids):
        clean = variant_id.split('_')[-1] if '_' in variant_id else variant_id
        if clean.startswith('p.'):
            clean = clean[2:]
        components.append(ctx.profile.protein_change_component(clean, variant_id))

    components.append({
        "code": {"coding": [{"system": LOINC, "code": "81254-5", "display": "Variant exact start-end"}]},
        "valueRange": {"low": {"value": record.POS},
                       "high": {"value": record.POS + max(len(ref), 1) - 1}}})

    text = f"Genomic variant at position {pos}: {ref}>{alt}"
    if gene != "unknown":
        text += f" in gene {gene}"
    if effect != "unknown":
        text += f" ({effect})"
    if variant_ids:
        text += f" - {', '.join(dict.fromkeys(variant_ids))}"
    if drug_pairs:
        text += " - WHO classification per drug - " + "; ".join(f"{d}: {g}" for d, g, _v in drug_pairs)

    observation = {
        "resourceType": "Observation",
        "id": f"{ctx.fid}-obs-{idx + 1}",
        "text": {"status": "generated", "div": div(text)},
        "status": "final",
        "category": [CATEGORY_LAB, CATEGORY_GE],
        "code": {"coding": [{"system": LOINC, "code": "69548-6", "display": "Genetic variant assessment"}]},
        "valueCodeableConcept": {"coding": [{"system": LOINC, "code": "LA9633-4", "display": "Present"}],
                                 "text": "Present"},
        **ctx.subject_refs(),
        "component": components,
    }
    return with_meta(observation, ctx.profile.meta("variant"))


def build_susceptibility_panel(ctx, resistant_drugs, coverage_info):
    coverage_assessed = bool(coverage_info and coverage_info.get('coverage_assessed'))
    coverage_drugs = {k.lower(): v for k, v in (coverage_info or {}).get('drugs', {}).items()}

    components = []
    counts = {'resistant': 0, 'no_resistance_detected': 0, 'indeterminate': 0}
    for drug_key, code, display in get_drug_panel_config():
        if drug_key in resistant_drugs:
            value_code, value_display, state = "LA6676-6", "Resistant", 'resistant'
        else:
            drug_coverage = coverage_drugs.get(drug_key)
            if coverage_assessed and drug_coverage and drug_coverage.get('assessable'):
                value_code, value_display, state = "LA24225-7", "Susceptible", 'no_resistance_detected'
            else:
                value_code, value_display, state = "LA9663-1", "Indeterminate", 'indeterminate'
        counts[state] += 1

        component = {
            "code": {"coding": [{"system": LOINC, "code": code, "display": display}]},
            "valueCodeableConcept": {"coding": [{"system": LOINC, "code": value_code,
                                                 "display": value_display}]},
        }
        if state == 'indeterminate':
            drug_coverage = coverage_drugs.get(drug_key) or {}
            if not coverage_assessed:
                reason = (coverage_info or {}).get('reason', 'no coverage data available for this sample')
            elif drug_coverage.get('inadequate_loci'):
                reason = "insufficient coverage at: " + ", ".join(drug_coverage['inadequate_loci'])
            elif drug_coverage.get('loci_without_coverage_data'):
                reason = "no coverage data for: " + ", ".join(drug_coverage['loci_without_coverage_data'])
            elif not drug_coverage:
                reason = ('the WHO mutation catalogue contains no graded variants for this '
                          'drug, so genotypic prediction is not possible; phenotypic DST is '
                          'required')
            else:
                reason = 'resistance loci for this drug are not represented in the coverage targets'
            component["valueCodeableConcept"]["text"] = (
                f"Indeterminate - no resistance-associated mutation detected, but "
                f"susceptibility cannot be confirmed ({reason})")
        components.append(component)

    summary = (f"{counts['resistant']} resistant, "
               f"{counts['no_resistance_detected']} no resistance detected, "
               f"{counts['indeterminate']} indeterminate")
    observation = {
        "resourceType": "Observation",
        "id": f"{ctx.fid}-susceptibility-panel",
        "text": {"status": "generated", "div": div(
            f"Mycobacterial susceptibility panel for {ctx.fid}: {summary}. Indeterminate means no "
            f"resistance-associated mutation was detected but coverage of the relevant loci could "
            f"not be confirmed; it is not a susceptible result.")},
        "status": "final",
        "category": [CATEGORY_LAB, CATEGORY_GE],
        "code": {"coding": [{"system": LOINC, "code": "89486-5",
                             "display": "Mycobacterial susceptibility panel Qualitative by Genotype method"}]},
        **ctx.subject_refs(),
        "component": components,
    }
    return with_meta(observation, ctx.profile.meta("panel")), summary


def build_region_observations(ctx, coverage_info, gene_map):
    if not coverage_info or not coverage_info.get('coverage_assessed'):
        return []
    panel_drugs = {key for key, _c, _d in get_drug_panel_config()}
    genes = coverage_info.get('genes', {})
    min_depth = coverage_info.get('min_depth', 10)
    observations = []
    for gene in sorted(genes):
        info = genes[gene]
        gene_drugs = [d for d in info.get('drugs', []) if d.lower() in panel_drugs]
        if not gene_drugs:
            continue
        breadth = info.get(f'breadth_at_{min_depth}x', 0.0)
        status = "adequately covered" if info.get('adequate', False) else "NOT adequately covered"
        observation = {
            "resourceType": "Observation",
            "id": f"{ctx.fid}-region-{sanitize_id(gene)}",
            "text": {"status": "generated", "div": div(
                f"Region studied: {gene} ({', '.join(sorted(gene_drugs))}) - {status}: "
                f"{float(breadth) * 100:.1f}% of bases at >= {min_depth}x, "
                f"mean depth {info.get('mean_depth', 0)}x")},
            "status": "final",
            "category": [CATEGORY_LAB, CATEGORY_GE],
            "code": {"coding": [{"system": LOINC, "code": "53041-0",
                                 "display": "DNA region of interest panel"}],
                     "text": f"Region studied: {gene}"},
            **ctx.subject_refs(),
            "component": [
                {"code": {"coding": [{"system": LOINC, "code": "48018-6", "display": "Gene studied [ID]"}]},
                 "valueCodeableConcept": gene_map.concept(gene)},
                REFSEQ_COMPONENT,
                {"code": {"coding": [{"system": LOINC, "code": "82121-5", "display": "Allelic read depth"}]},
                 "valueQuantity": {"value": info.get('mean_depth', 0), "unit": "reads per base pair",
                                   "system": "http://unitsofmeasure.org", "code": "[1]"}},
                {"code": tbfhir_code("coverage-breadth", f"Fraction of bases at >= {min_depth}x"),
                 "valueQuantity": {"value": round(float(breadth) * 100, 2), "unit": "%",
                                   "system": "http://unitsofmeasure.org", "code": "%"}},
            ],
        }
        observations.append(with_meta(observation, ctx.profile.meta("region")))
    return observations


def build_lineage_observation(ctx, lineage_info):
    if not lineage_info:
        return None
    lineage = lineage_info.get('lineage', 'unknown')
    if lineage == 'unknown':
        return None
    family = lineage_info.get('family', 'Unknown')
    confidence = lineage_info.get('confidence')
    score = lineage_info.get('score')
    matched = lineage_info.get('matched_snps')
    total = lineage_info.get('total_snps')

    text = f"Mycobacterial Lineage: {lineage}"
    if family != 'Unknown':
        text += f" ({family})"
    if matched is not None and total:
        text += f" - {matched}/{total} barcode SNPs"
    if confidence:
        text += f", confidence {confidence}"

    def attribute(code, display):
        return {"coding": [{"system": CS_LINEAGE_ATTRIBUTE, "code": code, "display": display}],
                "text": display}

    components = []
    if family and family != 'Unknown':
        components.append({"code": attribute("lineage-family", "Lineage family"), "valueString": family})
    if confidence:
        components.append({"code": attribute("lineage-confidence", "Lineage call confidence"),
                           "valueString": str(confidence)})
    for key, value, cast, code, display in (
            ("score", score, float, "lineage-score", "Lineage barcode score"),
            ("matched", matched, int, "barcode-snps-matched", "Barcode SNPs matched"),
            ("total", total, int, "barcode-snps-total", "Barcode SNPs examined")):
        if value is None:
            continue
        try:
            components.append({"code": attribute(code, display),
                               "valueQuantity": {"value": cast(value),
                                                 "system": "http://unitsofmeasure.org", "code": "1"}})
        except (TypeError, ValueError):
            pass

    codings = [{"system": CS_LINEAGE, "code": lineage, "display": f"M. tuberculosis {lineage}"}]
    if ctx.profile.lineage_legacy_coding:
        codings.append({"system": "http://tb-lineage.org", "code": lineage,
                        "display": f"TB Lineage {lineage}"})
    observation = {
        "resourceType": "Observation",
        "id": f"{ctx.fid}-lineage",
        "text": {"status": "generated", "div": div(text)},
        "status": "final",
        "category": [CATEGORY_LAB, CATEGORY_GE],
        "code": {"coding": [{"system": LOINC, "code": "614-8",
                             "display": "Mycobacterial strain [Type] in Isolate by Mycobacterial subtyping"}]},
        "valueCodeableConcept": {"coding": codings, "text": f"Lineage {lineage}"},
        **ctx.subject_refs(),
    }
    if components:
        observation["component"] = components
    return with_meta(observation, ctx.profile.meta("lineage"))

def build_diagnostic_report(ctx, observations):
    (resistance_class, resistance_description, resistant_genes,
     resistant_drugs, indeterminate_drugs) = classify_drug_resistance(observations)
    lineage_info = extract_lineage_info(observations)

    parts = [f"{resistance_class} ({resistance_description})"]
    if resistant_genes:
        parts.append(f"Detected resistance genes: {', '.join(resistant_genes)}")
    if resistant_drugs:
        parts.append(f"Detected drug resistance: {', '.join(resistant_drugs)}")
    if indeterminate_drugs:
        parts.append("NOT ASSESSABLE (coverage of resistance loci not confirmed; "
                     f"susceptibility cannot be inferred): {', '.join(indeterminate_drugs)}")
    if lineage_info:
        parts.append(f"TB {lineage_info} detected")
    parts.append("Reference genome: NC_000962.3")
    conclusion = ". ".join(parts)

    conclusion_codes = []
    if resistance_class:
        entry = {"text": resistance_class}
        coding = ctx.profile.conclusion_coding(resistance_class)
        if coding:
            entry["coding"] = [coding]
        conclusion_codes.append(entry)
    if lineage_info:
        conclusion_codes.append({"text": lineage_info})

    report_html = ("<div xmlns=\"http://www.w3.org/1999/xhtml\">\n"
                   "<h1>TB Genomic Analysis Report</h1>\n"
                   f"<p><strong>Sample ID:</strong> {xhtml_text(ctx.fid)}</p>\n"
                   f"<p><strong>Report Date:</strong> {ctx.timestamp}</p>\n"
                   f"<p><strong>Resistance Classification:</strong> {xhtml_text(resistance_class)}</p>\n"
                   f"<p><strong>Conclusion:</strong> {xhtml_text(conclusion)}</p>\n")
    if observations:
        report_html += "<h2>Detected Mutations</h2><ul>"
        for obs in observations:
            for component in obs.get('component', []):
                text = component.get('valueCodeableConcept', {}).get('text', '')
                if text and text != 'unknown':
                    report_html += f"<li>{xhtml_text(text)}</li>"
        report_html += "</ul>"
    if lineage_info:
        report_html += f"<p><strong>Mycobacterial Lineage:</strong> {xhtml_text(lineage_info)}</p>"
    report_html += "</div>"

    report = {
        "resourceType": "DiagnosticReport",
        "id": f"{ctx.fid}-genomic-report",
        "text": {"status": "generated",
                 "div": div(f"TB genomic analysis report for sample {ctx.fid}. {conclusion}")},
        "status": "final",
        "category": [CATEGORY_GE],
        "code": ctx.profile.report_code(),
        "subject": {"reference": f"Patient/{ctx.fid}-patient"},
        "effectiveDateTime": ctx.timestamp,
        "issued": ctx.timestamp,
        "result": [{"reference": f"Observation/{obs['id']}"} for obs in observations if obs.get('id')],
        "specimen": [{"reference": f"Specimen/{ctx.fid}-specimen"}],
        "conclusion": conclusion,
        "conclusionCode": conclusion_codes,
        "presentedForm": [{"contentType": "text/html", "language": "en-US",
                           "title": "TB Genomic Analysis Report",
                           "data": base64.b64encode(report_html.encode('utf-8')).decode('ascii')}],
    }
    org_id = ctx.performer["reference"].split("/", 1)[1] if ctx.performer else None
    identifier = ctx.profile.identifier("diagnostic-report", f"TB-GEN-{ctx.fid}-{today_compact()}", org_id)
    if identifier:
        report["identifier"] = [identifier]
    if ctx.performer:
        report["performer"] = [dict(ctx.performer)]
    return with_meta(report, ctx.profile.meta("report")), resistance_class


def build_device(profile, pipeline_info):
    version = pipeline_info.get('pipeline_version', 'unknown')
    device_id = f"tb-mutation-analysis-{sanitize_id(version) if version else 'unknown'}"
    properties = []

    def prop_type(code, display):
        return {"coding": [{"system": CS_PIPELINE_PROPERTY, "code": code, "display": display}],
                "text": display}

    if pipeline_info.get('source_repository'):
        properties.append({"type": prop_type('source-repository', 'Source repository'),
                           "valueCode": [{"text": str(pipeline_info['source_repository'])}]})
    for key, code, display in (
            ('filter_min_depth', 'filter-min-depth', 'Minimum read depth filter'),
            ('filter_min_quality', 'filter-min-quality', 'Minimum variant quality filter'),
            ('coverage_min_depth', 'coverage-min-depth', 'Coverage assessability depth threshold'),
            ('coverage_min_breadth', 'coverage-min-breadth', 'Coverage assessability breadth threshold')):
        value = pipeline_info.get(key)
        if value is None:
            continue
        properties.append({"type": prop_type(code, display),
                           "valueQuantity": [{"value": float(value), "system": "http://unitsofmeasure.org",
                                              "code": "1"}]})
    device = {
        "resourceType": "Device",
        "id": device_id,
        "status": "active",
        "deviceName": [{"name": "TBtoFHIR", "type": "manufacturer-name"}],
        "type": {"coding": [{"system": SNOMED, "code": "706689003",
                             "display": "Application program software"}],
                 "text": "Bioinformatics analysis pipeline"},
        "version": [{"type": {"text": "software-version"}, "value": str(version)}],
    }
    if properties:
        device["property"] = properties
    return with_meta(device, profile.meta("device"))


def build_provenance(ctx, targets, device_id, source_description):
    agents = [{
        "type": {"coding": [{"system": "http://terminology.hl7.org/CodeSystem/provenance-participant-type",
                             "code": "assembler", "display": "Assembler"}]},
        "who": {"reference": f"Device/{device_id}", "display": "TBtoFHIR"},
    }]
    if ctx.performer:
        agents.append({
            "type": {"coding": [{"system": "http://terminology.hl7.org/CodeSystem/provenance-participant-type",
                                 "code": "performer", "display": "Performer"}]},
            "who": dict(ctx.performer),
        })
    provenance = {
        "resourceType": "Provenance",
        "id": f"{ctx.fid}-provenance",
        "target": targets,
        "recorded": ctx.timestamp,
        "activity": {"coding": [{"system": "http://terminology.hl7.org/CodeSystem/v3-DataOperation",
                                 "code": "CREATE", "display": "create"}]},
        "agent": agents,
    }
    if source_description:
        provenance["entity"] = [{"role": "derivation", "what": {"display": source_description}}]
    return with_meta(provenance, ctx.profile.meta("provenance"))

def load_sample_json(path=None, directory=None, sample_id=None, suffix=None):
    """One sample's lineage/coverage summary, from a file or (legacy) a directory."""
    if path:
        with open(path, encoding="utf-8") as handle:
            return json.load(handle)
    if directory and os.path.isdir(directory):
        records = {}
        for name in sorted(os.listdir(directory)):
            if name.endswith(suffix):
                with open(os.path.join(directory, name), encoding="utf-8") as handle:
                    data = json.load(handle)
                records[data.get('sample_id', name[: -len(suffix)])] = data
        return lookup_sample(records, sample_id)
    return None


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--input', required=True, help='annotated VCF (.vcf or .vcf.gz)')
    p.add_argument('--sample_id', required=True, help='sample identifier (exact)')
    p.add_argument('--output', required=True, help='output Bundle JSON')
    p.add_argument('--lineage_json', help="this sample's <id>.lineage.json")
    p.add_argument('--coverage_json', help="this sample's <id>.coverage.json")
    p.add_argument('--lineage_dir', help=argparse.SUPPRESS)    # v1 compatibility
    p.add_argument('--coverage_dir', help=argparse.SUPPRESS)   # v1 compatibility
    p.add_argument('--gene_map', help='optional TSV gene<TAB>ncbi_gene_id extending the built-in map')
    p.add_argument('--fhir_profile', default='generic', choices=['generic', 'satusehat'])
    p.add_argument('--identifier_system_base', default=None,
                   help='base URI for business identifiers in the generic profile')
    p.add_argument('--performer_org_id', default=None,
                   help='FHIR id of the performing laboratory Organization (not PHI)')
    p.add_argument('--performer_org_name', default=None)
    p.add_argument('--pipeline_version', default='unknown')
    p.add_argument('--source_repository', default=SOURCE_REPOSITORY)
    p.add_argument('--filter_min_depth', type=float, default=None)
    p.add_argument('--filter_min_quality', type=float, default=None)
    p.add_argument('--coverage_min_depth', type=float, default=None)
    p.add_argument('--coverage_min_breadth', type=float, default=None)
    p.add_argument('--legacy-annotation-table', dest='legacy_annotation_table', action='store_true',
                   help='Tolerate a one-row-per-pair annotation table by applying the most severe '
                        'grade to every listed drug. Over-calls resistance.')
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    if not os.path.exists(args.input) or os.path.getsize(args.input) == 0:
        sys.exit(f"ERROR: input VCF missing or empty: {args.input}")

    profile = get_profile(args.fhir_profile, args.identifier_system_base)
    fid = fhir_id(args.sample_id)
    performer = None
    if args.performer_org_id:
        performer = {"reference": f"Organization/{args.performer_org_id}"}
        if args.performer_org_name:
            performer["display"] = args.performer_org_name
    ctx = SampleContext(args.sample_id, fid, profile, performer, now_iso())
    gene_map = GeneMap(args.gene_map)

    lineage_info = load_sample_json(args.lineage_json, args.lineage_dir, args.sample_id, '.lineage.json')
    coverage_info = load_sample_json(args.coverage_json, args.coverage_dir, args.sample_id, '.coverage.json')
    for label, info in (("lineage", lineage_info), ("coverage", coverage_info)):
        if info and info.get('sample_id') and fhir_id(info['sample_id']) != fid:
            sys.exit(f"ERROR: {label} file is for sample {info['sample_id']!r}, not {args.sample_id!r}")

    print(f"Lineage: {'found' if lineage_info else 'none for this sample'}")
    if coverage_info and coverage_info.get('coverage_assessed'):
        assessable = sum(1 for d in coverage_info.get('drugs', {}).values() if d.get('assessable'))
        print(f"Coverage data found: {assessable} drugs assessable")
    else:
        reason = (coverage_info or {}).get('reason', 'no coverage file for this sample')
        print(f"No usable coverage data ({reason}); non-resistant drugs will be reported as "
              f"INDETERMINATE, not susceptible")

    variant_observations = []
    resistant_drugs = set()
    skipped = []
    legacy_collapsed = 0
    for idx, record in enumerate(read_vcf(args.input)):
        try:
            drug_pairs = resolve_variant_drug_pairs(record, allow_legacy=args.legacy_annotation_table)
        except AnnotationError as exc:
            print(f"\nFATAL: {exc}\n")
            sys.exit(2)
        if args.legacy_annotation_table and len(drug_pairs) > 1 and \
                len(get_variant_list(record, 'WHO_CLASSIFICATION')) == 1:
            legacy_collapsed += 1
        for drug, grade, _v in drug_pairs:
            if is_resistance_grade(grade):
                name = drug.lower().strip()
                resistant_drugs.add(name)
                if 'fluoroquinolone' in name:
                    resistant_drugs.update({'moxifloxacin', 'levofloxacin', 'ciprofloxacin', 'ofloxacin'})
        try:
            variant_observations.append(build_variant_observation(ctx, idx, record, drug_pairs, gene_map))
        except (ValueError, KeyError, TypeError) as exc:
            skipped.append((idx, f"{record.CHROM}:{record.POS}", str(exc)))

    panel, panel_summary = build_susceptibility_panel(ctx, resistant_drugs, coverage_info)
    observations = variant_observations + [panel] + build_region_observations(ctx, coverage_info, gene_map)
    lineage_obs = build_lineage_observation(ctx, lineage_info)
    if lineage_obs:
        observations.append(lineage_obs)

    report, resistance_class = build_diagnostic_report(ctx, observations)
    device = build_device(profile, {
        'pipeline_version': args.pipeline_version, 'source_repository': args.source_repository,
        'filter_min_depth': args.filter_min_depth, 'filter_min_quality': args.filter_min_quality,
        'coverage_min_depth': args.coverage_min_depth, 'coverage_min_breadth': args.coverage_min_breadth})
    targets = [{"reference": f"DiagnosticReport/{report['id']}"}] + \
              [{"reference": f"Observation/{obs['id']}"} for obs in observations]
    provenance = build_provenance(ctx, targets, device['id'],
                                  f"Annotated variant calls for {fid} ({os.path.basename(args.input)})")

    resources = point_references_at_entries([report] + observations + [device, provenance])
    bundle = {
        "resourceType": "Bundle",
        "id": bundle_id("genomic", fid),
        "type": "collection",
        "timestamp": ctx.timestamp,
        "entry": [{"fullUrl": full_url(r["resourceType"], r["id"]), "resource": r} for r in resources],
    }
    with open(args.output, 'w', encoding='utf-8') as out:
        json.dump(bundle, out, indent=2, ensure_ascii=False)

    print(f"\nSummary for {args.sample_id} (FHIR id {fid}, profile {profile.name}):")
    print(f"  variants written        : {len(variant_observations)}")
    print(f"  records skipped (error) : {len(skipped)}")
    for idx, location, message in skipped[:10]:
        print(f"      - record {idx} at {location}: {message}")
    print(f"  resistant drugs detected: {', '.join(sorted(resistant_drugs)) or 'none'}")
    print(f"  susceptibility panel    : {panel_summary}")
    print(f"  classification          : {resistance_class}")
    if legacy_collapsed:
        print(f"  WARNING: {legacy_collapsed} multi-drug alleles had a single collapsed grade "
              f"applied to every listed drug (--legacy-annotation-table).")
    if gene_map.unmapped:
        print(f"  genes without an NCBI GeneID: {len(gene_map.unmapped)} (text only): "
              f"{', '.join(sorted(gene_map.unmapped))}")
    if skipped:
        sys.exit(f"ERROR: {len(skipped)} variant record(s) could not be converted")


if __name__ == "__main__":
    main()
