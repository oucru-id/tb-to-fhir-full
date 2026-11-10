
import re


def _drug_group_for(drug_display):
    if 'rifampicin' in drug_display or 'rifampin' in drug_display:
        return 'rifampicin'
    if 'isoniazid' in drug_display:
        return 'isoniazid'
    if 'ethambutol' in drug_display:
        return 'ethambutol'
    if 'pyrazinamide' in drug_display:
        return 'pyrazinamide'
    if 'streptomycin' in drug_display:
        return 'streptomycin'
    if 'ethionamide' in drug_display:
        return 'ethionamide'
    if any(fq in drug_display for fq in ['levofloxacin', 'moxifloxacin', 'ofloxacin', 'ciprofloxacin']):
        return 'fluoroquinolone'
    if any(sli in drug_display for sli in ['amikacin', 'kanamycin', 'capreomycin']):
        return 'second_line_injectable'
    if any(ga in drug_display for ga in ['bedaquiline', 'linezolid']):
        return 'group_a'
    return drug_display


def classify_drug_resistance(observations):

    resistant_drugs_groups = set()
    indeterminate_groups = set()
    susceptible_groups = set()
    detected_drugs = set()
    indeterminate_drugs = set()
    detected_genes = set()

    for obs in observations:
        codes = obs.get('code', {}).get('coding', [])
        is_panel = any(c.get('code') == '89486-5' for c in codes)

        if is_panel:
            components = obs.get('component', [])
            for comp in components:
                value_coding = comp.get('valueCodeableConcept', {}).get('coding', [])
                value_codes = {vc.get('code') for vc in value_coding}

                drug_display = comp.get('code', {}).get('coding', [{}])[0].get('display', '').lower()
                drug_name = re.sub(r'\s*\[susceptibility\]\s*', ' ', drug_display)
                drug_name = re.sub(r'\s*by genotype method\s*', ' ', drug_name).strip()
                group = _drug_group_for(drug_display)

                if 'LA6676-6' in value_codes:
                    detected_drugs.add(drug_name)
                    resistant_drugs_groups.add(group)
                elif 'LA24225-7' in value_codes:
                    susceptible_groups.add(group)
                else:
                    indeterminate_drugs.add(drug_name)
                    indeterminate_groups.add(group)

    for obs in observations:
        codes = obs.get('code', {}).get('coding', [])
        is_variant = any(c.get('code') == '69548-6' for c in codes)
        
        if is_variant:
            components = obs.get('component', [])

            current_gene = ""
            has_resistance_grade = False

            for component in components:
                code_display = component.get('code', {}).get('coding', [{}])[0].get('display', '').lower()

                if 'gene studied' in code_display:
                    current_gene = component.get('valueCodeableConcept', {}).get('text', '')

                if 'genetic variation clinical significance' in code_display:
                    significance = component.get('valueCodeableConcept', {}).get('text', '')
                    if significance and "Assoc w R" in significance:
                        has_resistance_grade = True

            if current_gene and has_resistance_grade:
                detected_genes.add(current_gene)
    
    has_rif = 'rifampicin' in resistant_drugs_groups
    has_inh = 'isoniazid' in resistant_drugs_groups
    has_fq = 'fluoroquinolone' in resistant_drugs_groups
    has_group_a = 'group_a' in resistant_drugs_groups

    rif_indeterminate = 'rifampicin' in indeterminate_groups
    inh_indeterminate = 'isoniazid' in indeterminate_groups

    is_mdr = has_rif and has_inh
    is_rr = has_rif

    if not resistant_drugs_groups:
        if not susceptible_groups:
            classification = "Indeterminate"
            description = ("No resistance-associated mutation detected, but no drug could be "
                           "assessed: coverage of the resistance loci was not confirmed. This "
                           "is NOT a susceptible result")
        elif indeterminate_groups:
            classification = "No resistance detected - partial"
            description = (f"No resistance-associated mutation detected in "
                           f"{len(susceptible_groups)} assessable drug group(s); "
                           f"{len(indeterminate_groups)} could not be assessed "
                           f"({', '.join(sorted(indeterminate_drugs))})")
        else:
            classification = "No resistance detected"
            description = ("No resistance-associated mutation detected across all assessed "
                           "drugs, with adequate coverage of the resistance loci")

        return (classification, description,
                sorted(list(detected_genes)), sorted(list(detected_drugs)),
                sorted(list(indeterminate_drugs)))

    if rif_indeterminate and not has_rif:
        classification = "Drug-resistant - rifampicin not assessable"
        description = (f"Resistance detected to: {', '.join(sorted(resistant_drugs_groups))}. "
                       f"Rifampicin could not be assessed, so RR/MDR/pre-XDR/XDR "
                       f"classification cannot be determined")
        return (classification, description,
                sorted(list(detected_genes)), sorted(list(detected_drugs)),
                sorted(list(indeterminate_drugs)))

    classification = "Drug-resistant"
    description = f"Resistance to: {', '.join(sorted(resistant_drugs_groups))}"

    if (is_mdr or is_rr) and has_fq and has_group_a:
        classification = "XDR-TB"
        description = "Extensively drug-resistant tuberculosis (MDR/RR + FQ + Group A)"
    elif (is_mdr or is_rr) and has_fq:
        classification = "Pre-XDR-TB"
        description = "Pre-extensively drug-resistant tuberculosis (MDR/RR + FQ)"
    elif has_rif and has_inh:
        classification = "MDR-TB"
        description = "Multidrug-resistant tuberculosis"
    elif has_rif and not has_inh:
        classification = "RR-TB"
        if inh_indeterminate:
            description = ("Rifampicin-resistant tuberculosis. Isoniazid could not be "
                           "assessed, so MDR-TB cannot be excluded")
        else:
            description = "Rifampicin-resistant tuberculosis"
    elif has_inh and not has_rif:
        classification = "HR-TB"
        description = "Isoniazid-resistant, rifampicin-susceptible tuberculosis"
    elif len(resistant_drugs_groups) == 1:
        if 'streptomycin' in resistant_drugs_groups:
            classification = "Streptomycin-resistant TB"
            description = "Streptomycin mono-resistant tuberculosis"
        elif 'ethionamide' in resistant_drugs_groups:
            classification = "Ethionamide-resistant TB"
            description = "Ethionamide mono-resistant tuberculosis"
        elif 'pyrazinamide' in resistant_drugs_groups:
            classification = "Pyrazinamide-resistant TB"
            description = "Pyrazinamide mono-resistant tuberculosis"
        elif 'ethambutol' in resistant_drugs_groups:
            classification = "Ethambutol-resistant TB"
            description = "Ethambutol mono-resistant tuberculosis"
        elif 'fluoroquinolone' in resistant_drugs_groups:
            has_cipro = any('ciprofloxacin' in d for d in detected_drugs)
            has_other_fq = any(fq in d for d in detected_drugs for fq in ['levofloxacin', 'moxifloxacin', 'ofloxacin'])
            
            if has_cipro and not has_other_fq:
                classification = "Ciprofloxacin-resistant TB"
                description = "Ciprofloxacin mono-resistant tuberculosis"
            else:
                classification = "Drug-resistant"
                description = f"Resistance to: {', '.join(sorted(resistant_drugs_groups))}"
        else:
            classification = "Drug-resistant"
            description = f"Resistance to: {', '.join(sorted(resistant_drugs_groups))}"
    else:
        classification = "Drug-resistant"
        description = f"Resistance to: {', '.join(sorted(resistant_drugs_groups))}"

    return (classification, description,
            sorted(list(detected_genes)), sorted(list(detected_drugs)),
            sorted(list(indeterminate_drugs)))

def extract_lineage_info(observations):
    for obs in observations:
        codes = obs.get('code', {}).get('coding', [])
        for code in codes:
            if code.get('code') == '614-8':
                lineage_text = obs.get('valueCodeableConcept', {}).get('text', '')
                if lineage_text and lineage_text != 'unknown':
                    return lineage_text

        components = obs.get('component', [])
        for component in components:
            code_display = component.get('code', {}).get('coding', [{}])[0].get('display', '').lower()
            if 'mycobacterial strain' in code_display or 'lineage' in code_display:
                lineage_text = component.get('valueCodeableConcept', {}).get('text', '')
                if lineage_text and lineage_text != 'unknown':
                    return lineage_text
    return None

RESISTANCE_CONCLUSION_SNOMED = {
    "Drug-resistant - rifampicin not assessable": ("413556004", "Antibiotic resistant tuberculosis"),
    "RR-TB": ("415345001", "Rifampicin resistant tuberculosis"),
    "HR-TB": ("414546009", "Isoniazid resistant tuberculosis"),
    "MDR-TB": ("423092005", "Multidrug resistant tuberculosis"),
    "XDR-TB": ("710106005", "Extensively drug resistant tuberculosis"),
    "Streptomycin-resistant TB": ("415622003", "Streptomycin resistant tuberculosis"),
    "Ethionamide-resistant TB": ("414149006", "Ethionamide resistant tuberculosis"),
    "Pyrazinamide-resistant TB": ("415222009", "Pyrazinamide resistant tuberculosis"),
    "Ciprofloxacin-resistant TB": ("413852006", "Ciprofloxacin resistant tuberculosis"),
    "Ethambutol-resistant TB": ("414146004", "Ethambutol resistant tuberculosis"),
    "Drug-resistant": ("413556004", "Antibiotic resistant tuberculosis"),
}


def snomed_conclusion_coding(resistance_class):
    entry = RESISTANCE_CONCLUSION_SNOMED.get(resistance_class)
    if not entry:
        return None
    return {"system": "http://snomed.info/sct", "code": entry[0], "display": entry[1]}
