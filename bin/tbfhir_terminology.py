
LOINC = "http://loinc.org"
SNOMED = "http://snomed.info/sct"
SEQUENCE_ONTOLOGY = "http://www.sequenceontology.org"
NCBI_GENE = "https://www.ncbi.nlm.nih.gov/gene"
REFSEQ = "http://www.ncbi.nlm.nih.gov/refseq"
HGVS = "http://varnomen.hgvs.org"
REFERENCE_SEQUENCE = "NC_000962.3"

_DRUGS = {
    'rifampicin': ('29175007', 'Product containing rifampicin'),
    'isoniazid': ('81335000', 'Product containing isoniazid'),
    'pyrazinamide': ('13592004', 'Product containing pyrazinamide'),
    'ethambutol': ('24450004', 'Product containing ethambutol'),
    'streptomycin': ('40877002', 'Product containing streptomycin'),
    'fluoroquinolone': ('1010205001', 'Medicinal product containing fluoroquinolone and acting as antibacterial agent'),
    'levofloxacin': ('96087006', 'Product containing levofloxacin'),
    'moxifloxacin': ('371296007', 'Product containing moxifloxacin'),
    'ofloxacin': ('96086002', 'Product containing ofloxacin'),
    'ciprofloxacin': ('7577004', 'Product containing ciprofloxacin'),
    'gatifloxacin': ('371238005', 'Product containing gatifloxacin'),
    'amikacin': ('48836000', 'Product containing amikacin'),
    'kanamycin': ('71451001', 'Product containing kanamycin'),
    'capreomycin': ('14170004', 'Product containing capreomycin'),
    'ethionamide': ('414148003', 'Product containing ethionamide'),
    'linezolid': ('125695009', 'Product containing linezolid'),
    'bedaquiline': ('714087008', 'Product containing bedaquiline'),
    'clofazimine': ('72924009', 'Product containing clofazimine'),
    'delamanid': ('714098009', 'Product containing delamanid'),
    'pretomanid': ('789321008', 'Product containing pretomanid'),
    'cycloserine': ('51334008', 'Product containing cycloserine'),
    'terizidone': ('1395908004', 'Product containing terizidone'),
    'para-aminosalicylic_acid': ('417238004', 'Product containing aminosalicylic acid'),
}


def get_drug_snomed_mapping(drug_name):
    entry = _DRUGS.get(drug_name.lower().strip().replace(' ', '_'))
    return {'code': entry[0], 'display': entry[1]} if entry else None


_EFFECTS = {
    'missense_variant': ('SO:0001583', 'missense_variant'),
    'synonymous_variant': ('SO:0001819', 'synonymous_variant'),
    'stop_gained': ('SO:0001587', 'stop_gained'),
    'stop_lost': ('SO:0001578', 'stop_lost'),
    'frameshift_variant': ('SO:0001589', 'frameshift_variant'),
    'inframe_insertion': ('SO:0001821', 'inframe_insertion'),
    'inframe_deletion': ('SO:0001822', 'inframe_deletion'),
    'splice_site_variant': ('SO:0001629', 'splice_site_variant'),
    'upstream_gene_variant': ('SO:0001631', 'upstream_gene_variant'),
    'downstream_gene_variant': ('SO:0001632', 'downstream_gene_variant'),
    'intergenic_variant': ('SO:0001628', 'intergenic_variant'),
    'intron_variant': ('SO:0001627', 'intron_variant'),
    '5_prime_utr_variant': ('SO:0001623', '5_prime_UTR_variant'),
    '3_prime_utr_variant': ('SO:0001624', '3_prime_UTR_variant'),
    'start_lost': ('SO:0002012', 'start_lost'),
    'stop_retained_variant': ('SO:0001567', 'stop_retained_variant'),
    'protein_altering_variant': ('SO:0001818', 'protein_altering_variant'),
    'coding_sequence_variant': ('SO:0001580', 'coding_sequence_variant'),
    'non_coding_transcript_variant': ('SO:0001619', 'non_coding_transcript_variant'),
    'regulatory_region_variant': ('SO:0001566', 'regulatory_region_variant'),
    'loss_of_function_variant': ('SO:0002054', 'loss_of_function_variant'),
    'non_coding_transcript_exon_variant': ('SO:0001792', 'non_coding_transcript_exon_variant'),
    'initiator_codon_variant': ('SO:0001582', 'initiator_codon_variant'),
    'feature_ablation': ('SO:0001879', 'feature_ablation'),
    'feature_truncation': ('SO:0001906', 'feature_truncation'),
    'disruptive_inframe_deletion': ('SO:0001826', 'disruptive_inframe_deletion'),
    'disruptive_inframe_insertion': ('SO:0001824', 'disruptive_inframe_insertion'),
    'splice_acceptor_variant': ('SO:0001574', 'splice_acceptor_variant'),
    'splice_donor_variant': ('SO:0001575', 'splice_donor_variant'),
    'splice_region_variant': ('SO:0001630', 'splice_region_variant'),
    'transcript_ablation': ('SO:0001893', 'transcript_ablation'),
    'rrna_variant': ('SO:0001637', 'rRNA_gene_variant'),
}


def get_effect_so_mapping(effect):
    entry = _EFFECTS.get(effect.lower().strip())
    return {'code': entry[0], 'display': entry[1]} if entry else None


_BUILTIN_NCBI_GENE_IDS = {
    'katG': '885638', 'rpoB': '888164', 'inhA': '886523', 'gyrA': '887105',
    'gyrB': '887081', 'embB': '886126', 'pncA': '888260', 'rpsL': '888259',
    'rrs': '2700429', 'ethA': '886175', 'ahpC': '885717', 'fabG1': '886551',
    'eis': '885903', 'tlyA': '885396', 'alr': '887634', 'ddlA': '888415',
    'rpoC': '888177', 'embA': '886123', 'embC': '886112', 'ubiA': '886129',
    'ndh': '885746', 'ndhA': '886430', 'Rv0678': '888235', 'gid': '886243',
}


class GeneMap:

    def __init__(self, extra_path=None):
        self.ids = dict(_BUILTIN_NCBI_GENE_IDS)
        self.unmapped = set()
        self.loaded_external = False
        if extra_path:
            self.load(extra_path)

    def load(self, path):
        with open(path, encoding="utf-8") as handle:
            for line in handle:
                if line.startswith('#') or not line.strip():
                    continue
                parts = line.rstrip('\n').split('\t')
                if len(parts) >= 2 and parts[0].strip() and parts[1].strip().isdigit():
                    self.ids[parts[0].strip()] = parts[1].strip()
        self.loaded_external = True

    def concept(self, gene):
        gene_id = self.ids.get(gene)
        if gene_id:
            return {"coding": [{"system": NCBI_GENE, "code": gene_id, "display": gene}],
                    "text": gene}
        self.unmapped.add(gene)
        return {"text": gene}


def get_drug_panel_config():
    return [
        ('rifampicin', '89489-9', 'rifAMPin [Susceptibility] by Genotype method'),
        ('isoniazid', '89488-1', 'Isoniazid [Susceptibility] by Genotype method'),
        ('ethambutol', '89491-5', 'Ethambutol [Susceptibility] by Genotype method'),
        ('pyrazinamide', '92242-7', 'Pyrazinamide [Susceptibility] by Genotype method'),
        ('moxifloxacin', '96112-8', 'Moxifloxacin [Susceptibility] by Genotype method'),
        ('levofloxacin', '20629-2', 'levoFLOXacin [Susceptibility]'),
        ('bedaquiline', '96107-8', 'Bedaquiline [Susceptibility] by Genotype method'),
        ('delamanid', '96109-4', 'Delamanid [Susceptibility] by Genotype method'),
        ('pretomanid', '93850-6', 'Pretomanid [Susceptibility]'),
        ('streptomycin', '96114-4', 'Streptomycin [Susceptibility] by Genotype method'),
        ('amikacin', '89484-0', 'Amikacin [Susceptibility] by Genotype method'),
        ('kanamycin', '89482-4', 'Kanamycin [Susceptibility] by Genotype method'),
        ('capreomycin', '89483-2', 'Capreomycin [Susceptibility] by Genotype method'),
        ('clofazimine', '96108-6', 'Clofazimine [Susceptibility] by Genotype method'),
        ('ethionamide', '96110-2', 'Ethionamide [Susceptibility] by Genotype method'),
        ('linezolid', '96111-0', 'Linezolid [Susceptibility] by Genotype method'),
        ('cycloserine', '103959-3', 'cycloSERINE [Susceptibility] by Genotype method'),
    ]
