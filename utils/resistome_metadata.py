from collections import defaultdict
import pandas as pd

VCF_ANN_DIR = '/storage/NFS/ANALYSIS/DNAseq/VCF_ann'
VCF_DIR = '/storage/NFS/ANALYSIS/DNAseq/VCF_ann'
BAM_DIR = '/storage/NFS/ANALYSIS/DNAseq/VCF_ann'

CURRENT_METADATA_FILEPATH = '/storage/NFS/ROTATION_PROJECT/daisy/Winzeler_databases/Sequenced_Clone_Metadata_12-15-2024.csv'

BAD_SAMPLE_PREFIXES = ['Christin-179', 'BMGF-Winzeler-080', 'Wirth-Thies', '2-','Fidock-Andrew', 'GHDDI-ActinomycinA-Dd2',
                       'Wirth-Piktransfect-2-1', 'fidock', 'clador', 'Wirth-Halofuginone-eIK1-I2P-revertant-clone-F7',
                       'Fidock-KAE609-P19.04', 'Winzeler-HT2-Toxin-2-3D7-T1', 'Stokes-V1S-C580Y-5X-WLW-R1', 'Wirth-DSM265-3D7-1']

def parse_metadata_table(metadata_filepath):

    df = pd.read_csv(metadata_filepath, converters={i: str for i in range(100)})

    vcf_clones_dict = defaultdict(list); clone_vcf_dict = {}
    clone_parent_dict = {}; clone_compound_dict = {}; clone_strain_dict = {}

    for row in df.iterrows():
        
        status, raw_dir, analyzed_dir, bam, lab, strain, compound, aliases, final_sample, cloned, clone, parent, \
            target, driver_known, driver_gene_str, driver_mutation_str, pmid, \
            seq_platform, _, parent_ec50, clone_ec50, ec50_foldshift, \
            vcf, vcf_date, vcf_vers, parent_vcf, parent_in_same_vcf = row[1].__array__()[:27]
        
        # Only include resistant samples with known parents, and their parents
        if parent == '' or parent == '?' or (compound in ['NON_SELECTION', 'FIELD_SAMPLE', '']):
            continue       
        
        # Further restrict to samples with grouped VCFs
        if vcf == '' or parent_in_same_vcf == 'FALSE':
            continue
        
        # Exclude "bad" samples (TODO)
        bad_sample = False
        for prefix in BAD_SAMPLE_PREFIXES:
            if clone.startswith(prefix):
                bad_sample = True
        
        if bad_sample is True:
            continue
        
        clone_parent_dict[clone] = parent
        clone_compound_dict[clone] = compound
        clone_strain_dict[clone] = strain 
        
        vcf_clones_dict[vcf].append(clone)
        clone_vcf_dict[clone] = vcf
    
    return clone_parent_dict, clone_compound_dict, clone_strain_dict, vcf_clones_dict, clone_vcf_dict

# =============================
# Metadata objects to be used
# =============================

SAMPLE_PARENT_DICT, SAMPLE_COMPOUND_DICT, SAMPLE_STRAIN_DICT, VCF_SAMPLES_DICT, SAMPLE_VCF_DICT = parse_metadata_table(CURRENT_METADATA_FILEPATH)

def get_sample_VCF_path(sample, annotated=True):
    vcf = SAMPLE_VCF_DICT[sample]
    if annotated:
        return f'{VCF_ANN_DIR}/{vcf}.ann.txt'
    else:
        return f'{VCF_DIR}/{vcf}'

def 