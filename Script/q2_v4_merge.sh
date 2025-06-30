#!/usr/bin/env bash

main(){
    TIME_START=`date +%s`
    echo -e "Starting ...\n  `date -d@${TIME_START}`\n"

    set_variables
    create_folders

    create_manifest_file
    create_metadata_file

#    prepare_taxonomic_classifier

    merge_denoised_datasets
    generate_phylogenetic_tree
    taxonomic_analysis
    export_qiime_artifacts

    filter_datasets
    export_qiime_artifacts_filt
    decontaminate_datasets
    export_qiime_artifacts_decont

    define_sample_groups

    ancombc2_differential_abundance_analysis_seed
    ancombc2_differential_abundance_analysis_gnotobiotic
    ancombc2_differential_abundance_analysis_field

    TIME_END=`date +%s`
    echo -e "\nJob successfully processed.\n  `date -d@${TIME_END}`"
    echo "  Execution time was `expr ${TIME_END} - ${TIME_START}` s."
}

set_variables(){
    echo '[ ] Setting variables'

    THREADS=40
    THREADS_CLASSIFIER=20

    BASE_FOLDER='PRJNA1282304'

    RUNDATE_R4='2022-05-09'
    RUNDATE_R5='2023-06-06'
    RUNDATE_R6='2025-04-07'

    RAW_READS_FOLDER_R4="${BASE_FOLDER}/run_220509"
    RAW_READS_FOLDER_R5="${BASE_FOLDER}/run_230606"
    RAW_READS_FOLDER_R6="${BASE_FOLDER}/run_250407"

    CONF_FOLDER="${BASE_FOLDER}/etc_v4_merg"
    MANIFEST_FILE='sample-manifest.csv'
    METADATA_FILE='sample-metadata.tsv'

    REGION='V4'
    PRIMER_NAME='U515F-Bakt_805R' # Ref.: ProbeBase (https://probebase.net)
    FORWARD_PRIMER='GTGCCAGCMGCCGCGGTAA' # Full name (Alm et al.): S-*-Univ-0515-a-S-19
    REVERSE_PRIMER='GACTACHVGGGTATCTAATCC' # Full name (Alm et al.): S-D-Bact-0785-a-A-21
    FORWARD_PRIMER_LEN=`echo ${FORWARD_PRIMER} | tr -d '\n' | wc -c`
    REVERSE_PRIMER_LEN=`echo ${REVERSE_PRIMER} | tr -d '\n' | wc -c`

    DATA_FOLDER_R4='data_v4_run4'
    DATA_FOLDER_R5='data_v4_run5'
    DATA_FOLDER_R6='data_v4_run6'

    DATA_FOLDER="${BASE_FOLDER}/data_v4_merg"
    BIOM_FOLDER="${DATA_FOLDER}/biom"

    DATABASES_FOLDER="${BASE_FOLDER}/databases"
    REF_DB_VERSION='silva138.2'
    REF_DB_FOLDER="${DATABASES_FOLDER}/${REF_DB_VERSION}"

    CLASSIFIER_TYPE='Default_weighted'   # Choices('Default', 'Default_weighted', 'Custom_16S', 'Custom_V4')
    # For custom classifiers:
    DEREP_MODE='uniq'   # Choices('uniq', 'lca', 'majority', 'super')

    echo -e '[X] Setting variables\n'
}

create_folders(){
    echo '[ ] Creating folders'

    mkdir -p "${CONF_FOLDER}"

    mkdir -p "${DATA_FOLDER}"
    mkdir -p "${BIOM_FOLDER}"

    mkdir -p "${DATABASES_FOLDER}"
    mkdir -p "${REF_DB_FOLDER}"

    echo -e '[X] Creating folders\n'
}

create_manifest_file(){
    echo '[ ] Creating sample-manifest.csv file'

    cat > "${CONF_FOLDER}/${MANIFEST_FILE}" <<EOF
sample-id,absolute-filepath,direction
SPB,${RAW_READS_FOLDER}/10_S10_L001_R1_001.fastq.gz,forward
SPB,${RAW_READS_FOLDER}/10_S10_L001_R2_001.fastq.gz,reverse
SPT,${RAW_READS_FOLDER}/12_S12_L001_R1_001.fastq.gz,forward
SPT,${RAW_READS_FOLDER}/12_S12_L001_R2_001.fastq.gz,reverse
SE1620,${RAW_READS_FOLDER}/15_S15_L001_R1_001.fastq.gz,forward
SE1620,${RAW_READS_FOLDER}/15_S15_L001_R2_001.fastq.gz,reverse
SPM20Q,${RAW_READS_FOLDER}/23_S23_L001_R1_001.fastq.gz,forward
SPM20Q,${RAW_READS_FOLDER}/23_S23_L001_R2_001.fastq.gz,reverse
TBR,${RAW_READS_FOLDER}/24_S24_L001_R1_001.fastq.gz,forward
TBR,${RAW_READS_FOLDER}/24_S24_L001_R2_001.fastq.gz,reverse
TDR,${RAW_READS_FOLDER}/25_S25_L001_R1_001.fastq.gz,forward
TDR,${RAW_READS_FOLDER}/25_S25_L001_R2_001.fastq.gz,reverse
TBS,${RAW_READS_FOLDER}/26_S26_L001_R1_001.fastq.gz,forward
TBS,${RAW_READS_FOLDER}/26_S26_L001_R2_001.fastq.gz,reverse
BSC2021_27,${RAW_READS_FOLDER}/27_S27_L001_R1_001.fastq.gz,forward
BSC2021_27,${RAW_READS_FOLDER}/27_S27_L001_R2_001.fastq.gz,reverse
C_ex_30,${RAW_READS_FOLDER}/30_S30_L001_R1_001.fastq.gz,forward
C_ex_30,${RAW_READS_FOLDER}/30_S30_L001_R2_001.fastq.gz,reverse
TKS,${RAW_READS_FOLDER}/33_S33_L001_R1_001.fastq.gz,forward
TKS,${RAW_READS_FOLDER}/33_S33_L001_R2_001.fastq.gz,reverse
AE_498,${RAW_READS_FOLDER}/34_S34_L001_R1_001.fastq.gz,forward
AE_498,${RAW_READS_FOLDER}/34_S34_L001_R2_001.fastq.gz,reverse
TRI_7121,${RAW_READS_FOLDER}/35_S35_L001_R1_001.fastq.gz,forward
TRI_7121,${RAW_READS_FOLDER}/35_S35_L001_R2_001.fastq.gz,reverse
TRI_445,${RAW_READS_FOLDER}/36_S36_L001_R1_001.fastq.gz,forward
TRI_445,${RAW_READS_FOLDER}/36_S36_L001_R2_001.fastq.gz,reverse
TRI_448,${RAW_READS_FOLDER}/37_S37_L001_R1_001.fastq.gz,forward
TRI_448,${RAW_READS_FOLDER}/37_S37_L001_R2_001.fastq.gz,reverse
TRI_6263,${RAW_READS_FOLDER}/38_S38_L001_R1_001.fastq.gz,forward
TRI_6263,${RAW_READS_FOLDER}/38_S38_L001_R2_001.fastq.gz,reverse
TRI_1998,${RAW_READS_FOLDER}/39_S39_L001_R1_001.fastq.gz,forward
TRI_1998,${RAW_READS_FOLDER}/39_S39_L001_R2_001.fastq.gz,reverse
TRI_4168,${RAW_READS_FOLDER}/47_S47_L001_R1_001.fastq.gz,forward
TRI_4168,${RAW_READS_FOLDER}/47_S47_L001_R2_001.fastq.gz,reverse
TRI_4323,${RAW_READS_FOLDER}/48_S48_L001_R1_001.fastq.gz,forward
TRI_4323,${RAW_READS_FOLDER}/48_S48_L001_R2_001.fastq.gz,reverse
TRI_10765,${RAW_READS_FOLDER}/49_S49_L001_R1_001.fastq.gz,forward
TRI_10765,${RAW_READS_FOLDER}/49_S49_L001_R2_001.fastq.gz,reverse
AE_246,${RAW_READS_FOLDER}/50_S50_L001_R1_001.fastq.gz,forward
AE_246,${RAW_READS_FOLDER}/50_S50_L001_R2_001.fastq.gz,reverse
TRI_4101,${RAW_READS_FOLDER}/51_S51_L001_R1_001.fastq.gz,forward
TRI_4101,${RAW_READS_FOLDER}/51_S51_L001_R2_001.fastq.gz,reverse
CH_SO,${RAW_READS_FOLDER}/52_S52_L001_R1_001.fastq.gz,forward
CH_SO,${RAW_READS_FOLDER}/52_S52_L001_R2_001.fastq.gz,reverse
CH_RO,${RAW_READS_FOLDER}/53_S53_L001_R1_001.fastq.gz,forward
CH_RO,${RAW_READS_FOLDER}/53_S53_L001_R2_001.fastq.gz,reverse
SE16_SO,${RAW_READS_FOLDER}/56_S56_L001_R1_001.fastq.gz,forward
SE16_SO,${RAW_READS_FOLDER}/56_S56_L001_R2_001.fastq.gz,reverse
SE16_RO,${RAW_READS_FOLDER}/57_S57_L001_R1_001.fastq.gz,forward
SE16_RO,${RAW_READS_FOLDER}/57_S57_L001_R2_001.fastq.gz,reverse
SE16_SH,${RAW_READS_FOLDER}/58_S58_L001_R1_001.fastq.gz,forward
SE16_SH,${RAW_READS_FOLDER}/58_S58_L001_R2_001.fastq.gz,reverse
SE16_SP,${RAW_READS_FOLDER}/59_S59_L001_R1_001.fastq.gz,forward
SE16_SP,${RAW_READS_FOLDER}/59_S59_L001_R2_001.fastq.gz,reverse
BRAJua,${RAW_READS_FOLDER}/5_S5_L001_R1_001.fastq.gz,forward
BRAJua,${RAW_READS_FOLDER}/5_S5_L001_R2_001.fastq.gz,reverse
TDR_RO,${RAW_READS_FOLDER}/61_S61_L001_R1_001.fastq.gz,forward
TDR_RO,${RAW_READS_FOLDER}/61_S61_L001_R2_001.fastq.gz,reverse
TDR_SH,${RAW_READS_FOLDER}/62_S62_L001_R1_001.fastq.gz,forward
TDR_SH,${RAW_READS_FOLDER}/62_S62_L001_R2_001.fastq.gz,reverse
TDR_SP,${RAW_READS_FOLDER}/63_S63_L001_R1_001.fastq.gz,forward
TDR_SP,${RAW_READS_FOLDER}/63_S63_L001_R2_001.fastq.gz,reverse
SPT_SO,${RAW_READS_FOLDER}/64_S64_L001_R1_001.fastq.gz,forward
SPT_SO,${RAW_READS_FOLDER}/64_S64_L001_R2_001.fastq.gz,reverse
SPT_RO,${RAW_READS_FOLDER}/65_S65_L001_R1_001.fastq.gz,forward
SPT_RO,${RAW_READS_FOLDER}/65_S65_L001_R2_001.fastq.gz,reverse
SPT_SH,${RAW_READS_FOLDER}/66_S66_L001_R1_001.fastq.gz,forward
SPT_SH,${RAW_READS_FOLDER}/66_S66_L001_R2_001.fastq.gz,reverse
SPT_SP,${RAW_READS_FOLDER}/67_S67_L001_R1_001.fastq.gz,forward
SPT_SP,${RAW_READS_FOLDER}/67_S67_L001_R2_001.fastq.gz,reverse
BRNVil,${RAW_READS_FOLDER}/6_S6_L001_R1_001.fastq.gz,forward
BRNVil,${RAW_READS_FOLDER}/6_S6_L001_R2_001.fastq.gz,reverse
CH4_RO,${RAW_READS_FOLDER}/70_S70_L001_R1_001.fastq.gz,forward
CH4_RO,${RAW_READS_FOLDER}/70_S70_L001_R2_001.fastq.gz,reverse
CH4_SH,${RAW_READS_FOLDER}/71_S71_L001_R1_001.fastq.gz,forward
CH4_SH,${RAW_READS_FOLDER}/71_S71_L001_R2_001.fastq.gz,reverse
CH4_SP,${RAW_READS_FOLDER}/72_S72_L001_R1_001.fastq.gz,forward
CH4_SP,${RAW_READS_FOLDER}/72_S72_L001_R2_001.fastq.gz,reverse
SE64_74,${RAW_READS_FOLDER}/74_S74_L001_R1_001.fastq.gz,forward
SE64_74,${RAW_READS_FOLDER}/74_S74_L001_R2_001.fastq.gz,reverse
CPCRW_77,${RAW_READS_FOLDER}/77_S77_L001_R1_001.fastq.gz,forward
CPCRW_77,${RAW_READS_FOLDER}/77_S77_L001_R2_001.fastq.gz,reverse
SE16,${RAW_READS_FOLDER}/7_S7_L001_R1_001.fastq.gz,forward
SE16,${RAW_READS_FOLDER}/7_S7_L001_R2_001.fastq.gz,reverse
C_ex_81,${RAW_READS_FOLDER}/81_S81_L001_R1_001.fastq.gz,forward
C_ex_81,${RAW_READS_FOLDER}/81_S81_L001_R2_001.fastq.gz,reverse
BSCMatR3_SH,${RAW_READS_FOLDER}/85_S85_L001_R1_001.fastq.gz,forward
BSCMatR3_SH,${RAW_READS_FOLDER}/85_S85_L001_R2_001.fastq.gz,reverse
BSCMatR3_SP,${RAW_READS_FOLDER}/86_S86_L001_R1_001.fastq.gz,forward
BSCMatR3_SP,${RAW_READS_FOLDER}/86_S86_L001_R2_001.fastq.gz,reverse
SE36,${RAW_READS_FOLDER}/8_S8_L001_R1_001.fastq.gz,forward
SE36,${RAW_READS_FOLDER}/8_S8_L001_R2_001.fastq.gz,reverse
SE64_9,${RAW_READS_FOLDER}/9_S9_L001_R1_001.fastq.gz,forward
SE64_9,${RAW_READS_FOLDER}/9_S9_L001_R2_001.fastq.gz,reverse
A1,${RAW_READS_FOLDER}/A1_S74_L001_R1_001.fastq.gz,forward
A1,${RAW_READS_FOLDER}/A1_S74_L001_R2_001.fastq.gz,reverse
A2,${RAW_READS_FOLDER}/A2_S75_L001_R1_001.fastq.gz,forward
A2,${RAW_READS_FOLDER}/A2_S75_L001_R2_001.fastq.gz,reverse
A3,${RAW_READS_FOLDER}/A3_S76_L001_R1_001.fastq.gz,forward
A3,${RAW_READS_FOLDER}/A3_S76_L001_R2_001.fastq.gz,reverse
R1,${RAW_READS_FOLDER}/R1_S77_L001_R1_001.fastq.gz,forward
R1,${RAW_READS_FOLDER}/R1_S77_L001_R2_001.fastq.gz,reverse
R2,${RAW_READS_FOLDER}/R2_S78_L001_R1_001.fastq.gz,forward
R2,${RAW_READS_FOLDER}/R2_S78_L001_R2_001.fastq.gz,reverse
R3,${RAW_READS_FOLDER}/R3_S79_L001_R1_001.fastq.gz,forward
R3,${RAW_READS_FOLDER}/R3_S79_L001_R2_001.fastq.gz,reverse
Negative_Control,${RAW_READS_FOLDER}/Negative-Control-PGEN_S57_L001_R1_001.fastq.gz,forward
Negative_Control,${RAW_READS_FOLDER}/Negative-Control-PGEN_S57_L001_R2_001.fastq.gz,reverse
RO_Ex,${RAW_READS_FOLDER}/RO-Ex_S49_L001_R1_001.fastq.gz,forward
RO_Ex,${RAW_READS_FOLDER}/RO-Ex_S49_L001_R2_001.fastq.gz,reverse
SP_Ex,${RAW_READS_FOLDER}/SP-Ex_S53_L001_R1_001.fastq.gz,forward
SP_Ex,${RAW_READS_FOLDER}/SP-Ex_S53_L001_R2_001.fastq.gz,reverse
Sto_Bag_2,${RAW_READS_FOLDER}/Sto-bag-2_S55_L001_R1_001.fastq.gz,forward
Sto_Bag_2,${RAW_READS_FOLDER}/Sto-bag-2_S55_L001_R2_001.fastq.gz,reverse
Sto_Bag,${RAW_READS_FOLDER}/Sto-bag_S54_L001_R1_001.fastq.gz,forward
Sto_Bag,${RAW_READS_FOLDER}/Sto-bag_S54_L001_R2_001.fastq.gz,reverse
Tae_RO_2,${RAW_READS_FOLDER}/Tae-RO-2_S35_L001_R1_001.fastq.gz,forward
Tae_RO_2,${RAW_READS_FOLDER}/Tae-RO-2_S35_L001_R2_001.fastq.gz,reverse
Tae_SE_1,${RAW_READS_FOLDER}/Tae-SE-1_S30_L001_R1_001.fastq.gz,forward
Tae_SE_1,${RAW_READS_FOLDER}/Tae-SE-1_S30_L001_R2_001.fastq.gz,reverse
Tae_SE_5,${RAW_READS_FOLDER}/Tae-SE-5_S47_L001_R1_001.fastq.gz,forward
Tae_SE_5,${RAW_READS_FOLDER}/Tae-SE-5_S47_L001_R2_001.fastq.gz,reverse
Tae_SH_2,${RAW_READS_FOLDER}/Tae-SH-2_S33_L001_R1_001.fastq.gz,forward
Tae_SH_2,${RAW_READS_FOLDER}/Tae-SH-2_S33_L001_R2_001.fastq.gz,reverse
Tae_SO_1,${RAW_READS_FOLDER}/Tae-SO-1_S28_L001_R1_001.fastq.gz,forward
Tae_SO_1,${RAW_READS_FOLDER}/Tae-SO-1_S28_L001_R2_001.fastq.gz,reverse
Tae_SO_2,${RAW_READS_FOLDER}/Tae-SO-2_S29_L001_R1_001.fastq.gz,forward
Tae_SO_2,${RAW_READS_FOLDER}/Tae-SO-2_S29_L001_R2_001.fastq.gz,reverse
Tae_SP_1,${RAW_READS_FOLDER}/Tae-SP-1_S26_L001_R1_001.fastq.gz,forward
Tae_SP_1,${RAW_READS_FOLDER}/Tae-SP-1_S26_L001_R2_001.fastq.gz,reverse
Tdu_RO_2,${RAW_READS_FOLDER}/Tdu-RO-2_S14_L001_R1_001.fastq.gz,forward
Tdu_RO_2,${RAW_READS_FOLDER}/Tdu-RO-2_S14_L001_R2_001.fastq.gz,reverse
Tdu_RO_3,${RAW_READS_FOLDER}/Tdu-RO-3_S15_L001_R1_001.fastq.gz,forward
Tdu_RO_3,${RAW_READS_FOLDER}/Tdu-RO-3_S15_L001_R2_001.fastq.gz,reverse
Tdu_SE_1,${RAW_READS_FOLDER}/Tdu-SE-1_S7_L001_R1_001.fastq.gz,forward
Tdu_SE_1,${RAW_READS_FOLDER}/Tdu-SE-1_S7_L001_R2_001.fastq.gz,reverse
Tdu_SE_3,${RAW_READS_FOLDER}/Tdu-SE-3_S9_L001_R1_001.fastq.gz,forward
Tdu_SE_3,${RAW_READS_FOLDER}/Tdu-SE-3_S9_L001_R2_001.fastq.gz,reverse
Tdu_SH_1,${RAW_READS_FOLDER}/Tdu-SH-1_S10_L001_R1_001.fastq.gz,forward
Tdu_SH_1,${RAW_READS_FOLDER}/Tdu-SH-1_S10_L001_R2_001.fastq.gz,reverse
Tdu_SH_2,${RAW_READS_FOLDER}/Tdu-SH-2_S11_L001_R1_001.fastq.gz,forward
Tdu_SH_2,${RAW_READS_FOLDER}/Tdu-SH-2_S11_L001_R2_001.fastq.gz,reverse
Tdu_SO_1,${RAW_READS_FOLDER}/Tdu-SO-1_S4_L001_R1_001.fastq.gz,forward
Tdu_SO_1,${RAW_READS_FOLDER}/Tdu-SO-1_S4_L001_R2_001.fastq.gz,reverse
Tdu_SO_2,${RAW_READS_FOLDER}/Tdu-SO-2_S5_L001_R1_001.fastq.gz,forward
Tdu_SO_2,${RAW_READS_FOLDER}/Tdu-SO-2_S5_L001_R2_001.fastq.gz,reverse
Tdu_SO_3,${RAW_READS_FOLDER}/Tdu-SO-3_S6_L001_R1_001.fastq.gz,forward
Tdu_SO_3,${RAW_READS_FOLDER}/Tdu-SO-3_S6_L001_R2_001.fastq.gz,reverse
Tdu_SP_1,${RAW_READS_FOLDER}/Tdu-SP-1_S1_L001_R1_001.fastq.gz,forward
Tdu_SP_1,${RAW_READS_FOLDER}/Tdu-SP-1_S1_L001_R2_001.fastq.gz,reverse
Tdu_SP_3,${RAW_READS_FOLDER}/Tdu-SP-3_S3_L001_R1_001.fastq.gz,forward
Tdu_SP_3,${RAW_READS_FOLDER}/Tdu-SP-3_S3_L001_R2_001.fastq.gz,reverse
Tpe_RO_1,${RAW_READS_FOLDER}/Tpe-RO-1_S44_L001_R1_001.fastq.gz,forward
Tpe_RO_1,${RAW_READS_FOLDER}/Tpe-RO-1_S44_L001_R2_001.fastq.gz,reverse
Tpe_RO_2,${RAW_READS_FOLDER}/Tpe-RO-2_S45_L001_R1_001.fastq.gz,forward
Tpe_RO_2,${RAW_READS_FOLDER}/Tpe-RO-2_S45_L001_R2_001.fastq.gz,reverse
Tpe_SH_1,${RAW_READS_FOLDER}/Tpe-SH-1_S42_L001_R1_001.fastq.gz,forward
Tpe_SH_1,${RAW_READS_FOLDER}/Tpe-SH-1_S42_L001_R2_001.fastq.gz,reverse
Tpe_SH_2,${RAW_READS_FOLDER}/Tpe-SH-2_S43_L001_R1_001.fastq.gz,forward
Tpe_SH_2,${RAW_READS_FOLDER}/Tpe-SH-2_S43_L001_R2_001.fastq.gz,reverse
Tpe_SO_1,${RAW_READS_FOLDER}/Tpe-SO-1_S38_L001_R1_001.fastq.gz,forward
Tpe_SO_1,${RAW_READS_FOLDER}/Tpe-SO-1_S38_L001_R2_001.fastq.gz,reverse
Tpe_SO_2,${RAW_READS_FOLDER}/Tpe-SO-2_S39_L001_R1_001.fastq.gz,forward
Tpe_SO_2,${RAW_READS_FOLDER}/Tpe-SO-2_S39_L001_R2_001.fastq.gz,reverse
Tpe_SP_1,${RAW_READS_FOLDER}/Tpe-SP-1_S36_L001_R1_001.fastq.gz,forward
Tpe_SP_1,${RAW_READS_FOLDER}/Tpe-SP-1_S36_L001_R2_001.fastq.gz,reverse
Tpe_SP_2,${RAW_READS_FOLDER}/Tpe-SP-2_S37_L001_R1_001.fastq.gz,forward
Tpe_SP_2,${RAW_READS_FOLDER}/Tpe-SP-2_S37_L001_R2_001.fastq.gz,reverse
Ttu_RO_1,${RAW_READS_FOLDER}/Ttu-RO-1_S24_L001_R1_001.fastq.gz,forward
Ttu_RO_1,${RAW_READS_FOLDER}/Ttu-RO-1_S24_L001_R2_001.fastq.gz,reverse
Ttu_RO_2,${RAW_READS_FOLDER}/Ttu-RO-2_S25_L001_R1_001.fastq.gz,forward
Ttu_RO_2,${RAW_READS_FOLDER}/Ttu-RO-2_S25_L001_R2_001.fastq.gz,reverse
Ttu_SE_1,${RAW_READS_FOLDER}/Ttu-SE-1_S20_L001_R1_001.fastq.gz,forward
Ttu_SE_1,${RAW_READS_FOLDER}/Ttu-SE-1_S20_L001_R2_001.fastq.gz,reverse
Ttu_SE_2,${RAW_READS_FOLDER}/Ttu-SE-2_S21_L001_R1_001.fastq.gz,forward
Ttu_SE_2,${RAW_READS_FOLDER}/Ttu-SE-2_S21_L001_R2_001.fastq.gz,reverse
Ttu_SH_1,${RAW_READS_FOLDER}/Ttu-SH-1_S22_L001_R1_001.fastq.gz,forward
Ttu_SH_1,${RAW_READS_FOLDER}/Ttu-SH-1_S22_L001_R2_001.fastq.gz,reverse
Ttu_SH_2,${RAW_READS_FOLDER}/Ttu-SH-2_S23_L001_R1_001.fastq.gz,forward
Ttu_SH_2,${RAW_READS_FOLDER}/Ttu-SH-2_S23_L001_R2_001.fastq.gz,reverse
Ttu_SO_1,${RAW_READS_FOLDER}/Ttu-SO-1_S18_L001_R1_001.fastq.gz,forward
Ttu_SO_1,${RAW_READS_FOLDER}/Ttu-SO-1_S18_L001_R2_001.fastq.gz,reverse
Ttu_SO_2,${RAW_READS_FOLDER}/Ttu-SO-2_S19_L001_R1_001.fastq.gz,forward
Ttu_SO_2,${RAW_READS_FOLDER}/Ttu-SO-2_S19_L001_R2_001.fastq.gz,reverse
Ttu_SP_1,${RAW_READS_FOLDER}/Ttu-SP-1_S16_L001_R1_001.fastq.gz,forward
Ttu_SP_1,${RAW_READS_FOLDER}/Ttu-SP-1_S16_L001_R2_001.fastq.gz,reverse
Ttu_SP_2,${RAW_READS_FOLDER}/Ttu-SP-2_S17_L001_R1_001.fastq.gz,forward
Ttu_SP_2,${RAW_READS_FOLDER}/Ttu-SP-2_S17_L001_R2_001.fastq.gz,reverse
EOF

    echo -e '[X] Creating sample-manifest.csv file\n'
}

create_metadata_file(){
    echo '[ ] Creating samples sample-metadata.tsv file'

    cat > "${CONF_FOLDER}/${METADATA_FILE}" <<EOF
sample-id	SampleName	ForwardPrimer	ReversePrimer	PrimerName	Region	RunDate	PNAs	Beads	Enzyme	SampleType	ControlType	PlantPart	WheatSpecies	PlantPart_WheatSpecies	PloidyLevel	Location	Domestication	Control	Seed	Gnotobiotic	Field
#q2:types	categorical	categorical	categorical	categorical	categorical	categorical	categorical	categorical	categorical	categorical	categorical	categorical	categorical	categorical	categorical	categorical	categorical	categorical	categorical	categorical	categorical
AE_246	Ata_RU	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R4}	No	Yes	Kappa	Seed	NA	Seed	A.tauschii	Seed_A.tauschii	2n	Asia	Ancestral	No	Yes	No	No
AE_498	Ata_UZ	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R4}	No	Yes	Kappa	Seed	NA	Seed	A.tauschii	Seed_A.tauschii	2n	Asia	Ancestral	No	Yes	No	No
BRAJua	Tae_SP5	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R4}	No	Yes	Kappa	Seed	NA	Seed	T.aestivum	Seed_T.aestivum	6n	Europe	Commercial	No	Yes	No	No
BRNVil	Tae_SP4	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R4}	No	Yes	Kappa	Seed	NA	Seed	T.aestivum	Seed_T.aestivum	6n	Europe	Commercial	No	Yes	No	No
BSC2021_27	BSC202127	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R4}	No	Yes	Kappa	Seed	NA	Seed	T.aestivum	Seed_T.aestivum	6n	Europe	Commercial	No	No	Yes	Yes
BSCMatR3_SH	Tae_Sh_C	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R4}	No	Yes	Kappa	Field	NA	Shoot	T.aestivum	Shoot_T.aestivum	6n	America	Commercial	No	No	No	Yes
BSCMatR3_SP	Tae_Sp_C	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R4}	No	Yes	Kappa	Field	NA	Spike	T.aestivum	Spike_T.aestivum	6n	America	Commercial	No	No	No	Yes
CH4_RO	Tae_Ro_A	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R4}	No	Yes	Kappa	Field	NA	Root	T.aestivum	Root_T.aestivum	6n	America	Commercial	No	No	No	Yes
CH4_SH	Tae_Sh_A	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R4}	No	Yes	Kappa	Field	NA	Shoot	T.aestivum	Shoot_T.aestivum	6n	America	Commercial	No	No	No	Yes
CH4_SP	Tae_Sp_A	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R4}	No	Yes	Kappa	Field	NA	Spike	T.aestivum	Spike_T.aestivum	6n	America	Commercial	No	No	No	Yes
CH_RO	Tae_Ro_B	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R4}	No	Yes	Kappa	Field	NA	Root	T.aestivum	Root_T.aestivum	6n	America	Commercial	No	No	No	Yes
CH_SO	Tae_So_B	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R4}	No	Yes	Kappa	Field	NA	Soil	T.aestivum	Soil_T.aestivum	6n	America	Commercial	No	No	No	Yes
CPCRW_77	CPCRW_77	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R4}	No	Yes	Kappa	Control	H2O	NA	NA	NA	NA	NA	NA	Yes	No	No	No
C_ex_30	C_ex_30	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R4}	No	Yes	Kappa	Control	Extraction	NA	NA	NA	NA	NA	NA	Yes	No	No	No
C_ex_81	C_ex_81	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R4}	Yes	Yes	Kappa	Control	Extraction	NA	NA	NA	NA	NA	NA	Yes	No	No	No
SE16	Tpe_SP5	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R4}	No	Yes	Kappa	Seed	NA	Seed	T.spelta	Seed_T.spelta	6n	Europe	Ancestral	No	Yes	No	No
SE1620	Tpe_SP2	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R4}	No	Yes	Kappa	Seed	NA	Seed	T.spelta	Seed_T.spelta	6n	Europe	Ancestral	No	Yes	No	Yes
SE16_RO	Tpe_RO	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R4}	No	Yes	Kappa	Field	NA	Root	T.spelta	Root_T.spelta	6n	Europe	Ancestral	No	No	No	Yes
SE16_SH	Tpe_SH	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R4}	No	Yes	Kappa	Field	NA	Shoot	T.spelta	Shoot_T.spelta	6n	Europe	Ancestral	No	No	No	Yes
SE16_SO	Tpe_SO	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R4}	No	Yes	Kappa	Field	NA	Soil	T.spelta	Soil_T.spelta	6n	Europe	Ancestral	No	No	No	Yes
SE16_SP	Tpe_SP	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R4}	No	Yes	Kappa	Field	NA	Spike	T.spelta	Spike_T.spelta	6n	Europe	Ancestral	No	No	No	Yes
SE36	Tpe_SP3	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R4}	No	Yes	Kappa	Seed	NA	Seed	T.spelta	Seed_T.spelta	6n	Europe	Ancestral	No	Yes	No	No
SE64_74	Tpe_SP4	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R4}	No	Yes	Kappa	Seed	NA	Seed	T.spelta	Seed_T.spelta	6n	Europe	Ancestral	No	Yes	No	Yes
SE64_9	Tpe_SPC	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R4}	No	Yes	Kappa	Seed	NA	Seed	T.spelta	Seed_T.spelta	6n	Europe	Ancestral	No	No	No	Yes
SPB	Tpe_SP1	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R4}	No	Yes	Kappa	Seed	NA	Seed	T.spelta	Seed_T.spelta	6n	Europe	Ancestral	No	Yes	No	No
SPM20Q	Tmo_SP	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R4}	No	Yes	Kappa	Seed	NA	Seed	T.monococcum	Seed_T.monococcum	2n	Europe	Ancestral	No	Yes	No	No
SPT	Ttu_SP	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R4}	No	Yes	Kappa	Seed	NA	Seed	T.turgidum	Seed_T.turgidum	4n	Europe	Ancestral	No	Yes	No	Yes
SPT_RO	SPT_RO	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R4}	No	Yes	Kappa	Field	NA	Root	T.turgidum	Root_T.turgidum	4n	Europe	Ancestral	No	No	No	Yes
SPT_SH	SPT_SH	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R4}	No	Yes	Kappa	Field	NA	Shoot	T.turgidum	Shoot_T.turgidum	4n	Europe	Ancestral	No	No	No	Yes
SPT_SO	SPT_SO	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R4}	No	Yes	Kappa	Field	NA	Soil	T.turgidum	Soil_T.turgidum	4n	Europe	Ancestral	No	No	No	Yes
SPT_SP	SPT_SP	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R4}	No	Yes	Kappa	Field	NA	Spike	T.turgidum	Spike_T.turgidum	4n	Europe	Ancestral	No	No	No	Yes
TBR	Tae_SP1	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R4}	No	Yes	Kappa	Seed	NA	Seed	T.aestivum	Seed_T.aestivum	6n	Europe	Commercial	No	Yes	No	No
TBS	Tae_SP2	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R4}	No	Yes	Kappa	Seed	NA	Seed	T.aestivum	Seed_T.aestivum	6n	Europe	Commercial	No	Yes	Yes	No
TDR	Tdu_SP	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R4}	No	Yes	Kappa	Seed	NA	Seed	T.durum	Seed_T.durum	4n	Europe	Commercial	No	Yes	No	Yes
TDR_RO	TDR_RO	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R4}	No	Yes	Kappa	Field	NA	Root	T.durum	Root_T.durum	4n	Europe	Commercial	No	No	No	Yes
TDR_SH	TDR_SH	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R4}	No	Yes	Kappa	Field	NA	Shoot	T.durum	Shoot_T.durum	4n	Europe	Commercial	No	No	No	Yes
TDR_SP	TDR_SP	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R4}	No	Yes	Kappa	Field	NA	Spike	T.durum	Spike_T.durum	4n	Europe	Commercial	No	No	No	Yes
TKS	Tae_GE	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R4}	No	Yes	Kappa	Seed	NA	Seed	T.aestivum	Seed_T.aestivum	6n	Europe	Commercial	No	Yes	No	No
TRI_10765	Tae_GR	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R4}	No	Yes	Kappa	Seed	NA	Seed	T.aestivum	Seed_T.aestivum	6n	Europe	Commercial	No	Yes	No	No
TRI_1998	Tmo_BU	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R4}	No	Yes	Kappa	Seed	NA	Seed	T.monococcum	Seed_T.monococcum	2n	Europe	Ancestral	No	Yes	No	No
TRI_4101	Tdu_AF	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R4}	No	Yes	Kappa	Seed	NA	Seed	T.durum	Seed_T.durum	4n	Asia	Commercial	No	Yes	No	No
TRI_4168	Tdi_UN	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R4}	No	Yes	Kappa	Seed	NA	Seed	T.dicoccum	Seed_T.dicoccum	4n	Unknown	Ancestral	No	Yes	No	No
TRI_4323	Tmo_UN	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R4}	No	Yes	Kappa	Seed	NA	Seed	T.monococcum	Seed_T.monococcum	2n	Unknown	Ancestral	No	Yes	No	No
TRI_445	Tdi_US	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R4}	No	Yes	Kappa	Seed	NA	Seed	T.dicoccum	Seed_T.dicoccum	4n	America	Ancestral	No	Yes	No	No
TRI_448	Tdu_AR	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R4}	No	Yes	Kappa	Seed	NA	Seed	T.durum	Seed_T.durum	4n	America	Commercial	No	Yes	No	No
TRI_6263	Tdu_IR	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R4}	No	Yes	Kappa	Seed	NA	Seed	T.durum	Seed_T.durum	4n	Asia	Commercial	No	Yes	No	No
TRI_7121	Tae_NZ	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R4}	No	Yes	Kappa	Seed	NA	Seed	T.aestivum	Seed_T.aestivum	6n	Oceania	Commercial	No	Yes	No	No
A1	Tae_Sh1	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R5}	No	Yes	Kappa	Gnotobiotic	NA	Shoot	T.aestivum	Shoot_T.aestivum	6n	Europe	Commercial	No	No	Yes	No
A2	Tae_Sh2	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R5}	No	Yes	Kappa	Gnotobiotic	NA	Shoot	T.aestivum	Shoot_T.aestivum	6n	Europe	Commercial	No	No	Yes	No
A3	Tae_Sh3	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R5}	No	Yes	Kappa	Gnotobiotic	NA	Shoot	T.aestivum	Shoot_T.aestivum	6n	Europe	Commercial	No	No	Yes	No
R1	Tae_Ro1	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R5}	No	Yes	Kappa	Gnotobiotic	NA	Root	T.aestivum	Root_T.aestivum	6n	Europe	Commercial	No	No	Yes	No
R2	Tae_Ro2	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R5}	No	Yes	Kappa	Gnotobiotic	NA	Root	T.aestivum	Root_T.aestivum	6n	Europe	Commercial	No	No	Yes	No
R3	Tae_Ro3	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R5}	No	Yes	Kappa	Gnotobiotic	NA	Root	T.aestivum	Root_T.aestivum	6n	Europe	Commercial	No	No	Yes	No
Negative_Control	Negative_Control	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R6}	No	Yes	Kappa	Control	H2O	NA	NA	NA	NA	NA	NA	Yes	No	No	No
RO_Ex	RO_Ex	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R6}	No	Yes	Kappa	Control	Extraction	NA	NA	NA	NA	NA	NA	Yes	No	No	No
SP_Ex	SP_Ex	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R6}	No	Yes	Kappa	Control	Extraction	NA	NA	NA	NA	NA	NA	Yes	No	No	No
Sto_Bag	Sto_Bag	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R6}	No	Yes	Kappa	Control	Sto_Bag	NA	NA	NA	NA	NA	NA	Yes	No	No	No
Sto_Bag_2	Sto_Bag_2	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R6}	No	Yes	Kappa	Control	Sto_Bag	NA	NA	NA	NA	NA	NA	Yes	No	No	No
Tae_RO_2	Tae_RO_2	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R6}	No	Yes	Kappa	Field	NA	Root	T.aestivum	Root_T.aestivum	6n	America	Commercial	No	No	No	Yes
Tae_SE_1	Tae_SE_1	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R6}	No	Yes	Kappa	Seed	NA	Seed	T.aestivum	Seed_T.aestivum	6n	Europe	Commercial	No	No	Yes	Yes
Tae_SE_5	Tae_SE_5	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R6}	No	Yes	Kappa	Seed	NA	Seed	T.aestivum	Seed_T.aestivum	6n	Europe	Commercial	No	No	No	Yes
Tae_SH_2	Tae_SH_2	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R6}	No	Yes	Kappa	Field	NA	Shoot	T.aestivum	Shoot_T.aestivum	6n	Europe	Commercial	No	No	No	Yes
Tae_SO_1	Tae_SO_1	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R6}	No	Yes	Kappa	Field	NA	Soil	T.aestivum	Soil_T.aestivum	6n	Europe	Commercial	No	No	No	Yes
Tae_SO_2	Tae_SO_2	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R6}	No	Yes	Kappa	Field	NA	Soil	T.aestivum	Soil_T.aestivum	6n	Europe	Commercial	No	No	No	Yes
Tae_SP_1	Tae_SP_1	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R6}	No	Yes	Kappa	Field	NA	Spike	T.aestivum	Spike_T.aestivum	6n	Europe	Commercial	No	No	No	Yes
Tdu_RO_2	Tdu_RO_2	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R6}	No	Yes	Kappa	Field	NA	Root	T.durum	Root_T.durum	4n	Europe	Commercial	No	No	No	Yes
Tdu_RO_3	Tdu_RO_3	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R6}	No	Yes	Kappa	Field	NA	Root	T.durum	Root_T.durum	4n	Europe	Commercial	No	No	No	Yes
Tdu_SE_1	Tdu_SE_1	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R6}	No	Yes	Kappa	Seed	NA	Seed	T.durum	Seed_T.durum	4n	Europe	Commercial	No	No	No	Yes
Tdu_SE_3	Tdu_SE_3	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R6}	No	Yes	Kappa	Seed	NA	Seed	T.durum	Seed_T.durum	4n	Europe	Commercial	No	No	No	Yes
Tdu_SH_1	Tdu_SH_1	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R6}	No	Yes	Kappa	Field	NA	Shoot	T.durum	Shoot_T.durum	4n	Europe	Commercial	No	No	No	Yes
Tdu_SH_2	Tdu_SH_2	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R6}	No	Yes	Kappa	Field	NA	Shoot	T.durum	Shoot_T.durum	4n	Europe	Commercial	No	No	No	Yes
Tdu_SO_1	Tdu_SO_1	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R6}	No	Yes	Kappa	Field	NA	Soil	T.durum	Soil_T.durum	4n	Europe	Commercial	No	No	No	Yes
Tdu_SO_2	Tdu_SO_2	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R6}	No	Yes	Kappa	Field	NA	Soil	T.durum	Soil_T.durum	4n	Europe	Commercial	No	No	No	Yes
Tdu_SO_3	Tdu_SO_3	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R6}	No	Yes	Kappa	Field	NA	Soil	T.durum	Soil_T.durum	4n	Europe	Commercial	No	No	No	Yes
Tdu_SP_1	Tdu_SP_1	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R6}	No	Yes	Kappa	Field	NA	Spike	T.durum	Spike_T.durum	4n	Europe	Commercial	No	No	No	Yes
Tdu_SP_3	Tdu_SP_3	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R6}	No	Yes	Kappa	Field	NA	Spike	T.durum	Spike_T.durum	4n	Europe	Commercial	No	No	No	Yes
Tpe_RO_1	Tpe_RO_1	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R6}	No	Yes	Kappa	Field	NA	Root	T.spelta	Root_T.spelta	6n	Europe	Ancestral	No	No	No	Yes
Tpe_RO_2	Tpe_RO_2	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R6}	No	Yes	Kappa	Field	NA	Root	T.spelta	Root_T.spelta	6n	Europe	Ancestral	No	No	No	Yes
Tpe_SH_1	Tpe_SH_1	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R6}	No	Yes	Kappa	Field	NA	Shoot	T.spelta	Shoot_T.spelta	6n	Europe	Ancestral	No	No	No	Yes
Tpe_SH_2	Tpe_SH_2	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R6}	No	Yes	Kappa	Field	NA	Shoot	T.spelta	Shoot_T.spelta	6n	Europe	Ancestral	No	No	No	Yes
Tpe_SO_1	Tpe_SO_1	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R6}	No	Yes	Kappa	Field	NA	Soil	T.spelta	Soil_T.spelta	6n	Europe	Ancestral	No	No	No	Yes
Tpe_SO_2	Tpe_SO_2	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R6}	No	Yes	Kappa	Field	NA	Soil	T.spelta	Soil_T.spelta	6n	Europe	Ancestral	No	No	No	Yes
Tpe_SP_1	Tpe_SP_1	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R6}	No	Yes	Kappa	Field	NA	Spike	T.spelta	Spike_T.spelta	6n	Europe	Ancestral	No	No	No	Yes
Tpe_SP_2	Tpe_SP_2	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R6}	No	Yes	Kappa	Field	NA	Spike	T.spelta	Spike_T.spelta	6n	Europe	Ancestral	No	No	No	Yes
Ttu_RO_1	Ttu_RO_1	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R6}	No	Yes	Kappa	Field	NA	Root	T.turgidum	Root_T.turgidum	4n	Europe	Ancestral	No	No	No	Yes
Ttu_RO_2	Ttu_RO_2	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R6}	No	Yes	Kappa	Field	NA	Root	T.turgidum	Root_T.turgidum	4n	Europe	Ancestral	No	No	No	Yes
Ttu_SE_1	Ttu_SE_1	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R6}	No	Yes	Kappa	Seed	NA	Seed	T.turgidum	Seed_T.turgidum	4n	Europe	Ancestral	No	No	No	Yes
Ttu_SE_2	Ttu_SE_2	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R6}	No	Yes	Kappa	Seed	NA	Seed	T.turgidum	Seed_T.turgidum	4n	Europe	Ancestral	No	No	No	Yes
Ttu_SH_1	Ttu_SH_1	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R6}	No	Yes	Kappa	Field	NA	Shoot	T.turgidum	Shoot_T.turgidum	4n	Europe	Ancestral	No	No	No	Yes
Ttu_SH_2	Ttu_SH_2	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R6}	No	Yes	Kappa	Field	NA	Shoot	T.turgidum	Shoot_T.turgidum	4n	Europe	Ancestral	No	No	No	Yes
Ttu_SO_1	Ttu_SO_1	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R6}	No	Yes	Kappa	Field	NA	Soil	T.turgidum	Soil_T.turgidum	4n	Europe	Ancestral	No	No	No	Yes
Ttu_SO_2	Ttu_SO_2	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R6}	No	Yes	Kappa	Field	NA	Soil	T.turgidum	Soil_T.turgidum	4n	Europe	Ancestral	No	No	No	Yes
Ttu_SP_1	Ttu_SP_1	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R6}	No	Yes	Kappa	Field	NA	Spike	T.turgidum	Spike_T.turgidum	4n	Europe	Ancestral	No	No	No	Yes
Ttu_SP_2	Ttu_SP_2	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE_R6}	No	Yes	Kappa	Field	NA	Spike	T.turgidum	Spike_T.turgidum	4n	Europe	Ancestral	No	No	No	Yes
EOF

    micromamba run q2-ampl-2025_4 \
        qiime metadata tabulate \
            --m-input-file "${CONF_FOLDER}/${METADATA_FILE}" \
            --o-visualization "${DATA_FOLDER}/sample-metadata.qzv"

    echo -e '[X] Creating samples sample-metadata.tsv file\n'
}

prepare_taxonomic_classifier(){
    echo '[ ] Prepare taxonomic sequence classifier'

    if [[ "${CLASSIFIER_TYPE}" == 'Default' ]] ; then
        echo '    [ ] Downloading precomputed Silva 138.2 default classifier'

        # Silva (16S/18S rRNA) version 138.2
        # ref.: https://library.qiime2.org/data-resources
        # Sklearn Version: 1.4.2
        # Date Trained: 2024-05-30
        SILVA_CLASSIFIER_URI='https://data.qiime2.org/classifiers/sklearn-1.4.2/silva/silva-138-99-nb-classifier.qza'
        SILVA_CLASSIFIER_FILE='silva-138-99-nb-classifier.qza'
        SILVA_CLASSIFIER_SHA256='c08a1aa4d56b449b511f7215543a43249ae9c54b57491428a7e5548a62613616'

        wget -P "${REF_DB_FOLDER}" "${SILVA_CLASSIFIER_URI}"
        if [[ "${SILVA_CLASSIFIER_SHA256}" != `sha256sum "${REF_DB_FOLDER}/${SILVA_CLASSIFIER_FILE}" | cut -d ' ' -f1` ]] ; then
            echo 'Error downloading classifier'
            exit
        fi

        echo -e '    [X] Downloading precomputed Silva 138.2 default classifier\n'
    elif [[ "${CLASSIFIER_TYPE}" == 'Default_weighted' ]] ; then
        echo '    [ ] Downloading precomputed Silva 138.2 default diverse weigthed classifier'

        # Silva (16S/18S rRNA) version 138 | Diverse weighted
        # ref.: https://library.qiime2.org/data-resources
        # Sklearn Version: 1.4.2
        # Date Trained: 2024-07-04
        # Weights created with 14 diverse environments:
        #   “Sediment (non-saline)”, “Plant corpus”, “Animal secretion”, “Sediment (saline)”, “Animal surface”,
        #   “Surface (saline)”, “Plant rhizosphere”, “Soil (non-saline)”, “Animal distal gut”, “Water (saline)”,
        #   “Animal proximal gut”, “Water (non-saline)”, “Animal corpus”, “Plant surface”
        SILVA_CLASSIFIER_URI='https://data.qiime2.org/classifiers/sklearn-1.4.2/silva/silva-138-99-nb-diverse-weighted-classifier.qza'
        SILVA_CLASSIFIER_FILE='silva-138-99-nb-diverse-weighted-classifier.qza'
        SILVA_CLASSIFIER_SHA256='decfae408061fab8ff2fec7dac1fe2a2e0041581589715062cc789bd4f9933db'

        wget -P "${REF_DB_FOLDER}" "${SILVA_CLASSIFIER_URI}"
        if [[ "${SILVA_CLASSIFIER_SHA256}" != `sha256sum "${REF_DB_FOLDER}/${SILVA_CLASSIFIER_FILE}" | cut -d ' ' -f1` ]] ; then
            echo 'Error downloading classifier'
            exit
        fi

        echo -e '    [X] Downloading precomputed Silva 138.2 default diverse weighted classifier\n'
    elif [[ "${CLASSIFIER_TYPE}" == 'Custom_16S' || "${CLASSIFIER_TYPE}" == 'Custom_V4' ]] ; then
        echo '    [ ] Downloading Silva 138.2 reference database'

        # wget -P "${REF_DB_FOLDER}" 'https://www.arb-silva.de/fileadmin/silva_databases/release_138.2/Exports/README.txt'

        wget -P "${REF_DB_FOLDER}" 'https://www.arb-silva.de/fileadmin/silva_databases/release_138.2/Exports/taxonomy/tax_slv_ssu_138.2.txt.gz'
        gunzip "${REF_DB_FOLDER}/tax_slv_ssu_138.2.txt.gz"

        wget -P "${REF_DB_FOLDER}" 'https://www.arb-silva.de/fileadmin/silva_databases/release_138.2/Exports/taxonomy/taxmap_slv_ssu_ref_nr_138.2.txt.gz'
        gunzip "${REF_DB_FOLDER}/taxmap_slv_ssu_ref_nr_138.2.txt.gz"

        wget -P "${REF_DB_FOLDER}" 'https://www.arb-silva.de/fileadmin/silva_databases/release_138.2/Exports/taxonomy/tax_slv_ssu_138.2.tre.gz'
        gunzip "${REF_DB_FOLDER}/tax_slv_ssu_138.2.tre.gz"

        wget -P "${REF_DB_FOLDER}" 'https://www.arb-silva.de/fileadmin/silva_databases/release_138.2/Exports/SILVA_138.2_SSURef_NR99_tax_silva_trunc.fasta.gz'
        gunzip "${REF_DB_FOLDER}/SILVA_138.2_SSURef_NR99_tax_silva_trunc.fasta.gz"
        # wget -P "${REF_DB_FOLDER}" 'https://www.arb-silva.de/fileadmin/silva_databases/release_138.2/Exports/SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz'
        # gunzip "${REF_DB_FOLDER}/SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz"

        micromamba run -p q2-ampl-2025_4 \
            qiime tools import \
                --type 'FeatureData[SILVATaxonomy]' \
                --input-path "${REF_DB_FOLDER}/tax_slv_ssu_138.2.txt" \
                --output-path "${REF_DB_FOLDER}/taxranks-silva-138.2-ssu-nr99.qza"

        micromamba run q2-ampl-2025_4 \
            qiime tools import \
                --type 'FeatureData[SILVATaxidMap]' \
                --input-path "${REF_DB_FOLDER}/taxmap_slv_ssu_ref_nr_138.2.txt" \
                --output-path "${REF_DB_FOLDER}/taxmap-silva-138.2-ssu-nr99.qza"

        micromamba run q2-ampl-2025_4 \
            qiime tools import \
                --type 'Phylogeny[Rooted]' \
                --input-path "${REF_DB_FOLDER}/tax_slv_ssu_138.2.tre" \
                --output-path "${REF_DB_FOLDER}/taxtree-silva-138.2-nr99.qza"

        micromamba run q2-ampl-2025_4 \
            qiime tools import \
                --type 'FeatureData[RNASequence]' \
                --input-path "${REF_DB_FOLDER}/SILVA_138.2_SSURef_NR99_tax_silva_trunc.fasta" \
                --output-path "${REF_DB_FOLDER}/silva-138.2-ssu-nr99-rna-seqs.qza"
        # micromamba run q2-ampl-2025_4 \
        #     qiime tools import \
        #         --type 'FeatureData[RNASequence]' \
        #         --input-path "${REF_DB_FOLDER}/SILVA_138.2_SSURef_NR99_tax_silva.fasta" \
        #         --output-path "${REF_DB_FOLDER}/silva-138.2-ssu-nr99-rna-seqs.qza"

        micromamba run q2-ampl-2025_4 \
            qiime rescript parse-silva-taxonomy \
                --i-taxonomy-tree "${REF_DB_FOLDER}/taxtree-silva-138.2-nr99.qza" \
                --i-taxonomy-map "${REF_DB_FOLDER}/taxmap-silva-138.2-ssu-nr99.qza" \
                --i-taxonomy-ranks "${REF_DB_FOLDER}/taxranks-silva-138.2-ssu-nr99.qza" \
                --o-taxonomy "${REF_DB_FOLDER}/silva-138.2-ssu-nr99-tax.qza"

        # Alternative qiime2 shortcut for all above commands
        # This command gets the non-truncated Silva 138.2 database (bases out of the aligned 16S have not been removed)
        # micromamba run q2-ampl-2025_4 \
        #     qiime rescript get-silva-data \
        #         --p-version '138.2' \
        #         --p-target 'SSURef_NR99' \
        #         --o-silva-sequences "${REF_DB_FOLDER}/silva-138.2-ssu-nr99-rna-seqs.qza" \
        #         --o-silva-taxonomy "${REF_DB_FOLDER}/silva-138.2-ssu-nr99-tax.qza"

        micromamba run q2-ampl-2025_4 \
            qiime rescript reverse-transcribe \
                --i-rna-sequences "${REF_DB_FOLDER}/silva-138.2-ssu-nr99-rna-seqs.qza" \
                --o-dna-sequences "${REF_DB_FOLDER}/silva-138.2-ssu-nr99-seqs.qza"

        echo -e '    [X] Downloading Silva 138.2 reference database \n'

        echo '    [ ] Pre-processing Silva 138.2 reference database'

        # TODO: In case of using the non-truncated database, check if the non-16S bases can be removed at this point
        micromamba run q2-ampl-2025_4 \
            qiime rescript cull-seqs \
                --i-sequences "${REF_DB_FOLDER}/silva-138.2-ssu-nr99-seqs.qza" \
                --o-clean-sequences "${REF_DB_FOLDER}/silva-138.2-ssu-nr99-seqs-cleaned.qza"

        micromamba run q2-ampl-2025_4 \
            qiime rescript filter-seqs-length-by-taxon \
                --i-sequences "${REF_DB_FOLDER}/silva-138.2-ssu-nr99-seqs-cleaned.qza" \
                --i-taxonomy "${REF_DB_FOLDER}/silva-138.2-ssu-nr99-tax.qza" \
                --p-labels Archaea Bacteria Eukaryota \
                --p-min-lens 900 1200 1400 \
                --o-filtered-seqs "${REF_DB_FOLDER}/silva-138.2-ssu-nr99-seqs-filt.qza" \
                --o-discarded-seqs "${REF_DB_FOLDER}/silva-138.2-ssu-nr99-seqs-discard.qza"

        echo -e '    [X] Pre-processing Silva 138.2 reference database\n'

        if [[ "${CLASSIFIER_TYPE}" == 'Custom_16S' ]] ; then
            echo '    [ ] Training 16S Silva 138.2 classifier'

            micromamba run q2-ampl-2025_4 \
                qiime rescript dereplicate \
                    --i-sequences "${REF_DB_FOLDER}/silva-138.2-ssu-nr99-seqs-filt.qza" \
                    --i-taxa "${REF_DB_FOLDER}/silva-138.2-ssu-nr99-tax.qza" \
                    --p-mode "${DEREP_MODE}" \
                    --o-dereplicated-sequences "${REF_DB_FOLDER}/silva-138.2-ssu-nr99-seqs-derep_${DEREP_MODE}.qza" \
                    --o-dereplicated-taxa "${REF_DB_FOLDER}/silva-138.2-ssu-nr99-tax-derep_${DEREP_MODE}.qza"

            micromamba run q2-ampl-2025_4 \
                qiime feature-classifier fit-classifier-naive-bayes \
                    --i-reference-reads "${REF_DB_FOLDER}/silva-138.2-ssu-nr99-seqs-derep_${DEREP_MODE}.qza" \
                    --i-reference-taxonomy "${REF_DB_FOLDER}/silva-138.2-ssu-nr99-tax-derep_${DEREP_MODE}.qza" \
                    --o-classifier "${REF_DB_FOLDER}/silva-138.2-ssu-nr99-derep_${DEREP_MODE}-classifier.qza"

            # TODO: Use these commands instead of the previous one to evaluate the performance of the generated classifier
            # micromamba run q2-ampl-2025_4 \
            #     qiime rescript evaluate-fit-classifier \
            #         --i-sequences "${REF_DB_FOLDER}/silva-138.2-ssu-nr99-seqs-derep_${DEREP_MODE}.qza" \
            #         --i-taxonomy "${REF_DB_FOLDER}/silva-138.2-ssu-nr99-tax-derep_${DEREP_MODE}.qza" \
            #         --o-classifier "${REF_DB_FOLDER}/silva-138.2-ssu-nr99-derep_${DEREP_MODE}-classifier.qza" \
            #         --o-observed-taxonomy "${REF_DB_FOLDER}/silva-138.2-ssu-nr99-derep_${DEREP_MODE}-predicted-taxonomy.qza" \
            #         --o-evaluation "${REF_DB_FOLDER}/silva-138.2-ssu-nr99-derep_${DEREP_MODE}-fit-classifier-evaluation.qzv"
            #
            # micromamba run q2-ampl-2025_4 \
            #     qiime rescript evaluate-taxonomy \
            #         --i-taxonomies "${REF_DB_FOLDER}/silva-138.2-ssu-nr99-tax-derep_${DEREP_MODE}.qza" "${REF_DB_FOLDER}/silva-138.2-ssu-nr99-derep_${DEREP_MODE}-predicted-taxonomy.qza" \
            #         --p-labels "silva-138.2-ssu-nr99-tax-derep_${DEREP_MODE}" "silva-138.2-ssu-nr99-derep_${DEREP_MODE}-predicted-taxonomy" \
            #         --o-taxonomy-stats "${REF_DB_FOLDER}/silva-138.2-ssu-nr99-tax-derep_${DEREP_MODE}-taxonomy-evaluation.qzv"

            echo -e '    [X] Training 16S Silva 138.2 classifier'
        else
            echo '    [ ] Training V4 specific Silva 138.2 classifier'

            micromamba run q2-ampl-2025_4 \
                qiime feature-classifier extract-reads \
                    --i-sequences "${REF_DB_FOLDER}/silva-138.2-ssu-nr99-seqs-filt.qza" \
                    --p-f-primer ${FORWARD_PRIMER} \
                    --p-r-primer ${REVERSE_PRIMER} \
                    --p-n-jobs 2 \
                    --p-read-orientation 'forward' \
                    --o-reads "${REF_DB_FOLDER}/silva-138.2-ssu-nr99-seqs-filt-${PRIMER_NAME}.qza"

            micromamba run q2-ampl-2025_4 \
                qiime rescript dereplicate \
                    --i-sequences "${REF_DB_FOLDER}/silva-138.2-ssu-nr99-seqs-filt-${PRIMER_NAME}.qza" \
                    --i-taxa "${REF_DB_FOLDER}/silva-138.2-ssu-nr99-tax.qza" \
                    --p-mode "${DEREP_MODE}" \
                    --o-dereplicated-sequences "${REF_DB_FOLDER}/silva-138.2-ssu-nr99-seqs-derep_${DEREP_MODE}-${PRIMER_NAME}.qza" \
                    --o-dereplicated-taxa "${REF_DB_FOLDER}/silva-138.2-ssu-nr99-tax-derep_${DEREP_MODE}-${PRIMER_NAME}.qza"

            micromamba run q2-ampl-2025_4 \
                qiime feature-classifier fit-classifier-naive-bayes \
                    --i-reference-reads "${REF_DB_FOLDER}/silva-138.2-ssu-nr99-seqs-derep_${DEREP_MODE}-${PRIMER_NAME}.qza" \
                    --i-reference-taxonomy "${REF_DB_FOLDER}/silva-138.2-ssu-nr99-tax-derep_${DEREP_MODE}-${PRIMER_NAME}.qza" \
                    --o-classifier "${REF_DB_FOLDER}/silva-138.2-ssu-nr99-derep_${DEREP_MODE}-${PRIMER_NAME}-classifier.qza"

            # TODO: Use these commands instead of the previous one to evaluate the performance of the generated classifier
            # micromamba run q2-ampl-2025_4 \
            #     qiime rescript evaluate-fit-classifier \
            #         --i-sequences "${REF_DB_FOLDER}/silva-138.2-ssu-nr99-seqs-derep_${DEREP_MODE}-${PRIMER_NAME}.qza" \
            #         --i-taxonomy "${REF_DB_FOLDER}/silva-138.2-ssu-nr99-tax-derep_${DEREP_MODE}-${PRIMER_NAME}.qza" \
            #         --o-classifier "${REF_DB_FOLDER}/silva-138.2-ssu-nr99-derep_${DEREP_MODE}-${PRIMER_NAME}-classifier.qza" \
            #         --o-observed-taxonomy "${REF_DB_FOLDER}/silva-138.2-ssu-nr99-derep_${DEREP_MODE}-${PRIMER_NAME}-predicted-taxonomy.qza" \
            #         --o-evaluation "${REF_DB_FOLDER}/silva-138.2-ssu-nr99-derep_${DEREP_MODE}-${PRIMER_NAME}-fit-classifier-evaluation.qzv"
            #
            # micromamba run q2-ampl-2025_4 \
            #     qiime rescript evaluate-taxonomy \
            #         --i-taxonomies "${REF_DB_FOLDER}/silva-138.2-ssu-nr99-tax-derep_${DEREP_MODE}-${PRIMER_NAME}.qza" "${REF_DB_FOLDER}/silva-138.2-ssu-nr99-derep_${DEREP_MODE}-${PRIMER_NAME}-predicted-taxonomy.qza" \
            #         --p-labels "silva-138.2-ssu-nr99-tax-derep_${DEREP_MODE}-${PRIMER_NAME}" "silva-138.2-ssu-nr99-derep_${DEREP_MODE}-${PRIMER_NAME}-predicted-taxonomy" \
            #         --o-taxonomy-stats "${REF_DB_FOLDER}/silva-138.2-ssu-nr99-tax-derep_${DEREP_MODE}-${PRIMER_NAME}-taxonomy-evaluation.qzv"

            echo -e '    [X] Training V4 specific Silva 138.2 classifier\n'
        fi
    fi

    echo '[X] Prepare taxonomic sequence classifier'
}

merge_denoised_datasets(){
    echo '[ ] Merging denoised datasets'

    micromamba run q2-ampl-2025_4 \
        qiime feature-table merge-seqs \
            --i-data "${DATA_FOLDER_R4}/asv-seqs.qza" \
            --i-data "${DATA_FOLDER_R5}/asv-seqs.qza" \
            --i-data "${DATA_FOLDER_R6}/asv-seqs.qza" \
            --o-merged-data "${DATA_FOLDER}/asv-seqs.qza"

    micromamba run q2-ampl-2025_4 \
        qiime feature-table merge \
            --i-tables "${DATA_FOLDER_R4}/asv-table.qza" \
            --i-tables "${DATA_FOLDER_R5}/asv-table.qza" \
            --i-tables "${DATA_FOLDER_R6}/asv-table.qza" \
            --o-merged-table "${DATA_FOLDER}/asv-table.qza"

    micromamba run q2-ampl-2025_4 \
        qiime feature-table summarize-plus \
            --i-table "${DATA_FOLDER}/asv-table.qza" \
            --m-metadata-file "${CONF_FOLDER}/${METADATA_FILE}" \
            --o-sample-frequencies "${DATA_FOLDER}/sample-frequencies.qza" \
            --o-feature-frequencies "${DATA_FOLDER}/asv-frequencies.qza" \
            --o-summary "${DATA_FOLDER}/asv-table.qzv"

    micromamba run q2-ampl-2025_4 \
        qiime feature-table tabulate-seqs \
            --i-data "${DATA_FOLDER}/asv-seqs.qza" \
            --m-metadata-file "${DATA_FOLDER}/asv-frequencies.qza" \
            --o-visualization "${DATA_FOLDER}/asv-seqs.qzv"

    echo -e '[X] Merging denoised datasets\n'
}

generate_phylogenetic_tree(){
    echo '[ ] Generating a tree for phylogenetic analyses'

    micromamba run q2-ampl-2025_4 \
        qiime phylogeny align-to-tree-mafft-fasttree \
            --i-sequences "${DATA_FOLDER}/asv-seqs.qza" \
            --p-n-threads ${THREADS} \
            --o-alignment "${DATA_FOLDER}/aligned-asv-seqs.qza" \
            --o-masked-alignment "${DATA_FOLDER}/masked-aligned-asv-seqs.qza" \
            --o-tree "${DATA_FOLDER}/unrooted-tree.qza" \
            --o-rooted-tree "${DATA_FOLDER}/rooted-tree.qza"

    echo -e '[X] Generating a tree for phylogenetic analyses\n'
}

taxonomic_analysis(){
    echo '[ ] Analyzing taxonomic composition'

    if [[ "${CLASSIFIER_TYPE}" == 'Default' ]] ; then
        CLASSIFIER_FILE='silva-138-99-nb-classifier.qza'
    elif [[ "${CLASSIFIER_TYPE}" == 'Default_weighted' ]] ; then
        CLASSIFIER_FILE='silva-138-99-nb-diverse-weighted-classifier.qza'
    elif [[ "${CLASSIFIER_TYPE}" == 'Custom_16S' ]] ; then
        CLASSIFIER_FILE='silva-138.2-ssu-nr99-derep_${DEREP_MODE}-classifier.qza'
    elif [[ "${CLASSIFIER_TYPE}" == 'Custom_V4' ]] ; then
        CLASSIFIER_FILE='silva-138.2-ssu-nr99-derep_${DEREP_MODE}-${PRIMER_NAME}-classifier.qza'
    fi

    micromamba run q2-ampl-2025_4 \
        qiime feature-classifier classify-sklearn \
            --i-classifier "${REF_DB_FOLDER}/${CLASSIFIER_FILE}" \
            --i-reads "${DATA_FOLDER}/asv-seqs.qza" \
            --p-n-jobs "${THREADS_CLASSIFIER}" \
            --o-classification "${DATA_FOLDER}/taxonomy.qza"

    micromamba run q2-ampl-2025_4 \
        qiime metadata tabulate \
            --m-input-file "${DATA_FOLDER}/taxonomy.qza" \
            --m-input-file "${DATA_FOLDER}/asv-seqs.qza" \
            --o-visualization "${DATA_FOLDER}/taxonomy-with-seqs.qzv"

    micromamba run q2-ampl-2025_4 \
        qiime taxa barplot \
            --i-table "${DATA_FOLDER}/asv-table.qza" \
            --i-taxonomy "${DATA_FOLDER}/taxonomy.qza" \
            --m-metadata-file "${CONF_FOLDER}/${METADATA_FILE}" \
            --o-visualization "${DATA_FOLDER}/taxa-bar-plots.qzv"

    echo -e '[X] Analyzing taxonomic composition\n'
}

export_qiime_artifacts(){
    echo '[ ] Exporting denoised Qiime artifacts'

    cp "${CONF_FOLDER}/{METADATA_FILE}" "${BIOM_FOLDER}"

    micromamba run q2-ampl-2025_4 \
        qiime tools export \
            --input-path "${DATA_FOLDER}/asv-table.qza" \
            --output-path "${BIOM_FOLDER}"
    mv "${BIOM_FOLDER}/feature-table.biom" "${BIOM_FOLDER}/asv-table.biom"

    micromamba run q2-ampl-2025_4 \
        biom convert --to-tsv \
            --input-fp "${BIOM_FOLDER}/asv-table.biom" \
            --output-fp "${BIOM_FOLDER}/asv-table.tsv"
    sed -i '1d' "${BIOM_FOLDER}/asv-table.tsv" # Delete comment line: '# Constructed from biom file'

    TREE_FILE=`unzip -l "${DATA_FOLDER}/unrooted-tree.qza" | grep 'tree.nwk' | sed -E 's/ +/\t/g' | cut -f5`
    unzip -p "${DATA_FOLDER}/unrooted-tree.qza" "${TREE_FILE}" > "${BIOM_FOLDER}/unrooted-tree.nwk"

    TREE_FILE=`unzip -l "${DATA_FOLDER}/rooted-tree.qza" | grep 'tree.nwk' | sed -E 's/ +/\t/g' | cut -f5`
    unzip -p "${DATA_FOLDER}/rooted-tree.qza" "${TREE_FILE}" > "${BIOM_FOLDER}/rooted-tree.nwk"

    micromamba run q2-ampl-2025_4 \
        qiime tools export \
            --input-path "${DATA_FOLDER}/taxonomy.qza" \
            --output-path "${BIOM_FOLDER}"

    echo -e '[X] Exporting denoised Qiime artifacts\n'
}

filter_datasets(){
    echo '[ ] Filtering datasets'

    micromamba run q2-ampl-2025_4 \
        qiime taxa filter-seqs \
            --i-sequences "${DATA_FOLDER}/asv-seqs.qza" \
            --i-taxonomy "${DATA_FOLDER}/taxonomy.qza" \
            --p-include 'p__' \
            --p-exclude 'p__;,Chloroplast,Mitochondria' \
            --o-filtered-sequences "${DATA_FOLDER}/asv-seqs-filt.qza"

    micromamba run q2-ampl-2025_4 \
        qiime taxa filter-table \
            --i-table "${DATA_FOLDER}/asv-table.qza" \
            --i-taxonomy "${DATA_FOLDER}/taxonomy.qza" \
            --p-include 'p__' \
            --p-exclude 'p__;,Chloroplast,Mitochondria' \
            --o-filtered-table "${DATA_FOLDER}/asv-table-filt.qza"

    micromamba run q2-ampl-2025_4 \
        qiime feature-table summarize-plus \
            --i-table "${DATA_FOLDER}/asv-table-filt.qza" \
            --m-metadata-file "${CONF_FOLDER}/${METADATA_FILE}" \
            --o-sample-frequencies "${DATA_FOLDER}/sample-frequencies-filt.qza" \
            --o-feature-frequencies "${DATA_FOLDER}/asv-frequencies-filt.qza" \
            --o-summary "${DATA_FOLDER}/asv-table-filt.qzv"

    micromamba run q2-ampl-2025_4 \
        qiime feature-table tabulate-seqs \
            --i-data "${DATA_FOLDER}/asv-seqs-filt.qza" \
            --m-metadata-file "${DATA_FOLDER}/asv-frequencies-filt.qza" \
            --o-visualization "${DATA_FOLDER}/asv-seqs-filt.qzv"

    micromamba run q2-ampl-2025_4 \
        qiime phylogeny align-to-tree-mafft-fasttree \
            --i-sequences "${DATA_FOLDER}/asv-seqs-filt.qza" \
            --p-n-threads ${THREADS} \
            --o-alignment "${DATA_FOLDER}/aligned-asv-seqs-filt.qza" \
            --o-masked-alignment "${DATA_FOLDER}/masked-aligned-asv-seqs-filt.qza" \
            --o-tree "${DATA_FOLDER}/unrooted-tree-filt.qza" \
            --o-rooted-tree "${DATA_FOLDER}/rooted-tree-filt.qza"

    micromamba run q2-ampl-2025_4 \
        qiime rescript filter-taxa \
            --i-taxonomy "${DATA_FOLDER}/taxonomy.qza" \
            --m-ids-to-keep-file "${DATA_FOLDER}/asv-frequencies-filt.qza" \
            --o-filtered-taxonomy "${DATA_FOLDER}/taxonomy-filt.qza"

    micromamba run q2-ampl-2025_4 \
        qiime metadata tabulate \
            --m-input-file "${DATA_FOLDER}/taxonomy-filt.qza" \
            --m-input-file "${DATA_FOLDER}/asv-seqs-filt.qza" \
            --o-visualization "${DATA_FOLDER}/taxonomy-with-seqs-filt.qzv"

    micromamba run q2-ampl-2025_4 \
        qiime taxa barplot \
            --i-table "${DATA_FOLDER}/asv-table-filt.qza" \
            --i-taxonomy "${DATA_FOLDER}/taxonomy-filt.qza" \
            --m-metadata-file "${CONF_FOLDER}/${METADATA_FILE}" \
            --o-visualization "${DATA_FOLDER}/taxa-bar-plots-filt.qzv"

    echo -e '[X] Filtering datasets\n'
}

export_qiime_artifacts_filt(){
    echo '[ ] Exporting filtered Qiime artifacts'

    micromamba run q2-ampl-2025_4 \
        qiime tools export \
            --input-path "${DATA_FOLDER}/asv-table-filt.qza" \
            --output-path "${BIOM_FOLDER}"
    mv "${BIOM_FOLDER}/feature-table.biom" "${BIOM_FOLDER}/asv-table-filt.biom"

    micromamba run q2-ampl-2025_4 \
        biom convert --to-tsv \
            --input-fp "${BIOM_FOLDER}/asv-table-filt.biom" \
            --output-fp "${BIOM_FOLDER}/asv-table-filt.tsv"
    sed -i '1d' "${BIOM_FOLDER}/asv-table-filt.tsv" # Delete comment line: '# Constructed from biom file'

    TREE_FILE=`unzip -l "${DATA_FOLDER}/unrooted-tree-filt.qza" | grep 'tree.nwk' | sed -E 's/ +/\t/g' | cut -f5`
    unzip -p "${DATA_FOLDER}/unrooted-tree-filt.qza" "${TREE_FILE}" > "${BIOM_FOLDER}/unrooted-tree-filt.nwk"

    TREE_FILE=`unzip -l "${DATA_FOLDER}/rooted-tree-filt.qza" | grep 'tree.nwk' | sed -E 's/ +/\t/g' | cut -f5`
    unzip -p "${DATA_FOLDER}/rooted-tree-filt.qza" "${TREE_FILE}" > "${BIOM_FOLDER}/rooted-tree-filt.nwk"

    mv "${BIOM_FOLDER}/taxonomy.tsv" "${BIOM_FOLDER}/taxonomy-full.tsv"
    micromamba run q2-ampl-2025_4 \
        qiime tools export \
            --input-path "${DATA_FOLDER}/taxonomy-filt.qza" \
            --output-path "${BIOM_FOLDER}"
    mv "${BIOM_FOLDER}/taxonomy.tsv" "${BIOM_FOLDER}/taxonomy-filt.tsv"
    mv "${BIOM_FOLDER}/taxonomy-full.tsv" "${BIOM_FOLDER}/taxonomy.tsv"

    echo -e '[X] Exporting filtered Qiime artifacts\n'
}

decontaminate_datasets(){
    echo '[ ] Decontaminating datasets'

    micromamba run q2-ampl-2025_4 \
        qiime quality-control decontam-identify \
            --i-table "${DATA_FOLDER}/asv-table-filt.qza" \
            --m-metadata-file "${CONF_FOLDER}/${METADATA_FILE}" \
            --p-method 'prevalence' \
            --p-prev-control-column 'Control' \
            --p-prev-control-indicator 'Yes' \
            --o-decontam-scores "${DATA_FOLDER}/decontam-scores.qza"

    micromamba run q2-ampl-2025_4 \
        qiime quality-control decontam-score-viz \
            --i-decontam-scores "${DATA_FOLDER}/decontam-scores.qza" \
            --i-table "${DATA_FOLDER}/asv-table-filt.qza" \
            --i-rep-seqs "${DATA_FOLDER}/asv-seqs-filt.qza" \
            --p-threshold 0.1 \
            --p-weighted \
            --p-bin-size 0.05 \
            --o-visualization "${DATA_FOLDER}/decontam-scores.qzv"

    micromamba run q2-ampl-2025_4 \
        qiime quality-control decontam-remove \
            --i-decontam-scores "${DATA_FOLDER}/decontam-scores.qza" \
            --i-table "${DATA_FOLDER}/asv-table-filt.qza" \
            --i-rep-seqs "${DATA_FOLDER}/asv-seqs-filt.qza" \
            --p-threshold 0.1 \
            --o-filtered-table "${DATA_FOLDER}/asv-table-decont.qza" \
            --o-filtered-rep-seqs "${DATA_FOLDER}/asv-seqs-decont.qza"

    micromamba run q2-ampl-2025_4 \
        qiime feature-table summarize-plus \
            --i-table "${DATA_FOLDER}/asv-table-decont.qza" \
            --m-metadata-file "${CONF_FOLDER}/${METADATA_FILE}" \
            --o-sample-frequencies "${DATA_FOLDER}/sample-frequencies-decont.qza" \
            --o-feature-frequencies "${DATA_FOLDER}/asv-frequencies-decont.qza" \
            --o-summary "${DATA_FOLDER}/asv-table-decont.qzv"

    micromamba run q2-ampl-2025_4 \
        qiime feature-table tabulate-seqs \
            --i-data "${DATA_FOLDER}/asv-seqs-decont.qza" \
            --m-metadata-file "${DATA_FOLDER}/asv-frequencies-decont.qza" \
            --o-visualization "${DATA_FOLDER}/asv-seqs-decont.qzv"

    micromamba run q2-ampl-2025_4 \
        qiime phylogeny align-to-tree-mafft-fasttree \
            --i-sequences "${DATA_FOLDER}/asv-seqs-decont.qza" \
            --p-n-threads ${THREADS} \
            --o-alignment "${DATA_FOLDER}/aligned-asv-seqs-decont.qza" \
            --o-masked-alignment "${DATA_FOLDER}/masked-aligned-asv-seqs-decont.qza" \
            --o-tree "${DATA_FOLDER}/unrooted-tree-decont.qza" \
            --o-rooted-tree "${DATA_FOLDER}/rooted-tree-decont.qza"

    micromamba run q2-ampl-2025_4 \
        qiime rescript filter-taxa \
            --i-taxonomy "${DATA_FOLDER}/taxonomy.qza" \
            --m-ids-to-keep-file "${DATA_FOLDER}/asv-frequencies-decont.qza" \
            --o-filtered-taxonomy "${DATA_FOLDER}/taxonomy-decont.qza"

    micromamba run q2-ampl-2025_4 \
        qiime metadata tabulate \
            --m-input-file "${DATA_FOLDER}/taxonomy-decont.qza" \
            --m-input-file "${DATA_FOLDER}/asv-seqs-decont.qza" \
            --o-visualization "${DATA_FOLDER}/taxonomy-with-seqs-decont.qzv"

    micromamba run q2-ampl-2025_4 \
        qiime taxa barplot \
            --i-table "${DATA_FOLDER}/asv-table-decont.qza" \
            --i-taxonomy "${DATA_FOLDER}/taxonomy-decont.qza" \
            --m-metadata-file "${CONF_FOLDER}/${METADATA_FILE}" \
            --o-visualization "${DATA_FOLDER}/taxa-bar-plots-decont.qzv"

    echo -e '[X] Decontaminating datasets\n'
}

export_qiime_artifacts_decont(){
    echo '[ ] Exporting decontaminated Qiime artifacts'

    micromamba run q2-ampl-2025_4 \
        qiime tools export \
            --input-path "${DATA_FOLDER}/asv-table-decont.qza" \
            --output-path "${BIOM_FOLDER}"
    mv "${BIOM_FOLDER}/feature-table.biom" "${BIOM_FOLDER}/asv-table-decont.biom"

    micromamba run q2-ampl-2025_4 \
        biom convert --to-tsv \
            --input-fp "${BIOM_FOLDER}/asv-table-decont.biom" \
            --output-fp "${BIOM_FOLDER}/asv-table-decont.tsv"
    sed -i '1d' "${BIOM_FOLDER}/asv-table-decont.tsv" # Delete comment line: '# Constructed from biom file'

    TREE_FILE=`unzip -l "${DATA_FOLDER}/unrooted-tree-decont.qza" | grep 'tree.nwk' | sed -E 's/ +/\t/g' | cut -f5`
    unzip -p "${DATA_FOLDER}/unrooted-tree-decont.qza" "${TREE_FILE}" > "${BIOM_FOLDER}/unrooted-tree-decont.nwk"

    TREE_FILE=`unzip -l "${DATA_FOLDER}/rooted-tree-decont.qza" | grep 'tree.nwk' | sed -E 's/ +/\t/g' | cut -f5`
    unzip -p "${DATA_FOLDER}/rooted-tree-decont.qza" "${TREE_FILE}" > "${BIOM_FOLDER}/rooted-tree-decont.nwk"

    mv "${BIOM_FOLDER}/taxonomy.tsv" "${BIOM_FOLDER}/taxonomy-full.tsv"
    micromamba run q2-ampl-2025_4 \
        qiime tools export \
            --input-path "${DATA_FOLDER}/taxonomy-decont.qza" \
            --output-path "${BIOM_FOLDER}"
    mv "${BIOM_FOLDER}/taxonomy.tsv" "${BIOM_FOLDER}/taxonomy-decont.tsv"
    mv "${BIOM_FOLDER}/taxonomy-full.tsv" "${BIOM_FOLDER}/taxonomy.tsv"

    echo -e '[X] Exporting decontaminated Qiime artifacts\n'
}

define_sample_groups(){
    echo '[ ] Creating groups of samples'

    # Seed samples grouping
    micromamba run q2-ampl-2025_4 \
        qiime feature-table filter-samples \
            --i-table "${DATA_FOLDER}/asv-table-decont.qza" \
            --m-metadata-file "${CONF_FOLDER}/${METADATA_FILE}" \
            --p-where '[Seed]="Yes"' \
            --o-filtered-table "${DATA_FOLDER}/asv-table-seed.qza"

    micromamba run q2-ampl-2025_4 \
        qiime taxa collapse \
            --i-table "${DATA_FOLDER}/asv-table-seed.qza" \
            --i-taxonomy "${DATA_FOLDER}/taxonomy.qza" \
            --p-level 6 \
            --o-collapsed-table "${DATA_FOLDER}/genus-table-seed.qza"

    micromamba run q2-ampl-2025_4 \
        qiime tools export \
            --input-path "${DATA_FOLDER}/asv-table-seed.qza" \
            --output-path "${BIOM_FOLDER}"
    mv "${BIOM_FOLDER}/feature-table.biom" "${BIOM_FOLDER}/asv-table-seed.biom"

    micromamba run q2-ampl-2025_4 \
        biom convert --to-tsv \
            --input-fp "${BIOM_FOLDER}/asv-table-seed.biom" \
            --output-fp "${BIOM_FOLDER}/asv-table-seed.tsv"
    sed -i '1d' "${BIOM_FOLDER}/asv-table-seed.tsv" # Delete comment line: '# Constructed from biom file'

    micromamba run q2-ampl-2025_4 \
        qiime feature-table summarize-plus \
            --i-table "${DATA_FOLDER}/asv-table-seed.qza" \
            --m-metadata-file "${CONF_FOLDER}/${METADATA_FILE}" \
            --o-sample-frequencies "${DATA_FOLDER}/sample-frequencies-seed.qza" \
            --o-feature-frequencies "${DATA_FOLDER}/asv-frequencies-seed.qza" \
            --o-summary "${DATA_FOLDER}/asv-table-seed.qzv"

    micromamba run q2-ampl-2025_4 \
        qiime rescript filter-taxa \
            --i-taxonomy "${DATA_FOLDER}/taxonomy.qza" \
            --m-ids-to-keep-file "${DATA_FOLDER}/asv-frequencies-seed.qza" \
            --o-filtered-taxonomy "${DATA_FOLDER}/taxonomy-seed.qza"

    mv "${BIOM_FOLDER}/taxonomy.tsv" "${BIOM_FOLDER}/taxonomy-full.tsv"
    micromamba run q2-ampl-2025_4 \
        qiime tools export \
            --input-path "${DATA_FOLDER}/taxonomy-seed.qza" \
            --output-path "${BIOM_FOLDER}"
    mv "${BIOM_FOLDER}/taxonomy.tsv" "${BIOM_FOLDER}/taxonomy-seed.tsv"
    mv "${BIOM_FOLDER}/taxonomy-full.tsv" "${BIOM_FOLDER}/taxonomy.tsv"

    micromamba run q2-ampl-2025_4 \
        qiime taxa barplot \
            --i-table "${DATA_FOLDER}/asv-table-seed.qza" \
            --i-taxonomy "${DATA_FOLDER}/taxonomy.qza" \
            --m-metadata-file "${CONF_FOLDER}/${METADATA_FILE}" \
            --o-visualization "${DATA_FOLDER}/taxa-bar-plots-seed.qzv"

    # Gnotobiotic samples grouping
    micromamba run q2-ampl-2025_4 \
        qiime feature-table filter-samples \
            --i-table "${DATA_FOLDER}/asv-table-decont.qza" \
            --m-metadata-file "${CONF_FOLDER}/${METADATA_FILE}" \
            --p-where '[Gnotobiotic]="Yes"' \
            --o-filtered-table "${DATA_FOLDER}/asv-table-gnotobiotic.qza"

    micromamba run q2-ampl-2025_4 \
        qiime taxa collapse \
            --i-table "${DATA_FOLDER}/asv-table-gnotobiotic.qza" \
            --i-taxonomy "${DATA_FOLDER}/taxonomy.qza" \
            --p-level 6 \
            --o-collapsed-table "${DATA_FOLDER}/genus-table-gnotobiotic.qza"

    micromamba run q2-ampl-2025_4 \
        qiime tools export \
            --input-path "${DATA_FOLDER}/asv-table-gnotobiotic.qza" \
            --output-path "${BIOM_FOLDER}"
    mv "${BIOM_FOLDER}/feature-table.biom" "${BIOM_FOLDER}/asv-table-gnotobiotic.biom"

    micromamba run q2-ampl-2025_4 \
        biom convert --to-tsv \
            --input-fp "${BIOM_FOLDER}/asv-table-gnotobiotic.biom" \
            --output-fp "${BIOM_FOLDER}/asv-table-gnotobiotic.tsv"
    sed -i '1d' "${BIOM_FOLDER}/asv-table-gnotobiotic.tsv" # Delete comment line: '# Constructed from biom file'

    micromamba run q2-ampl-2025_4 \
        qiime feature-table summarize-plus \
            --i-table "${DATA_FOLDER}/asv-table-gnotobiotic.qza" \
            --m-metadata-file "${CONF_FOLDER}/${METADATA_FILE}" \
            --o-sample-frequencies "${DATA_FOLDER}/sample-frequencies-gnotobiotic.qza" \
            --o-feature-frequencies "${DATA_FOLDER}/asv-frequencies-gnotobiotic.qza" \
            --o-summary "${DATA_FOLDER}/asv-table-gnotobiotic.qzv"

    micromamba run q2-ampl-2025_4 \
        qiime rescript filter-taxa \
            --i-taxonomy "${DATA_FOLDER}/taxonomy.qza" \
            --m-ids-to-keep-file "${DATA_FOLDER}/asv-frequencies-gnotobiotic.qza" \
            --o-filtered-taxonomy "${DATA_FOLDER}/taxonomy-gnotobiotic.qza"

    mv "${BIOM_FOLDER}/taxonomy.tsv" "${BIOM_FOLDER}/taxonomy-full.tsv"
    micromamba run q2-ampl-2025_4 \
        qiime tools export \
            --input-path "${DATA_FOLDER}/taxonomy-gnotobiotic.qza" \
            --output-path "${BIOM_FOLDER}"
    mv "${BIOM_FOLDER}/taxonomy.tsv" "${BIOM_FOLDER}/taxonomy-gnotobiotic.tsv"
    mv "${BIOM_FOLDER}/taxonomy-full.tsv" "${BIOM_FOLDER}/taxonomy.tsv"

    micromamba run q2-ampl-2025_4 \
        qiime taxa barplot \
            --i-table "${DATA_FOLDER}/asv-table-gnotobiotic.qza" \
            --i-taxonomy "${DATA_FOLDER}/taxonomy.qza" \
            --m-metadata-file "${CONF_FOLDER}/${METADATA_FILE}" \
            --o-visualization "${DATA_FOLDER}/taxa-bar-plots-gnotobiotic.qzv"

    # Field samples grouping
    micromamba run q2-ampl-2025_4 \
        qiime feature-table filter-samples \
            --i-table "${DATA_FOLDER}/asv-table-decont.qza" \
            --m-metadata-file "${CONF_FOLDER}/${METADATA_FILE}" \
            --p-where '[Field]="Yes"' \
            --o-filtered-table "${DATA_FOLDER}/asv-table-field.qza"

    micromamba run q2-ampl-2025_4 \
        qiime taxa collapse \
            --i-table "${DATA_FOLDER}/asv-table-field.qza" \
            --i-taxonomy "${DATA_FOLDER}/taxonomy.qza" \
            --p-level 6 \
            --o-collapsed-table "${DATA_FOLDER}/genus-table-field.qza"

    micromamba run q2-ampl-2025_4 \
        qiime tools export \
            --input-path "${DATA_FOLDER}/asv-table-field.qza" \
            --output-path "${BIOM_FOLDER}"
    mv "${BIOM_FOLDER}/feature-table.biom" "${BIOM_FOLDER}/asv-table-field.biom"

    micromamba run q2-ampl-2025_4 \
        biom convert --to-tsv \
            --input-fp "${BIOM_FOLDER}/asv-table-field.biom" \
            --output-fp "${BIOM_FOLDER}/asv-table-field.tsv"
    sed -i '1d' "${BIOM_FOLDER}/asv-table-field.tsv" # Delete comment line: '# Constructed from biom file'

    micromamba run q2-ampl-2025_4 \
        qiime feature-table summarize-plus \
            --i-table "${DATA_FOLDER}/asv-table-field.qza" \
            --m-metadata-file "${CONF_FOLDER}/${METADATA_FILE}" \
            --o-sample-frequencies "${DATA_FOLDER}/sample-frequencies-field.qza" \
            --o-feature-frequencies "${DATA_FOLDER}/asv-frequencies-field.qza" \
            --o-summary "${DATA_FOLDER}/asv-table-field.qzv"

    micromamba run q2-ampl-2025_4 \
        qiime rescript filter-taxa \
            --i-taxonomy "${DATA_FOLDER}/taxonomy.qza" \
            --m-ids-to-keep-file "${DATA_FOLDER}/asv-frequencies-field.qza" \
            --o-filtered-taxonomy "${DATA_FOLDER}/taxonomy-field.qza"

    mv "${BIOM_FOLDER}/taxonomy.tsv" "${BIOM_FOLDER}/taxonomy-full.tsv"
    micromamba run q2-ampl-2025_4 \
        qiime tools export \
            --input-path "${DATA_FOLDER}/taxonomy-field.qza" \
            --output-path "${BIOM_FOLDER}"
    mv "${BIOM_FOLDER}/taxonomy.tsv" "${BIOM_FOLDER}/taxonomy-field.tsv"
    mv "${BIOM_FOLDER}/taxonomy-full.tsv" "${BIOM_FOLDER}/taxonomy.tsv"

    micromamba run q2-ampl-2025_4 \
        qiime taxa barplot \
            --i-table "${DATA_FOLDER}/asv-table-field.qza" \
            --i-taxonomy "${DATA_FOLDER}/taxonomy.qza" \
            --m-metadata-file "${CONF_FOLDER}/${METADATA_FILE}" \
            --o-visualization "${DATA_FOLDER}/taxa-bar-plots-field.qzv"

    echo -e '[X] Creating groups of samples\n'
}

ancombc2_differential_abundance_analysis_seed(){
    echo '[ ] Analyzing seed differential abundance with ANCOM-BC2'

    METRICS=( 'lfc' 'se' 'W' 'p' 'q' 'diff' 'passed_ss' )

    # Check effect of PloidyLevel using '2n' as reference
    micromamba run q2-ampl-2025_4 \
        qiime composition ancombc2 \
            --i-table "${DATA_FOLDER}/asv-table-seed.qza" \
            --m-metadata-file "${CONF_FOLDER}/${METADATA_FILE}" \
            --p-fixed-effects-formula 'PloidyLevel' \
            --p-reference-levels 'PloidyLevel::2n' \
            --o-ancombc2-output "${DATA_FOLDER}/ancombc2-seed-ploidylevel@2n-asv.qza"

    mkdir -p "${BIOM_FOLDER}/tmp"
    micromamba run q2-ampl-2025_4 \
        qiime tools export \
            --input-path "${DATA_FOLDER}/ancombc2-seed-ploidylevel@2n-asv.qza" \
            --output-path "${BIOM_FOLDER}/tmp"
    for m in ${METRICS[@]} ; do
        sed '1d' "${BIOM_FOLDER}/tmp/${m}.jsonl" | \
            jq '.taxon,."(Intercept)",."PloidyLevel::4n",."PloidyLevel::6n"' | \
            sed -e 's/^"//' -e 's/"$//' -e 's/null/NaN/' | paste - - - - | \
            sed "1iTaxon\t(Intercept)|${m}\tPloidyLevel::4n|${m}\tPloidyLevel::6n|${m}" > "${BIOM_FOLDER}/tmp/${m}.tsv"
    done
    csvjoin -t -I -y 0 -c 'Taxon' "${BIOM_FOLDER}/tmp/lfc.tsv" "${BIOM_FOLDER}/tmp/se.tsv" "${BIOM_FOLDER}/tmp/W.tsv" "${BIOM_FOLDER}/tmp/p.tsv" "${BIOM_FOLDER}/tmp/q.tsv" "${BIOM_FOLDER}/tmp/diff.tsv" "${BIOM_FOLDER}/tmp/passed_ss.tsv" | \
        csvformat -T > "${BIOM_FOLDER}/ancombc2-seed-ploidylevel@2n-asv.tsv"
    rm -r "${BIOM_FOLDER}/tmp"

    micromamba run q2-ampl-2025_4 \
        qiime composition tabulate \
            --i-data "${DATA_FOLDER}/ancombc2-seed-ploidylevel@2n-asv.qza" \
            --o-visualization "${DATA_FOLDER}/ancombc2-seed-ploidylevel@2n-asv.qzv"

    micromamba run q2-ampl-2025_4 \
        qiime composition ancombc2-visualizer \
            --i-data "${DATA_FOLDER}/ancombc2-seed-ploidylevel@2n-asv.qza" \
            --i-taxonomy "${DATA_FOLDER}/taxonomy.qza" \
            --o-visualization "${DATA_FOLDER}/ancombc2-barplot-seed-ploidylevel@2n-asv.qzv"

    micromamba run q2-ampl-2025_4 \
        qiime composition ancombc2 \
            --i-table "${DATA_FOLDER}/genus-table-seed.qza" \
            --m-metadata-file "${CONF_FOLDER}/${METADATA_FILE}" \
            --p-fixed-effects-formula 'PloidyLevel' \
            --p-reference-levels 'PloidyLevel::2n' \
            --o-ancombc2-output "${DATA_FOLDER}/ancombc2-seed-ploidylevel@2n-genus.qza"

    mkdir -p "${BIOM_FOLDER}/tmp"
    micromamba run q2-ampl-2025_4 \
        qiime tools export \
            --input-path "${DATA_FOLDER}/ancombc2-seed-ploidylevel@2n-genus.qza" \
            --output-path "${BIOM_FOLDER}/tmp"
    for m in ${METRICS[@]} ; do
        sed '1d' "${BIOM_FOLDER}/tmp/${m}.jsonl" | \
            jq '.taxon,."(Intercept)",."PloidyLevel::4n",."PloidyLevel::6n"' | \
            sed -e 's/^"//' -e 's/"$//' -e 's/null/NaN/' | paste - - - - | \
            sed "1iTaxon\t(Intercept)|${m}\tPloidyLevel::4n|${m}\tPloidyLevel::6n|${m}" > "${BIOM_FOLDER}/tmp/${m}.tsv"
    done
    csvjoin -t -I -y 0 -c 'Taxon' "${BIOM_FOLDER}/tmp/lfc.tsv" "${BIOM_FOLDER}/tmp/se.tsv" "${BIOM_FOLDER}/tmp/W.tsv" "${BIOM_FOLDER}/tmp/p.tsv" "${BIOM_FOLDER}/tmp/q.tsv" "${BIOM_FOLDER}/tmp/diff.tsv" "${BIOM_FOLDER}/tmp/passed_ss.tsv" | \
        csvformat -T > "${BIOM_FOLDER}/ancombc2-seed-ploidylevel@2n-genus.tsv"
    rm -r "${BIOM_FOLDER}/tmp"

    micromamba run q2-ampl-2025_4 \
        qiime composition tabulate \
            --i-data "${DATA_FOLDER}/ancombc2-seed-ploidylevel@2n-genus.qza" \
            --o-visualization "${DATA_FOLDER}/ancombc2-seed-ploidylevel@2n-genus.qzv"

    micromamba run q2-ampl-2025_4 \
        qiime composition ancombc2-visualizer \
            --i-data "${DATA_FOLDER}/ancombc2-seed-ploidylevel@2n-genus.qza" \
            --o-visualization "${DATA_FOLDER}/ancombc2-barplot-seed-ploidylevel@2n-genus.qzv"

    # Check effect of PloidyLevel using '4n' as reference
    micromamba run q2-ampl-2025_4 \
        qiime composition ancombc2 \
            --i-table "${DATA_FOLDER}/asv-table-seed.qza" \
            --m-metadata-file "${CONF_FOLDER}/${METADATA_FILE}" \
            --p-fixed-effects-formula 'PloidyLevel' \
            --p-reference-levels 'PloidyLevel::4n' \
            --o-ancombc2-output "${DATA_FOLDER}/ancombc2-seed-ploidylevel@4n-asv.qza"

    mkdir -p "${BIOM_FOLDER}/tmp"
    micromamba run q2-ampl-2025_4 \
        qiime tools export \
            --input-path "${DATA_FOLDER}/ancombc2-seed-ploidylevel@4n-asv.qza" \
            --output-path "${BIOM_FOLDER}/tmp"
    for m in ${METRICS[@]} ; do
        sed '1d' "${BIOM_FOLDER}/tmp/${m}.jsonl" | \
            jq '.taxon,."(Intercept)",."PloidyLevel::2n",."PloidyLevel::6n"' | \
            sed -e 's/^"//' -e 's/"$//' -e 's/null/NaN/' | paste - - - - | \
            sed "1iTaxon\t(Intercept)|${m}\tPloidyLevel::2n|${m}\tPloidyLevel::6n|${m}" > "${BIOM_FOLDER}/tmp/${m}.tsv"
    done
    csvjoin -t -I -y 0 -c 'Taxon' "${BIOM_FOLDER}/tmp/lfc.tsv" "${BIOM_FOLDER}/tmp/se.tsv" "${BIOM_FOLDER}/tmp/W.tsv" "${BIOM_FOLDER}/tmp/p.tsv" "${BIOM_FOLDER}/tmp/q.tsv" "${BIOM_FOLDER}/tmp/diff.tsv" "${BIOM_FOLDER}/tmp/passed_ss.tsv" | \
        csvformat -T > "${BIOM_FOLDER}/ancombc2-seed-ploidylevel@4n-asv.tsv"
    rm -r "${BIOM_FOLDER}/tmp"

    micromamba run q2-ampl-2025_4 \
        qiime composition tabulate \
            --i-data "${DATA_FOLDER}/ancombc2-seed-ploidylevel@4n-asv.qza" \
            --o-visualization "${DATA_FOLDER}/ancombc2-seed-ploidylevel@4n-asv.qzv"

    micromamba run q2-ampl-2025_4 \
        qiime composition ancombc2-visualizer \
            --i-data "${DATA_FOLDER}/ancombc2-seed-ploidylevel@4n-asv.qza" \
            --i-taxonomy "${DATA_FOLDER}/taxonomy.qza" \
            --o-visualization "${DATA_FOLDER}/ancombc2-barplot-seed-ploidylevel@4n-asv.qzv"

    micromamba run q2-ampl-2025_4 \
        qiime composition ancombc2 \
            --i-table "${DATA_FOLDER}/genus-table-seed.qza" \
            --m-metadata-file "${CONF_FOLDER}/${METADATA_FILE}" \
            --p-fixed-effects-formula 'PloidyLevel' \
            --p-reference-levels 'PloidyLevel::4n' \
            --o-ancombc2-output "${DATA_FOLDER}/ancombc2-seed-ploidylevel@4n-genus.qza"

    mkdir -p "${BIOM_FOLDER}/tmp"
    micromamba run q2-ampl-2025_4 \
        qiime tools export \
            --input-path "${DATA_FOLDER}/ancombc2-seed-ploidylevel@4n-genus.qza" \
            --output-path "${BIOM_FOLDER}/tmp"
    for m in ${METRICS[@]} ; do
        sed '1d' "${BIOM_FOLDER}/tmp/${m}.jsonl" | \
            jq '.taxon,."(Intercept)",."PloidyLevel::2n",."PloidyLevel::6n"' | \
            sed -e 's/^"//' -e 's/"$//' -e 's/null/NaN/' | paste - - - - | \
            sed "1iTaxon\t(Intercept)|${m}\tPloidyLevel::2n|${m}\tPloidyLevel::6n|${m}" > "${BIOM_FOLDER}/tmp/${m}.tsv"
    done
    csvjoin -t -I -y 0 -c 'Taxon' "${BIOM_FOLDER}/tmp/lfc.tsv" "${BIOM_FOLDER}/tmp/se.tsv" "${BIOM_FOLDER}/tmp/W.tsv" "${BIOM_FOLDER}/tmp/p.tsv" "${BIOM_FOLDER}/tmp/q.tsv" "${BIOM_FOLDER}/tmp/diff.tsv" "${BIOM_FOLDER}/tmp/passed_ss.tsv" | \
        csvformat -T > "${BIOM_FOLDER}/ancombc2-seed-ploidylevel@4n-genus.tsv"
    rm -r "${BIOM_FOLDER}/tmp"

    micromamba run q2-ampl-2025_4 \
        qiime composition tabulate \
            --i-data "${DATA_FOLDER}/ancombc2-seed-ploidylevel@4n-genus.qza" \
            --o-visualization "${DATA_FOLDER}/ancombc2-seed-ploidylevel@4n-genus.qzv"

    micromamba run q2-ampl-2025_4 \
        qiime composition ancombc2-visualizer \
            --i-data "${DATA_FOLDER}/ancombc2-seed-ploidylevel@4n-genus.qza" \
            --o-visualization "${DATA_FOLDER}/ancombc2-barplot-seed-ploidylevel@4n-genus.qzv"

    # Check effect of PloidyLevel using '6n' as reference
    micromamba run q2-ampl-2025_4 \
        qiime composition ancombc2 \
            --i-table "${DATA_FOLDER}/asv-table-seed.qza" \
            --m-metadata-file "${CONF_FOLDER}/${METADATA_FILE}" \
            --p-fixed-effects-formula 'PloidyLevel' \
            --p-reference-levels 'PloidyLevel::6n' \
            --o-ancombc2-output "${DATA_FOLDER}/ancombc2-seed-ploidylevel@6n-asv.qza"

    mkdir -p "${BIOM_FOLDER}/tmp"
    micromamba run q2-ampl-2025_4 \
        qiime tools export \
            --input-path "${DATA_FOLDER}/ancombc2-seed-ploidylevel@6n-asv.qza" \
            --output-path "${BIOM_FOLDER}/tmp"
    for m in ${METRICS[@]} ; do
        sed '1d' "${BIOM_FOLDER}/tmp/${m}.jsonl" | \
            jq '.taxon,."(Intercept)",."PloidyLevel::2n",."PloidyLevel::4n"' | \
            sed -e 's/^"//' -e 's/"$//' -e 's/null/NaN/' | paste - - - - | \
            sed "1iTaxon\t(Intercept)|${m}\tPloidyLevel::2n|${m}\tPloidyLevel::4n|${m}" > "${BIOM_FOLDER}/tmp/${m}.tsv"
    done
    csvjoin -t -I -y 0 -c 'Taxon' "${BIOM_FOLDER}/tmp/lfc.tsv" "${BIOM_FOLDER}/tmp/se.tsv" "${BIOM_FOLDER}/tmp/W.tsv" "${BIOM_FOLDER}/tmp/p.tsv" "${BIOM_FOLDER}/tmp/q.tsv" "${BIOM_FOLDER}/tmp/diff.tsv" "${BIOM_FOLDER}/tmp/passed_ss.tsv" | \
        csvformat -T > "${BIOM_FOLDER}/ancombc2-seed-ploidylevel@6n-asv.tsv"
    rm -r "${BIOM_FOLDER}/tmp"

    micromamba run q2-ampl-2025_4 \
        qiime composition tabulate \
            --i-data "${DATA_FOLDER}/ancombc2-seed-ploidylevel@6n-asv.qza" \
            --o-visualization "${DATA_FOLDER}/ancombc2-seed-ploidylevel@6n-asv.qzv"

    micromamba run q2-ampl-2025_4 \
        qiime composition ancombc2-visualizer \
            --i-data "${DATA_FOLDER}/ancombc2-seed-ploidylevel@6n-asv.qza" \
            --i-taxonomy "${DATA_FOLDER}/taxonomy.qza" \
            --o-visualization "${DATA_FOLDER}/ancombc2-barplot-seed-ploidylevel@6n-asv.qzv"

    micromamba run q2-ampl-2025_4 \
        qiime composition ancombc2 \
            --i-table "${DATA_FOLDER}/genus-table-seed.qza" \
            --m-metadata-file "${CONF_FOLDER}/${METADATA_FILE}" \
            --p-fixed-effects-formula 'PloidyLevel' \
            --p-reference-levels 'PloidyLevel::6n' \
            --o-ancombc2-output "${DATA_FOLDER}/ancombc2-seed-ploidylevel@6n-genus.qza"

    mkdir -p "${BIOM_FOLDER}/tmp"
    micromamba run q2-ampl-2025_4 \
        qiime tools export \
            --input-path "${DATA_FOLDER}/ancombc2-seed-ploidylevel@6n-genus.qza" \
            --output-path "${BIOM_FOLDER}/tmp"
    for m in ${METRICS[@]} ; do
        sed '1d' "${BIOM_FOLDER}/tmp/${m}.jsonl" | \
            jq '.taxon,."(Intercept)",."PloidyLevel::2n",."PloidyLevel::4n"' | \
            sed -e 's/^"//' -e 's/"$//' -e 's/null/NaN/' | paste - - - - | \
            sed "1iTaxon\t(Intercept)|${m}\tPloidyLevel::2n|${m}\tPloidyLevel::4n|${m}" > "${BIOM_FOLDER}/tmp/${m}.tsv"
    done
    csvjoin -t -I -y 0 -c 'Taxon' "${BIOM_FOLDER}/tmp/lfc.tsv" "${BIOM_FOLDER}/tmp/se.tsv" "${BIOM_FOLDER}/tmp/W.tsv" "${BIOM_FOLDER}/tmp/p.tsv" "${BIOM_FOLDER}/tmp/q.tsv" "${BIOM_FOLDER}/tmp/diff.tsv" "${BIOM_FOLDER}/tmp/passed_ss.tsv" | \
        csvformat -T > "${BIOM_FOLDER}/ancombc2-seed-ploidylevel@6n-genus.tsv"
    rm -r "${BIOM_FOLDER}/tmp"

    micromamba run q2-ampl-2025_4 \
        qiime composition tabulate \
            --i-data "${DATA_FOLDER}/ancombc2-seed-ploidylevel@6n-genus.qza" \
            --o-visualization "${DATA_FOLDER}/ancombc2-seed-ploidylevel@6n-genus.qzv"

    micromamba run q2-ampl-2025_4 \
        qiime composition ancombc2-visualizer \
            --i-data "${DATA_FOLDER}/ancombc2-seed-ploidylevel@6n-genus.qza" \
            --o-visualization "${DATA_FOLDER}/ancombc2-barplot-seed-ploidylevel@6n-genus.qzv"

    # Check effect of Domestication using 'Ancestral' as reference
    micromamba run q2-ampl-2025_4 \
        qiime composition ancombc2 \
            --i-table "${DATA_FOLDER}/asv-table-seed.qza" \
            --m-metadata-file "${CONF_FOLDER}/${METADATA_FILE}" \
            --p-fixed-effects-formula 'Domestication' \
            --p-reference-levels 'Domestication::Ancestral' \
            --o-ancombc2-output "${DATA_FOLDER}/ancombc2-seed-domestication@ancestral-asv.qza"

    mkdir -p "${BIOM_FOLDER}/tmp"
    micromamba run q2-ampl-2025_4 \
        qiime tools export \
            --input-path "${DATA_FOLDER}/ancombc2-seed-domestication@ancestral-asv.qza" \
            --output-path "${BIOM_FOLDER}/tmp"
    for m in ${METRICS[@]} ; do
        sed '1d' "${BIOM_FOLDER}/tmp/${m}.jsonl" | \
            jq '.taxon,."(Intercept)",."Domestication::Commercial"' | \
            sed -e 's/^"//' -e 's/"$//' -e 's/null/NaN/' | paste - - - | \
            sed "1iTaxon\t(Intercept)|${m}\tDomestication::Commercial|${m}" > "${BIOM_FOLDER}/tmp/${m}.tsv"
    done
    csvjoin -t -I -y 0 -c 'Taxon' "${BIOM_FOLDER}/tmp/lfc.tsv" "${BIOM_FOLDER}/tmp/se.tsv" "${BIOM_FOLDER}/tmp/W.tsv" "${BIOM_FOLDER}/tmp/p.tsv" "${BIOM_FOLDER}/tmp/q.tsv" "${BIOM_FOLDER}/tmp/diff.tsv" "${BIOM_FOLDER}/tmp/passed_ss.tsv" | \
        csvformat -T > "${BIOM_FOLDER}/ancombc2-seed-domestication@ancestral-asv.tsv"
    rm -r "${BIOM_FOLDER}/tmp"

    micromamba run q2-ampl-2025_4 \
        qiime composition tabulate \
            --i-data "${DATA_FOLDER}/ancombc2-seed-domestication@ancestral-asv.qza" \
            --o-visualization "${DATA_FOLDER}/ancombc2-seed-domestication@ancestral-asv.qzv"

    micromamba run q2-ampl-2025_4 \
        qiime composition ancombc2-visualizer \
            --i-data "${DATA_FOLDER}/ancombc2-seed-domestication@ancestral-asv.qza" \
            --i-taxonomy "${DATA_FOLDER}/taxonomy.qza" \
            --o-visualization "${DATA_FOLDER}/ancombc2-barplot-seed-domestication@ancestral-asv.qzv"

    micromamba run q2-ampl-2025_4 \
        qiime composition ancombc2 \
            --i-table "${DATA_FOLDER}/genus-table-seed.qza" \
            --m-metadata-file "${CONF_FOLDER}/${METADATA_FILE}" \
            --p-fixed-effects-formula 'Domestication' \
            --p-reference-levels 'Domestication::Ancestral' \
            --o-ancombc2-output "${DATA_FOLDER}/ancombc2-seed-domestication@ancestral-genus.qza"

    mkdir -p "${BIOM_FOLDER}/tmp"
    micromamba run q2-ampl-2025_4 \
        qiime tools export \
            --input-path "${DATA_FOLDER}/ancombc2-seed-domestication@ancestral-genus.qza" \
            --output-path "${BIOM_FOLDER}/tmp"
    for m in ${METRICS[@]} ; do
        sed '1d' "${BIOM_FOLDER}/tmp/${m}.jsonl" | \
            jq '.taxon,."(Intercept)",."Domestication::Commercial"' | \
            sed -e 's/^"//' -e 's/"$//' -e 's/null/NaN/' | paste - - - | \
            sed "1iTaxon\t(Intercept)|${m}\tDomestication::Commercial|${m}" > "${BIOM_FOLDER}/tmp/${m}.tsv"
    done
    csvjoin -t -I -y 0 -c 'Taxon' "${BIOM_FOLDER}/tmp/lfc.tsv" "${BIOM_FOLDER}/tmp/se.tsv" "${BIOM_FOLDER}/tmp/W.tsv" "${BIOM_FOLDER}/tmp/p.tsv" "${BIOM_FOLDER}/tmp/q.tsv" "${BIOM_FOLDER}/tmp/diff.tsv" "${BIOM_FOLDER}/tmp/passed_ss.tsv" | \
        csvformat -T > "${BIOM_FOLDER}/ancombc2-seed-domestication@ancestral-genus.tsv"
    rm -r "${BIOM_FOLDER}/tmp"

    micromamba run q2-ampl-2025_4 \
        qiime composition tabulate \
            --i-data "${DATA_FOLDER}/ancombc2-seed-domestication@ancestral-genus.qza" \
            --o-visualization "${DATA_FOLDER}/ancombc2-seed-domestication@ancestral-genus.qzv"

    micromamba run q2-ampl-2025_4 \
        qiime composition ancombc2-visualizer \
            --i-data "${DATA_FOLDER}/ancombc2-seed-domestication@ancestral-genus.qza" \
            --o-visualization "${DATA_FOLDER}/ancombc2-barplot-seed-domestication@ancestral-genus.qzv"

    # Check effect of Location using 'Europe' as reference
    micromamba run q2-ampl-2025_4 \
        qiime composition ancombc2 \
            --i-table "${DATA_FOLDER}/asv-table-seed.qza" \
            --m-metadata-file "${CONF_FOLDER}/${METADATA_FILE}" \
            --p-fixed-effects-formula 'Location' \
            --p-reference-levels 'Location::Europe' \
            --o-ancombc2-output "${DATA_FOLDER}/ancombc2-seed-location@europe-asv.qza"

    mkdir -p "${BIOM_FOLDER}/tmp"
    micromamba run q2-ampl-2025_4 \
        qiime tools export \
            --input-path "${DATA_FOLDER}/ancombc2-seed-location@europe-asv.qza" \
            --output-path "${BIOM_FOLDER}/tmp"
    for m in ${METRICS[@]} ; do
        sed '1d' "${BIOM_FOLDER}/tmp/${m}.jsonl" | \
            jq '.taxon,."(Intercept)",."Location::America",."Location::Asia",."Location::Oceania",."Location::Unknown"' | \
            sed -e 's/^"//' -e 's/"$//' -e 's/null/NaN/' | paste - - - - - - | \
            sed "1iTaxon\t(Intercept)|${m}\tLocation::America|${m}\tLocation::Asia|${m}\tLocation::Oceania|${m}\tLocation::Unknown|${m}" > "${BIOM_FOLDER}/tmp/${m}.tsv"
    done
    csvjoin -t -I -y 0 -c 'Taxon' "${BIOM_FOLDER}/tmp/lfc.tsv" "${BIOM_FOLDER}/tmp/se.tsv" "${BIOM_FOLDER}/tmp/W.tsv" "${BIOM_FOLDER}/tmp/p.tsv" "${BIOM_FOLDER}/tmp/q.tsv" "${BIOM_FOLDER}/tmp/diff.tsv" "${BIOM_FOLDER}/tmp/passed_ss.tsv" | \
        csvformat -T > "${BIOM_FOLDER}/ancombc2-seed-location@europe-asv.tsv"
    rm -r "${BIOM_FOLDER}/tmp"

    micromamba run q2-ampl-2025_4 \
        qiime composition tabulate \
            --i-data "${DATA_FOLDER}/ancombc2-seed-location@europe-asv.qza" \
            --o-visualization "${DATA_FOLDER}/ancombc2-seed-location@europe-asv.qzv"

    micromamba run q2-ampl-2025_4 \
        qiime composition ancombc2-visualizer \
            --i-data "${DATA_FOLDER}/ancombc2-seed-location@europe-asv.qza" \
            --i-taxonomy "${DATA_FOLDER}/taxonomy.qza" \
            --o-visualization "${DATA_FOLDER}/ancombc2-barplot-seed-location@europe-asv.qzv"

    micromamba run q2-ampl-2025_4 \
        qiime composition ancombc2 \
            --i-table "${DATA_FOLDER}/genus-table-seed.qza" \
            --m-metadata-file "${CONF_FOLDER}/${METADATA_FILE}" \
            --p-fixed-effects-formula 'Location' \
            --p-reference-levels 'Location::Europe' \
            --o-ancombc2-output "${DATA_FOLDER}/ancombc2-seed-location@europe-genus.qza"

    mkdir -p "${BIOM_FOLDER}/tmp"
    micromamba run q2-ampl-2025_4 \
        qiime tools export \
            --input-path "${DATA_FOLDER}/ancombc2-seed-location@europe-genus.qza" \
            --output-path "${BIOM_FOLDER}/tmp"
    for m in ${METRICS[@]} ; do
        sed '1d' "${BIOM_FOLDER}/tmp/${m}.jsonl" | \
            jq '.taxon,."(Intercept)",."Location::America",."Location::Asia",."Location::Unknown"' | \
            sed -e 's/^"//' -e 's/"$//' -e 's/null/NaN/' | paste - - - - - | \
            sed "1iTaxon\t(Intercept)|${m}\tLocation::America|${m}\tLocation::Asia|${m}\tLocation::Unknown|${m}" > "${BIOM_FOLDER}/tmp/${m}.tsv"
    done
    csvjoin -t -I -y 0 -c 'Taxon' "${BIOM_FOLDER}/tmp/lfc.tsv" "${BIOM_FOLDER}/tmp/se.tsv" "${BIOM_FOLDER}/tmp/W.tsv" "${BIOM_FOLDER}/tmp/p.tsv" "${BIOM_FOLDER}/tmp/q.tsv" "${BIOM_FOLDER}/tmp/diff.tsv" "${BIOM_FOLDER}/tmp/passed_ss.tsv" | \
        csvformat -T > "${BIOM_FOLDER}/ancombc2-seed-location@europe-genus.tsv"
    rm -r "${BIOM_FOLDER}/tmp"

    micromamba run q2-ampl-2025_4 \
        qiime composition tabulate \
            --i-data "${DATA_FOLDER}/ancombc2-seed-location@europe-genus.qza" \
            --o-visualization "${DATA_FOLDER}/ancombc2-seed-location@europe-genus.qzv"

    micromamba run q2-ampl-2025_4 \
        qiime composition ancombc2-visualizer \
            --i-data "${DATA_FOLDER}/ancombc2-seed-location@europe-genus.qza" \
            --o-visualization "${DATA_FOLDER}/ancombc2-barplot-seed-location@europe-genus.qzv"

    echo -e '[X] Analyzing seed differential abundance with ANCOM-BC2\n'
}

ancombc2_differential_abundance_analysis_gnotobiotic(){
    echo '[ ] Analyzing gnotobiotic differential abundance with ANCOM-BC2'

    METRICS=( 'lfc' 'se' 'W' 'p' 'q' 'diff' 'passed_ss' )

    # Check effect of PlantPart using 'Seed' as reference
    micromamba run q2-ampl-2025_4 \
        qiime composition ancombc2 \
            --i-table "${DATA_FOLDER}/asv-table-gnotobiotic.qza" \
            --m-metadata-file "${CONF_FOLDER}/${METADATA_FILE}" \
            --p-fixed-effects-formula 'PlantPart' \
            --p-reference-levels 'PlantPart::Seed' \
            --o-ancombc2-output "${DATA_FOLDER}/ancombc2-gnotobiotic-plantpart@seed-asv.qza"

    mkdir -p "${BIOM_FOLDER}/tmp"
    micromamba run q2-ampl-2025_4 \
        qiime tools export \
            --input-path "${DATA_FOLDER}/ancombc2-gnotobiotic-plantpart@seed-asv.qza" \
            --output-path "${BIOM_FOLDER}/tmp"
    for m in ${METRICS[@]} ; do
        sed '1d' "${BIOM_FOLDER}/tmp/${m}.jsonl" | \
            jq '.taxon,."(Intercept)",."PlantPart::Root",."PlantPart::Shoot"' | \
            sed -e 's/^"//' -e 's/"$//' -e 's/null/NaN/' | paste - - - - | \
            sed "1iTaxon\t(Intercept)|${m}\tPlantPart::Root|${m}\tPlantPart::Shoot|${m}" > "${BIOM_FOLDER}/tmp/${m}.tsv"
    done
    csvjoin -t -I -y 0 -c 'Taxon' "${BIOM_FOLDER}/tmp/lfc.tsv" "${BIOM_FOLDER}/tmp/se.tsv" "${BIOM_FOLDER}/tmp/W.tsv" "${BIOM_FOLDER}/tmp/p.tsv" "${BIOM_FOLDER}/tmp/q.tsv" "${BIOM_FOLDER}/tmp/diff.tsv" "${BIOM_FOLDER}/tmp/passed_ss.tsv" | \
        csvformat -T > "${BIOM_FOLDER}/ancombc2-gnotobiotic-plantpart@seed-asv.tsv"
    rm -r "${BIOM_FOLDER}/tmp"

    micromamba run q2-ampl-2025_4 \
        qiime composition tabulate \
            --i-data "${DATA_FOLDER}/ancombc2-gnotobiotic-plantpart@seed-asv.qza" \
            --o-visualization "${DATA_FOLDER}/ancombc2-gnotobiotic-plantpart@seed-asv.qzv"

    micromamba run q2-ampl-2025_4 \
        qiime composition ancombc2-visualizer \
            --i-data "${DATA_FOLDER}/ancombc2-gnotobiotic-plantpart@seed-asv.qza" \
            --i-taxonomy "${DATA_FOLDER}/taxonomy.qza" \
            --o-visualization "${DATA_FOLDER}/ancombc2-barplot-gnotobiotic-plantpart@seed-asv.qzv"

    micromamba run q2-ampl-2025_4 \
        qiime composition ancombc2 \
            --i-table "${DATA_FOLDER}/genus-table-gnotobiotic.qza" \
            --m-metadata-file "${CONF_FOLDER}/${METADATA_FILE}" \
            --p-fixed-effects-formula 'PlantPart' \
            --p-reference-levels 'PlantPart::Seed' \
            --o-ancombc2-output "${DATA_FOLDER}/ancombc2-gnotobiotic-plantpart@seed-genus.qza"

    mkdir -p "${BIOM_FOLDER}/tmp"
    micromamba run q2-ampl-2025_4 \
        qiime tools export \
            --input-path "${DATA_FOLDER}/ancombc2-gnotobiotic-plantpart@seed-genus.qza" \
            --output-path "${BIOM_FOLDER}/tmp"
    for m in ${METRICS[@]} ; do
        sed '1d' "${BIOM_FOLDER}/tmp/${m}.jsonl" | \
            jq '.taxon,."(Intercept)",."PlantPart::Root",."PlantPart::Shoot"' | \
            sed -e 's/^"//' -e 's/"$//' -e 's/null/NaN/' | paste - - - - | \
            sed "1iTaxon\t(Intercept)|${m}\tPlantPart::Root|${m}\tPlantPart::Shoot|${m}" > "${BIOM_FOLDER}/tmp/${m}.tsv"
    done
    csvjoin -t -I -y 0 -c 'Taxon' "${BIOM_FOLDER}/tmp/lfc.tsv" "${BIOM_FOLDER}/tmp/se.tsv" "${BIOM_FOLDER}/tmp/W.tsv" "${BIOM_FOLDER}/tmp/p.tsv" "${BIOM_FOLDER}/tmp/q.tsv" "${BIOM_FOLDER}/tmp/diff.tsv" "${BIOM_FOLDER}/tmp/passed_ss.tsv" | \
        csvformat -T > "${BIOM_FOLDER}/ancombc2-gnotobiotic-plantpart@seed-genus.tsv"
    rm -r "${BIOM_FOLDER}/tmp"

    micromamba run q2-ampl-2025_4 \
        qiime composition tabulate \
            --i-data "${DATA_FOLDER}/ancombc2-gnotobiotic-plantpart@seed-genus.qza" \
            --o-visualization "${DATA_FOLDER}/ancombc2-gnotobiotic-plantpart@seed-genus.qzv"

    micromamba run q2-ampl-2025_4 \
        qiime composition ancombc2-visualizer \
            --i-data "${DATA_FOLDER}/ancombc2-gnotobiotic-plantpart@seed-genus.qza" \
            --o-visualization "${DATA_FOLDER}/ancombc2-barplot-gnotobiotic-plantpart@seed-genus.qzv"

    # Check effect of PlantPart using 'Shoot' as reference
    micromamba run q2-ampl-2025_4 \
        qiime composition ancombc2 \
            --i-table "${DATA_FOLDER}/asv-table-gnotobiotic.qza" \
            --m-metadata-file "${CONF_FOLDER}/${METADATA_FILE}" \
            --p-fixed-effects-formula 'PlantPart' \
            --p-reference-levels 'PlantPart::Shoot' \
            --o-ancombc2-output "${DATA_FOLDER}/ancombc2-gnotobiotic-plantpart@shoot-asv.qza"

    mkdir -p "${BIOM_FOLDER}/tmp"
    micromamba run q2-ampl-2025_4 \
        qiime tools export \
            --input-path "${DATA_FOLDER}/ancombc2-gnotobiotic-plantpart@shoot-asv.qza" \
            --output-path "${BIOM_FOLDER}/tmp"
    for m in ${METRICS[@]} ; do
        sed '1d' "${BIOM_FOLDER}/tmp/${m}.jsonl" | \
            jq '.taxon,."(Intercept)",."PlantPart::Root",."PlantPart::Seed"' | \
            sed -e 's/^"//' -e 's/"$//' -e 's/null/NaN/' | paste - - - - | \
            sed "1iTaxon\t(Intercept)|${m}\tPlantPart::Root|${m}\tPlantPart::Seed|${m}" > "${BIOM_FOLDER}/tmp/${m}.tsv"
    done
    csvjoin -t -I -y 0 -c 'Taxon' "${BIOM_FOLDER}/tmp/lfc.tsv" "${BIOM_FOLDER}/tmp/se.tsv" "${BIOM_FOLDER}/tmp/W.tsv" "${BIOM_FOLDER}/tmp/p.tsv" "${BIOM_FOLDER}/tmp/q.tsv" "${BIOM_FOLDER}/tmp/diff.tsv" "${BIOM_FOLDER}/tmp/passed_ss.tsv" | \
        csvformat -T > "${BIOM_FOLDER}/ancombc2-gnotobiotic-plantpart@shoot-asv.tsv"
    rm -r "${BIOM_FOLDER}/tmp"

    micromamba run q2-ampl-2025_4 \
        qiime composition tabulate \
            --i-data "${DATA_FOLDER}/ancombc2-gnotobiotic-plantpart@shoot-asv.qza" \
            --o-visualization "${DATA_FOLDER}/ancombc2-gnotobiotic-plantpart@shoot-asv.qzv"

    micromamba run q2-ampl-2025_4 \
        qiime composition ancombc2-visualizer \
            --i-data "${DATA_FOLDER}/ancombc2-gnotobiotic-plantpart@shoot-asv.qza" \
            --i-taxonomy "${DATA_FOLDER}/taxonomy.qza" \
            --o-visualization "${DATA_FOLDER}/ancombc2-barplot-gnotobiotic-plantpart@shoot-asv.qzv"

    micromamba run q2-ampl-2025_4 \
        qiime composition ancombc2 \
            --i-table "${DATA_FOLDER}/genus-table-gnotobiotic.qza" \
            --m-metadata-file "${CONF_FOLDER}/${METADATA_FILE}" \
            --p-fixed-effects-formula 'PlantPart' \
            --p-reference-levels 'PlantPart::Shoot' \
            --o-ancombc2-output "${DATA_FOLDER}/ancombc2-gnotobiotic-plantpart@shoot-genus.qza"

    mkdir -p "${BIOM_FOLDER}/tmp"
    micromamba run q2-ampl-2025_4 \
        qiime tools export \
            --input-path "${DATA_FOLDER}/ancombc2-gnotobiotic-plantpart@shoot-genus.qza" \
            --output-path "${BIOM_FOLDER}/tmp"
    for m in ${METRICS[@]} ; do
        sed '1d' "${BIOM_FOLDER}/tmp/${m}.jsonl" | \
            jq '.taxon,."(Intercept)",."PlantPart::Root",."PlantPart::Seed"' | \
            sed -e 's/^"//' -e 's/"$//' -e 's/null/NaN/' | paste - - - - | \
            sed "1iTaxon\t(Intercept)|${m}\tPlantPart::Root|${m}\tPlantPart::Seed|${m}" > "${BIOM_FOLDER}/tmp/${m}.tsv"
    done
    csvjoin -t -I -y 0 -c 'Taxon' "${BIOM_FOLDER}/tmp/lfc.tsv" "${BIOM_FOLDER}/tmp/se.tsv" "${BIOM_FOLDER}/tmp/W.tsv" "${BIOM_FOLDER}/tmp/p.tsv" "${BIOM_FOLDER}/tmp/q.tsv" "${BIOM_FOLDER}/tmp/diff.tsv" "${BIOM_FOLDER}/tmp/passed_ss.tsv" | \
        csvformat -T > "${BIOM_FOLDER}/ancombc2-gnotobiotic-plantpart@shoot-genus.tsv"
    rm -r "${BIOM_FOLDER}/tmp"

    micromamba run q2-ampl-2025_4 \
        qiime composition tabulate \
            --i-data "${DATA_FOLDER}/ancombc2-gnotobiotic-plantpart@shoot-genus.qza" \
            --o-visualization "${DATA_FOLDER}/ancombc2-gnotobiotic-plantpart@shoot-genus.qzv"

    micromamba run q2-ampl-2025_4 \
        qiime composition ancombc2-visualizer \
            --i-data "${DATA_FOLDER}/ancombc2-gnotobiotic-plantpart@shoot-genus.qza" \
            --o-visualization "${DATA_FOLDER}/ancombc2-barplot-gnotobiotic-plantpart@shoot-genus.qzv"

    # Check effect of PlantPart using 'Root' as reference
    micromamba run q2-ampl-2025_4 \
        qiime composition ancombc2 \
            --i-table "${DATA_FOLDER}/asv-table-gnotobiotic.qza" \
            --m-metadata-file "${CONF_FOLDER}/${METADATA_FILE}" \
            --p-fixed-effects-formula 'PlantPart' \
            --p-reference-levels 'PlantPart::Root' \
            --o-ancombc2-output "${DATA_FOLDER}/ancombc2-gnotobiotic-plantpart@root-asv.qza"

    mkdir -p "${BIOM_FOLDER}/tmp"
    micromamba run q2-ampl-2025_4 \
        qiime tools export \
            --input-path "${DATA_FOLDER}/ancombc2-gnotobiotic-plantpart@root-asv.qza" \
            --output-path "${BIOM_FOLDER}/tmp"
    for m in ${METRICS[@]} ; do
        sed '1d' "${BIOM_FOLDER}/tmp/${m}.jsonl" | \
            jq '.taxon,."(Intercept)",."PlantPart::Seed",."PlantPart::Shoot"' | \
            sed -e 's/^"//' -e 's/"$//' -e 's/null/NaN/' | paste - - - - | \
            sed "1iTaxon\t(Intercept)|${m}\tPlantPart::Seed|${m}\tPlantPart::Shoot|${m}" > "${BIOM_FOLDER}/tmp/${m}.tsv"
    done
    csvjoin -t -I -y 0 -c 'Taxon' "${BIOM_FOLDER}/tmp/lfc.tsv" "${BIOM_FOLDER}/tmp/se.tsv" "${BIOM_FOLDER}/tmp/W.tsv" "${BIOM_FOLDER}/tmp/p.tsv" "${BIOM_FOLDER}/tmp/q.tsv" "${BIOM_FOLDER}/tmp/diff.tsv" "${BIOM_FOLDER}/tmp/passed_ss.tsv" | \
        csvformat -T > "${BIOM_FOLDER}/ancombc2-gnotobiotic-plantpart@root-asv.tsv"
    rm -r "${BIOM_FOLDER}/tmp"

    micromamba run q2-ampl-2025_4 \
        qiime composition tabulate \
            --i-data "${DATA_FOLDER}/ancombc2-gnotobiotic-plantpart@root-asv.qza" \
            --o-visualization "${DATA_FOLDER}/ancombc2-gnotobiotic-plantpart@root-asv.qzv"

    micromamba run q2-ampl-2025_4 \
        qiime composition ancombc2-visualizer \
            --i-data "${DATA_FOLDER}/ancombc2-gnotobiotic-plantpart@root-asv.qza" \
            --i-taxonomy "${DATA_FOLDER}/taxonomy.qza" \
            --o-visualization "${DATA_FOLDER}/ancombc2-barplot-gnotobiotic-plantpart@root-asv.qzv"

    micromamba run q2-ampl-2025_4 \
        qiime composition ancombc2 \
            --i-table "${DATA_FOLDER}/genus-table-gnotobiotic.qza" \
            --m-metadata-file "${CONF_FOLDER}/${METADATA_FILE}" \
            --p-fixed-effects-formula 'PlantPart' \
            --p-reference-levels 'PlantPart::Root' \
            --o-ancombc2-output "${DATA_FOLDER}/ancombc2-gnotobiotic-plantpart@root-genus.qza"

    mkdir -p "${BIOM_FOLDER}/tmp"
    micromamba run q2-ampl-2025_4 \
        qiime tools export \
            --input-path "${DATA_FOLDER}/ancombc2-gnotobiotic-plantpart@root-genus.qza" \
            --output-path "${BIOM_FOLDER}/tmp"
    for m in ${METRICS[@]} ; do
        sed '1d' "${BIOM_FOLDER}/tmp/${m}.jsonl" | \
            jq '.taxon,."(Intercept)",."PlantPart::Seed",."PlantPart::Shoot"' | \
            sed -e 's/^"//' -e 's/"$//' -e 's/null/NaN/' | paste - - - - | \
            sed "1iTaxon\t(Intercept)|${m}\tPlantPart::Seed|${m}\tPlantPart::Shoot|${m}" > "${BIOM_FOLDER}/tmp/${m}.tsv"
    done
    csvjoin -t -I -y 0 -c 'Taxon' "${BIOM_FOLDER}/tmp/lfc.tsv" "${BIOM_FOLDER}/tmp/se.tsv" "${BIOM_FOLDER}/tmp/W.tsv" "${BIOM_FOLDER}/tmp/p.tsv" "${BIOM_FOLDER}/tmp/q.tsv" "${BIOM_FOLDER}/tmp/diff.tsv" "${BIOM_FOLDER}/tmp/passed_ss.tsv" | \
        csvformat -T > "${BIOM_FOLDER}/ancombc2-gnotobiotic-plantpart@root-genus.tsv"
    rm -r "${BIOM_FOLDER}/tmp"

    micromamba run q2-ampl-2025_4 \
        qiime composition tabulate \
            --i-data "${DATA_FOLDER}/ancombc2-gnotobiotic-plantpart@root-genus.qza" \
            --o-visualization "${DATA_FOLDER}/ancombc2-gnotobiotic-plantpart@root-genus.qzv"

    micromamba run q2-ampl-2025_4 \
        qiime composition ancombc2-visualizer \
            --i-data "${DATA_FOLDER}/ancombc2-gnotobiotic-plantpart@root-genus.qza" \
            --o-visualization "${DATA_FOLDER}/ancombc2-barplot-gnotobiotic-plantpart@root-genus.qzv"

    echo -e '[X] Analyzing gnotobiotic differential abundance with ANCOM-BC2\n'
}

ancombc2_differential_abundance_analysis_field(){
    echo '[ ] Analyzing field differential abundance with ANCOM-BC2'

    METRICS=( 'lfc' 'se' 'W' 'p' 'q' 'diff' 'passed_ss' )

    # Check effect of PlantPart using 'Seed' as reference
    micromamba run q2-ampl-2025_4 \
        qiime composition ancombc2 \
            --i-table "${DATA_FOLDER}/asv-table-field.qza" \
            --m-metadata-file "${CONF_FOLDER}/${METADATA_FILE}" \
            --p-fixed-effects-formula 'PlantPart' \
            --p-reference-levels 'PlantPart::Seed' \
            --o-ancombc2-output "${DATA_FOLDER}/ancombc2-field-plantpart@seed-asv.qza"

    mkdir -p "${BIOM_FOLDER}/tmp"
    micromamba run q2-ampl-2025_4 \
        qiime tools export \
            --input-path "${DATA_FOLDER}/ancombc2-field-plantpart@seed-asv.qza" \
            --output-path "${BIOM_FOLDER}/tmp"
    for m in ${METRICS[@]} ; do
        sed '1d' "${BIOM_FOLDER}/tmp/${m}.jsonl" | \
            jq '.taxon,."(Intercept)",."PlantPart::Root",."PlantPart::Shoot",."PlantPart::Soil",."PlantPart::Spike"' | \
            sed -e 's/^"//' -e 's/"$//' -e 's/null/NaN/' | paste - - - - - - | \
            sed "1iTaxon\t(Intercept)|${m}\tPlantPart::Root|${m}\tPlantPart::Shoot|${m}\tPlantPart::Soil|${m}\tPlantPart::Spike|${m}" > "${BIOM_FOLDER}/tmp/${m}.tsv"
    done
    csvjoin -t -I -y 0 -c 'Taxon' "${BIOM_FOLDER}/tmp/lfc.tsv" "${BIOM_FOLDER}/tmp/se.tsv" "${BIOM_FOLDER}/tmp/W.tsv" "${BIOM_FOLDER}/tmp/p.tsv" "${BIOM_FOLDER}/tmp/q.tsv" "${BIOM_FOLDER}/tmp/diff.tsv" "${BIOM_FOLDER}/tmp/passed_ss.tsv" | \
        csvformat -T > "${BIOM_FOLDER}/ancombc2-field-plantpart@seed-asv.tsv"
    rm -r "${BIOM_FOLDER}/tmp"

    micromamba run q2-ampl-2025_4 \
        qiime composition tabulate \
            --i-data "${DATA_FOLDER}/ancombc2-field-plantpart@seed-asv.qza" \
            --o-visualization "${DATA_FOLDER}/ancombc2-field-plantpart@seed-asv.qzv"

    micromamba run q2-ampl-2025_4 \
        qiime composition ancombc2-visualizer \
            --i-data "${DATA_FOLDER}/ancombc2-field-plantpart@seed-asv.qza" \
            --i-taxonomy "${DATA_FOLDER}/taxonomy.qza" \
            --o-visualization "${DATA_FOLDER}/ancombc2-barplot-field-plantpart@seed-asv.qzv"

    micromamba run q2-ampl-2025_4 \
        qiime composition ancombc2 \
            --i-table "${DATA_FOLDER}/genus-table-field.qza" \
            --m-metadata-file "${CONF_FOLDER}/${METADATA_FILE}" \
            --p-fixed-effects-formula 'PlantPart' \
            --p-reference-levels 'PlantPart::Seed' \
            --o-ancombc2-output "${DATA_FOLDER}/ancombc2-field-plantpart@seed-genus.qza"

    mkdir -p "${BIOM_FOLDER}/tmp"
    micromamba run q2-ampl-2025_4 \
        qiime tools export \
            --input-path "${DATA_FOLDER}/ancombc2-field-plantpart@seed-genus.qza" \
            --output-path "${BIOM_FOLDER}/tmp"
    for m in ${METRICS[@]} ; do
        sed '1d' "${BIOM_FOLDER}/tmp/${m}.jsonl" | \
            jq '.taxon,."(Intercept)",."PlantPart::Root",."PlantPart::Shoot",."PlantPart::Soil",."PlantPart::Spike"' | \
            sed -e 's/^"//' -e 's/"$//' -e 's/null/NaN/' | paste - - - - - - | \
            sed "1iTaxon\t(Intercept)|${m}\tPlantPart::Root|${m}\tPlantPart::Shoot|${m}\tPlantPart::Soil|${m}\tPlantPart::Spike|${m}" > "${BIOM_FOLDER}/tmp/${m}.tsv"
    done
    csvjoin -t -I -y 0 -c 'Taxon' "${BIOM_FOLDER}/tmp/lfc.tsv" "${BIOM_FOLDER}/tmp/se.tsv" "${BIOM_FOLDER}/tmp/W.tsv" "${BIOM_FOLDER}/tmp/p.tsv" "${BIOM_FOLDER}/tmp/q.tsv" "${BIOM_FOLDER}/tmp/diff.tsv" "${BIOM_FOLDER}/tmp/passed_ss.tsv" | \
        csvformat -T > "${BIOM_FOLDER}/ancombc2-field-plantpart@seed-genus.tsv"
    rm -r "${BIOM_FOLDER}/tmp"

    micromamba run q2-ampl-2025_4 \
        qiime composition tabulate \
            --i-data "${DATA_FOLDER}/ancombc2-field-plantpart@seed-genus.qza" \
            --o-visualization "${DATA_FOLDER}/ancombc2-field-plantpart@seed-genus.qzv"

    micromamba run q2-ampl-2025_4 \
        qiime composition ancombc2-visualizer \
            --i-data "${DATA_FOLDER}/ancombc2-field-plantpart@seed-genus.qza" \
            --o-visualization "${DATA_FOLDER}/ancombc2-barplot-field-plantpart@seed-genus.qzv"

    # Check effect of PlantPart using 'Root' as reference
    micromamba run q2-ampl-2025_4 \
        qiime composition ancombc2 \
            --i-table "${DATA_FOLDER}/asv-table-field.qza" \
            --m-metadata-file "${CONF_FOLDER}/${METADATA_FILE}" \
            --p-fixed-effects-formula 'PlantPart' \
            --p-reference-levels 'PlantPart::Root' \
            --o-ancombc2-output "${DATA_FOLDER}/ancombc2-field-plantpart@root-asv.qza"

    mkdir -p "${BIOM_FOLDER}/tmp"
    micromamba run q2-ampl-2025_4 \
        qiime tools export \
            --input-path "${DATA_FOLDER}/ancombc2-field-plantpart@root-asv.qza" \
            --output-path "${BIOM_FOLDER}/tmp"
    for m in ${METRICS[@]} ; do
        sed '1d' "${BIOM_FOLDER}/tmp/${m}.jsonl" | \
            jq '.taxon,."(Intercept)",."PlantPart::Seed",."PlantPart::Shoot",."PlantPart::Soil",."PlantPart::Spike"' | \
            sed -e 's/^"//' -e 's/"$//' -e 's/null/NaN/' | paste - - - - - - | \
            sed "1iTaxon\t(Intercept)|${m}\tPlantPart::Seed|${m}\tPlantPart::Shoot|${m}\tPlantPart::Soil|${m}\tPlantPart::Spike|${m}" > "${BIOM_FOLDER}/tmp/${m}.tsv"
    done
    csvjoin -t -I -y 0 -c 'Taxon' "${BIOM_FOLDER}/tmp/lfc.tsv" "${BIOM_FOLDER}/tmp/se.tsv" "${BIOM_FOLDER}/tmp/W.tsv" "${BIOM_FOLDER}/tmp/p.tsv" "${BIOM_FOLDER}/tmp/q.tsv" "${BIOM_FOLDER}/tmp/diff.tsv" "${BIOM_FOLDER}/tmp/passed_ss.tsv" | \
        csvformat -T > "${BIOM_FOLDER}/ancombc2-field-plantpart@root-asv.tsv"
    rm -r "${BIOM_FOLDER}/tmp"

    micromamba run q2-ampl-2025_4 \
        qiime composition tabulate \
            --i-data "${DATA_FOLDER}/ancombc2-field-plantpart@root-asv.qza" \
            --o-visualization "${DATA_FOLDER}/ancombc2-field-plantpart@root-asv.qzv"

    micromamba run q2-ampl-2025_4 \
        qiime composition ancombc2-visualizer \
            --i-data "${DATA_FOLDER}/ancombc2-field-plantpart@root-asv.qza" \
            --i-taxonomy "${DATA_FOLDER}/taxonomy.qza" \
            --o-visualization "${DATA_FOLDER}/ancombc2-barplot-field-plantpart@root-asv.qzv"

    micromamba run q2-ampl-2025_4 \
        qiime composition ancombc2 \
            --i-table "${DATA_FOLDER}/genus-table-field.qza" \
            --m-metadata-file "${CONF_FOLDER}/${METADATA_FILE}" \
            --p-fixed-effects-formula 'PlantPart' \
            --p-reference-levels 'PlantPart::Root' \
            --o-ancombc2-output "${DATA_FOLDER}/ancombc2-field-plantpart@root-genus.qza"

    mkdir -p "${BIOM_FOLDER}/tmp"
    micromamba run q2-ampl-2025_4 \
        qiime tools export \
            --input-path "${DATA_FOLDER}/ancombc2-field-plantpart@root-genus.qza" \
            --output-path "${BIOM_FOLDER}/tmp"
    for m in ${METRICS[@]} ; do
        sed '1d' "${BIOM_FOLDER}/tmp/${m}.jsonl" | \
            jq '.taxon,."(Intercept)",."PlantPart::Seed",."PlantPart::Shoot",."PlantPart::Soil",."PlantPart::Spike"' | \
            sed -e 's/^"//' -e 's/"$//' -e 's/null/NaN/' | paste - - - - - - | \
            sed "1iTaxon\t(Intercept)|${m}\tPlantPart::Seed|${m}\tPlantPart::Shoot|${m}\tPlantPart::Soil|${m}\tPlantPart::Spike|${m}" > "${BIOM_FOLDER}/tmp/${m}.tsv"
    done
    csvjoin -t -I -y 0 -c 'Taxon' "${BIOM_FOLDER}/tmp/lfc.tsv" "${BIOM_FOLDER}/tmp/se.tsv" "${BIOM_FOLDER}/tmp/W.tsv" "${BIOM_FOLDER}/tmp/p.tsv" "${BIOM_FOLDER}/tmp/q.tsv" "${BIOM_FOLDER}/tmp/diff.tsv" "${BIOM_FOLDER}/tmp/passed_ss.tsv" | \
        csvformat -T > "${BIOM_FOLDER}/ancombc2-field-plantpart@root-genus.tsv"
    rm -r "${BIOM_FOLDER}/tmp"

    micromamba run q2-ampl-2025_4 \
        qiime composition tabulate \
            --i-data "${DATA_FOLDER}/ancombc2-field-plantpart@root-genus.qza" \
            --o-visualization "${DATA_FOLDER}/ancombc2-field-plantpart@root-genus.qzv"

    micromamba run q2-ampl-2025_4 \
        qiime composition ancombc2-visualizer \
            --i-data "${DATA_FOLDER}/ancombc2-field-plantpart@root-genus.qza" \
            --o-visualization "${DATA_FOLDER}/ancombc2-barplot-field-plantpart@root-genus.qzv"

    # Check effect of PlantPart using 'Shoot' as reference
    micromamba run q2-ampl-2025_4 \
        qiime composition ancombc2 \
            --i-table "${DATA_FOLDER}/asv-table-field.qza" \
            --m-metadata-file "${CONF_FOLDER}/${METADATA_FILE}" \
            --p-fixed-effects-formula 'PlantPart' \
            --p-reference-levels 'PlantPart::Shoot' \
            --o-ancombc2-output "${DATA_FOLDER}/ancombc2-field-plantpart@shoot-asv.qza"

    mkdir -p "${BIOM_FOLDER}/tmp"
    micromamba run q2-ampl-2025_4 \
        qiime tools export \
            --input-path "${DATA_FOLDER}/ancombc2-field-plantpart@shoot-asv.qza" \
            --output-path "${BIOM_FOLDER}/tmp"
    for m in ${METRICS[@]} ; do
        sed '1d' "${BIOM_FOLDER}/tmp/${m}.jsonl" | \
            jq '.taxon,."(Intercept)",."PlantPart::Seed",."PlantPart::Root",."PlantPart::Soil",."PlantPart::Spike"' | \
            sed -e 's/^"//' -e 's/"$//' -e 's/null/NaN/' | paste - - - - - - | \
            sed "1iTaxon\t(Intercept)|${m}\tPlantPart::Seed|${m}\tPlantPart::Root|${m}\tPlantPart::Soil|${m}\tPlantPart::Spike|${m}" > "${BIOM_FOLDER}/tmp/${m}.tsv"
    done
    csvjoin -t -I -y 0 -c 'Taxon' "${BIOM_FOLDER}/tmp/lfc.tsv" "${BIOM_FOLDER}/tmp/se.tsv" "${BIOM_FOLDER}/tmp/W.tsv" "${BIOM_FOLDER}/tmp/p.tsv" "${BIOM_FOLDER}/tmp/q.tsv" "${BIOM_FOLDER}/tmp/diff.tsv" "${BIOM_FOLDER}/tmp/passed_ss.tsv" | \
        csvformat -T > "${BIOM_FOLDER}/ancombc2-field-plantpart@shoot-asv.tsv"
    rm -r "${BIOM_FOLDER}/tmp"

    micromamba run q2-ampl-2025_4 \
        qiime composition tabulate \
            --i-data "${DATA_FOLDER}/ancombc2-field-plantpart@shoot-asv.qza" \
            --o-visualization "${DATA_FOLDER}/ancombc2-field-plantpart@shoot-asv.qzv"

    micromamba run q2-ampl-2025_4 \
        qiime composition ancombc2-visualizer \
            --i-data "${DATA_FOLDER}/ancombc2-field-plantpart@shoot-asv.qza" \
            --i-taxonomy "${DATA_FOLDER}/taxonomy.qza" \
            --o-visualization "${DATA_FOLDER}/ancombc2-barplot-field-plantpart@shoot-asv.qzv"

    micromamba run q2-ampl-2025_4 \
        qiime composition ancombc2 \
            --i-table "${DATA_FOLDER}/genus-table-field.qza" \
            --m-metadata-file "${CONF_FOLDER}/${METADATA_FILE}" \
            --p-fixed-effects-formula 'PlantPart' \
            --p-reference-levels 'PlantPart::Shoot' \
            --o-ancombc2-output "${DATA_FOLDER}/ancombc2-field-plantpart@shoot-genus.qza"

    mkdir -p "${BIOM_FOLDER}/tmp"
    micromamba run q2-ampl-2025_4 \
        qiime tools export \
            --input-path "${DATA_FOLDER}/ancombc2-field-plantpart@shoot-genus.qza" \
            --output-path "${BIOM_FOLDER}/tmp"
    for m in ${METRICS[@]} ; do
        sed '1d' "${BIOM_FOLDER}/tmp/${m}.jsonl" | \
            jq '.taxon,."(Intercept)",."PlantPart::Seed",."PlantPart::Root",."PlantPart::Soil",."PlantPart::Spike"' | \
            sed -e 's/^"//' -e 's/"$//' -e 's/null/NaN/' | paste - - - - - - | \
            sed "1iTaxon\t(Intercept)|${m}\tPlantPart::Seed|${m}\tPlantPart::Root|${m}\tPlantPart::Soil|${m}\tPlantPart::Spike|${m}" > "${BIOM_FOLDER}/tmp/${m}.tsv"
    done
    csvjoin -t -I -y 0 -c 'Taxon' "${BIOM_FOLDER}/tmp/lfc.tsv" "${BIOM_FOLDER}/tmp/se.tsv" "${BIOM_FOLDER}/tmp/W.tsv" "${BIOM_FOLDER}/tmp/p.tsv" "${BIOM_FOLDER}/tmp/q.tsv" "${BIOM_FOLDER}/tmp/diff.tsv" "${BIOM_FOLDER}/tmp/passed_ss.tsv" | \
        csvformat -T > "${BIOM_FOLDER}/ancombc2-field-plantpart@shoot-genus.tsv"
    rm -r "${BIOM_FOLDER}/tmp"

    micromamba run q2-ampl-2025_4 \
        qiime composition tabulate \
            --i-data "${DATA_FOLDER}/ancombc2-field-plantpart@shoot-genus.qza" \
            --o-visualization "${DATA_FOLDER}/ancombc2-field-plantpart@shoot-genus.qzv"

    micromamba run q2-ampl-2025_4 \
        qiime composition ancombc2-visualizer \
            --i-data "${DATA_FOLDER}/ancombc2-field-plantpart@shoot-genus.qza" \
            --o-visualization "${DATA_FOLDER}/ancombc2-barplot-field-plantpart@shoot-genus.qzv"

    # Check effect of PlantPart using 'Soil' as reference
    micromamba run q2-ampl-2025_4 \
        qiime composition ancombc2 \
            --i-table "${DATA_FOLDER}/asv-table-field.qza" \
            --m-metadata-file "${CONF_FOLDER}/${METADATA_FILE}" \
            --p-fixed-effects-formula 'PlantPart' \
            --p-reference-levels 'PlantPart::Soil' \
            --o-ancombc2-output "${DATA_FOLDER}/ancombc2-field-plantpart@soil-asv.qza"

    mkdir -p "${BIOM_FOLDER}/tmp"
    micromamba run q2-ampl-2025_4 \
        qiime tools export \
            --input-path "${DATA_FOLDER}/ancombc2-field-plantpart@soil-asv.qza" \
            --output-path "${BIOM_FOLDER}/tmp"
    for m in ${METRICS[@]} ; do
        sed '1d' "${BIOM_FOLDER}/tmp/${m}.jsonl" | \
            jq '.taxon,."(Intercept)",."PlantPart::Seed",."PlantPart::Root",."PlantPart::Shoot",."PlantPart::Spike"' | \
            sed -e 's/^"//' -e 's/"$//' -e 's/null/NaN/' | paste - - - - - - | \
            sed "1iTaxon\t(Intercept)|${m}\tPlantPart::Seed|${m}\tPlantPart::Root|${m}\tPlantPart::Shoot|${m}\tPlantPart::Spike|${m}" > "${BIOM_FOLDER}/tmp/${m}.tsv"
    done
    csvjoin -t -I -y 0 -c 'Taxon' "${BIOM_FOLDER}/tmp/lfc.tsv" "${BIOM_FOLDER}/tmp/se.tsv" "${BIOM_FOLDER}/tmp/W.tsv" "${BIOM_FOLDER}/tmp/p.tsv" "${BIOM_FOLDER}/tmp/q.tsv" "${BIOM_FOLDER}/tmp/diff.tsv" "${BIOM_FOLDER}/tmp/passed_ss.tsv" | \
        csvformat -T > "${BIOM_FOLDER}/ancombc2-field-plantpart@soil-asv.tsv"
    rm -r "${BIOM_FOLDER}/tmp"

    micromamba run q2-ampl-2025_4 \
        qiime composition tabulate \
            --i-data "${DATA_FOLDER}/ancombc2-field-plantpart@soil-asv.qza" \
            --o-visualization "${DATA_FOLDER}/ancombc2-field-plantpart@soil-asv.qzv"

    micromamba run q2-ampl-2025_4 \
        qiime composition ancombc2-visualizer \
            --i-data "${DATA_FOLDER}/ancombc2-field-plantpart@soil-asv.qza" \
            --i-taxonomy "${DATA_FOLDER}/taxonomy.qza" \
            --o-visualization "${DATA_FOLDER}/ancombc2-barplot-field-plantpart@soil-asv.qzv"

    micromamba run q2-ampl-2025_4 \
        qiime composition ancombc2 \
            --i-table "${DATA_FOLDER}/genus-table-field.qza" \
            --m-metadata-file "${CONF_FOLDER}/${METADATA_FILE}" \
            --p-fixed-effects-formula 'PlantPart' \
            --p-reference-levels 'PlantPart::Soil' \
            --o-ancombc2-output "${DATA_FOLDER}/ancombc2-field-plantpart@soil-genus.qza"

    mkdir -p "${BIOM_FOLDER}/tmp"
    micromamba run q2-ampl-2025_4 \
        qiime tools export \
            --input-path "${DATA_FOLDER}/ancombc2-field-plantpart@soil-genus.qza" \
            --output-path "${BIOM_FOLDER}/tmp"
    for m in ${METRICS[@]} ; do
        sed '1d' "${BIOM_FOLDER}/tmp/${m}.jsonl" | \
            jq '.taxon,."(Intercept)",."PlantPart::Seed",."PlantPart::Root",."PlantPart::Shoot",."PlantPart::Spike"' | \
            sed -e 's/^"//' -e 's/"$//' -e 's/null/NaN/' | paste - - - - - - | \
            sed "1iTaxon\t(Intercept)|${m}\tPlantPart::Seed|${m}\tPlantPart::Root|${m}\tPlantPart::Shoot|${m}\tPlantPart::Spike|${m}" > "${BIOM_FOLDER}/tmp/${m}.tsv"
    done
    csvjoin -t -I -y 0 -c 'Taxon' "${BIOM_FOLDER}/tmp/lfc.tsv" "${BIOM_FOLDER}/tmp/se.tsv" "${BIOM_FOLDER}/tmp/W.tsv" "${BIOM_FOLDER}/tmp/p.tsv" "${BIOM_FOLDER}/tmp/q.tsv" "${BIOM_FOLDER}/tmp/diff.tsv" "${BIOM_FOLDER}/tmp/passed_ss.tsv" | \
        csvformat -T > "${BIOM_FOLDER}/ancombc2-field-plantpart@soil-genus.tsv"
    rm -r "${BIOM_FOLDER}/tmp"

    micromamba run q2-ampl-2025_4 \
        qiime composition tabulate \
            --i-data "${DATA_FOLDER}/ancombc2-field-plantpart@soil-genus.qza" \
            --o-visualization "${DATA_FOLDER}/ancombc2-field-plantpart@soil-genus.qzv"

    micromamba run q2-ampl-2025_4 \
        qiime composition ancombc2-visualizer \
            --i-data "${DATA_FOLDER}/ancombc2-field-plantpart@soil-genus.qza" \
            --o-visualization "${DATA_FOLDER}/ancombc2-barplot-field-plantpart@soil-genus.qzv"

    # Check effect of PlantPart using 'Spike' as reference
    micromamba run q2-ampl-2025_4 \
        qiime composition ancombc2 \
            --i-table "${DATA_FOLDER}/asv-table-field.qza" \
            --m-metadata-file "${CONF_FOLDER}/${METADATA_FILE}" \
            --p-fixed-effects-formula 'PlantPart' \
            --p-reference-levels 'PlantPart::Spike' \
            --o-ancombc2-output "${DATA_FOLDER}/ancombc2-field-plantpart@spike-asv.qza"

    mkdir -p "${BIOM_FOLDER}/tmp"
    micromamba run q2-ampl-2025_4 \
        qiime tools export \
            --input-path "${DATA_FOLDER}/ancombc2-field-plantpart@spike-asv.qza" \
            --output-path "${BIOM_FOLDER}/tmp"
    for m in ${METRICS[@]} ; do
        sed '1d' "${BIOM_FOLDER}/tmp/${m}.jsonl" | \
            jq '.taxon,."(Intercept)",."PlantPart::Seed",."PlantPart::Root",."PlantPart::Shoot",."PlantPart::Soil"' | \
            sed -e 's/^"//' -e 's/"$//' -e 's/null/NaN/' | paste - - - - - - | \
            sed "1iTaxon\t(Intercept)|${m}\tPlantPart::Seed|${m}\tPlantPart::Root|${m}\tPlantPart::Shoot|${m}\tPlantPart::Soil|${m}" > "${BIOM_FOLDER}/tmp/${m}.tsv"
    done
    csvjoin -t -I -y 0 -c 'Taxon' "${BIOM_FOLDER}/tmp/lfc.tsv" "${BIOM_FOLDER}/tmp/se.tsv" "${BIOM_FOLDER}/tmp/W.tsv" "${BIOM_FOLDER}/tmp/p.tsv" "${BIOM_FOLDER}/tmp/q.tsv" "${BIOM_FOLDER}/tmp/diff.tsv" "${BIOM_FOLDER}/tmp/passed_ss.tsv" | \
        csvformat -T > "${BIOM_FOLDER}/ancombc2-field-plantpart@spike-asv.tsv"
    rm -r "${BIOM_FOLDER}/tmp"

    micromamba run q2-ampl-2025_4 \
        qiime composition tabulate \
            --i-data "${DATA_FOLDER}/ancombc2-field-plantpart@spike-asv.qza" \
            --o-visualization "${DATA_FOLDER}/ancombc2-field-plantpart@spike-asv.qzv"

    micromamba run q2-ampl-2025_4 \
        qiime composition ancombc2-visualizer \
            --i-data "${DATA_FOLDER}/ancombc2-field-plantpart@spike-asv.qza" \
            --i-taxonomy "${DATA_FOLDER}/taxonomy.qza" \
            --o-visualization "${DATA_FOLDER}/ancombc2-barplot-field-plantpart@spike-asv.qzv"

    micromamba run q2-ampl-2025_4 \
        qiime composition ancombc2 \
            --i-table "${DATA_FOLDER}/genus-table-field.qza" \
            --m-metadata-file "${CONF_FOLDER}/${METADATA_FILE}" \
            --p-fixed-effects-formula 'PlantPart' \
            --p-reference-levels 'PlantPart::Spike' \
            --o-ancombc2-output "${DATA_FOLDER}/ancombc2-field-plantpart@spike-genus.qza"

    mkdir -p "${BIOM_FOLDER}/tmp"
    micromamba run q2-ampl-2025_4 \
        qiime tools export \
            --input-path "${DATA_FOLDER}/ancombc2-field-plantpart@spike-genus.qza" \
            --output-path "${BIOM_FOLDER}/tmp"
    for m in ${METRICS[@]} ; do
        sed '1d' "${BIOM_FOLDER}/tmp/${m}.jsonl" | \
            jq '.taxon,."(Intercept)",."PlantPart::Seed",."PlantPart::Root",."PlantPart::Shoot",."PlantPart::Soil"' | \
            sed -e 's/^"//' -e 's/"$//' -e 's/null/NaN/' | paste - - - - - - | \
            sed "1iTaxon\t(Intercept)|${m}\tPlantPart::Seed|${m}\tPlantPart::Root|${m}\tPlantPart::Shoot|${m}\tPlantPart::Soil|${m}" > "${BIOM_FOLDER}/tmp/${m}.tsv"
    done
    csvjoin -t -I -y 0 -c 'Taxon' "${BIOM_FOLDER}/tmp/lfc.tsv" "${BIOM_FOLDER}/tmp/se.tsv" "${BIOM_FOLDER}/tmp/W.tsv" "${BIOM_FOLDER}/tmp/p.tsv" "${BIOM_FOLDER}/tmp/q.tsv" "${BIOM_FOLDER}/tmp/diff.tsv" "${BIOM_FOLDER}/tmp/passed_ss.tsv" | \
        csvformat -T > "${BIOM_FOLDER}/ancombc2-field-plantpart@spike-genus.tsv"
    rm -r "${BIOM_FOLDER}/tmp"

    micromamba run q2-ampl-2025_4 \
        qiime composition tabulate \
            --i-data "${DATA_FOLDER}/ancombc2-field-plantpart@spike-genus.qza" \
            --o-visualization "${DATA_FOLDER}/ancombc2-field-plantpart@spike-genus.qzv"

    micromamba run q2-ampl-2025_4 \
        qiime composition ancombc2-visualizer \
            --i-data "${DATA_FOLDER}/ancombc2-field-plantpart@spike-genus.qza" \
            --o-visualization "${DATA_FOLDER}/ancombc2-barplot-field-plantpart@spike-genus.qzv"

    # Check effect of WheatSpecies using 'T.aestivum' as reference
    micromamba run q2-ampl-2025_4 \
        qiime composition ancombc2 \
            --i-table "${DATA_FOLDER}/asv-table-field.qza" \
            --m-metadata-file "${CONF_FOLDER}/${METADATA_FILE}" \
            --p-fixed-effects-formula 'WheatSpecies' \
            --p-reference-levels 'WheatSpecies::T.aestivum' \
            --o-ancombc2-output "${DATA_FOLDER}/ancombc2-field-wheatspecies@aestivum-asv.qza"

    mkdir -p "${BIOM_FOLDER}/tmp"
    micromamba run q2-ampl-2025_4 \
        qiime tools export \
            --input-path "${DATA_FOLDER}/ancombc2-field-wheatspecies@aestivum-asv.qza" \
            --output-path "${BIOM_FOLDER}/tmp"
    for m in ${METRICS[@]} ; do
        sed '1d' "${BIOM_FOLDER}/tmp/${m}.jsonl" | \
            jq '.taxon,."(Intercept)",."WheatSpecies::T.durum",."WheatSpecies::T.spelta",."WheatSpecies::T.turgidum"' | \
            sed -e 's/^"//' -e 's/"$//' -e 's/null/NaN/' | paste - - - - - | \
            sed "1iTaxon\t(Intercept)|${m}\tWheatSpecies::T.durum|${m}\tWheatSpecies::T.spelta|${m}\tWheatSpecies::T.turgidum|${m}" > "${BIOM_FOLDER}/tmp/${m}.tsv"
    done
    csvjoin -t -I -y 0 -c 'Taxon' "${BIOM_FOLDER}/tmp/lfc.tsv" "${BIOM_FOLDER}/tmp/se.tsv" "${BIOM_FOLDER}/tmp/W.tsv" "${BIOM_FOLDER}/tmp/p.tsv" "${BIOM_FOLDER}/tmp/q.tsv" "${BIOM_FOLDER}/tmp/diff.tsv" "${BIOM_FOLDER}/tmp/passed_ss.tsv" | \
        csvformat -T > "${BIOM_FOLDER}/ancombc2-field-wheatspecies@aestivum-asv.tsv"
    rm -r "${BIOM_FOLDER}/tmp"

    micromamba run q2-ampl-2025_4 \
        qiime composition tabulate \
            --i-data "${DATA_FOLDER}/ancombc2-field-wheatspecies@aestivum-asv.qza" \
            --o-visualization "${DATA_FOLDER}/ancombc2-field-wheatspecies@aestivum-asv.qzv"

    micromamba run q2-ampl-2025_4 \
        qiime composition ancombc2-visualizer \
            --i-data "${DATA_FOLDER}/ancombc2-field-wheatspecies@aestivum-asv.qza" \
            --i-taxonomy "${DATA_FOLDER}/taxonomy.qza" \
            --o-visualization "${DATA_FOLDER}/ancombc2-barplot-field-wheatspecies@aestivum-asv.qzv"

    micromamba run q2-ampl-2025_4 \
        qiime composition ancombc2 \
            --i-table "${DATA_FOLDER}/genus-table-field.qza" \
            --m-metadata-file "${CONF_FOLDER}/${METADATA_FILE}" \
            --p-fixed-effects-formula 'WheatSpecies' \
            --p-reference-levels 'WheatSpecies::T.aestivum' \
            --o-ancombc2-output "${DATA_FOLDER}/ancombc2-field-wheatspecies@aestivum-genus.qza"

    mkdir -p "${BIOM_FOLDER}/tmp"
    micromamba run q2-ampl-2025_4 \
        qiime tools export \
            --input-path "${DATA_FOLDER}/ancombc2-field-wheatspecies@aestivum-genus.qza" \
            --output-path "${BIOM_FOLDER}/tmp"
    for m in ${METRICS[@]} ; do
        sed '1d' "${BIOM_FOLDER}/tmp/${m}.jsonl" | \
            jq '.taxon,."(Intercept)",."WheatSpecies::T.durum",."WheatSpecies::T.spelta",."WheatSpecies::T.turgidum"' | \
            sed -e 's/^"//' -e 's/"$//' -e 's/null/NaN/' | paste - - - - - | \
            sed "1iTaxon\t(Intercept)|${m}\tWheatSpecies::T.durum|${m}\tWheatSpecies::T.spelta|${m}\tWheatSpecies::T.turgidum|${m}" > "${BIOM_FOLDER}/tmp/${m}.tsv"
    done
    csvjoin -t -I -y 0 -c 'Taxon' "${BIOM_FOLDER}/tmp/lfc.tsv" "${BIOM_FOLDER}/tmp/se.tsv" "${BIOM_FOLDER}/tmp/W.tsv" "${BIOM_FOLDER}/tmp/p.tsv" "${BIOM_FOLDER}/tmp/q.tsv" "${BIOM_FOLDER}/tmp/diff.tsv" "${BIOM_FOLDER}/tmp/passed_ss.tsv" | \
        csvformat -T > "${BIOM_FOLDER}/ancombc2-field-wheatspecies@aestivum-genus.tsv"
    rm -r "${BIOM_FOLDER}/tmp"

    micromamba run q2-ampl-2025_4 \
        qiime composition tabulate \
            --i-data "${DATA_FOLDER}/ancombc2-field-wheatspecies@aestivum-genus.qza" \
            --o-visualization "${DATA_FOLDER}/ancombc2-field-wheatspecies@aestivum-genus.qzv"

    micromamba run q2-ampl-2025_4 \
        qiime composition ancombc2-visualizer \
            --i-data "${DATA_FOLDER}/ancombc2-field-wheatspecies@aestivum-genus.qza" \
            --o-visualization "${DATA_FOLDER}/ancombc2-barplot-field-wheatspecies@aestivum-genus.qzv"

    # Check effect of (mixed column) PlantPart_WheatSpecies using 'Seed_T.aestivum' as reference
    micromamba run q2-ampl-2025_4 \
        qiime composition ancombc2 \
            --i-table "${DATA_FOLDER}/asv-table-field.qza" \
            --m-metadata-file "${CONF_FOLDER}/${METADATA_FILE}" \
            --p-fixed-effects-formula 'PlantPart_WheatSpecies' \
            --p-reference-levels 'PlantPart_WheatSpecies::Seed_T.aestivum' \
            --o-ancombc2-output "${DATA_FOLDER}/ancombc2-field-plantpart_wheatspecies@seed_aestivum-asv.qza"

    mkdir -p "${BIOM_FOLDER}/tmp"
    micromamba run q2-ampl-2025_4 \
        qiime tools export \
            --input-path "${DATA_FOLDER}/ancombc2-field-plantpart_wheatspecies@seed_aestivum-asv.qza" \
            --output-path "${BIOM_FOLDER}/tmp"
    for m in ${METRICS[@]} ; do
        sed '1d' "${BIOM_FOLDER}/tmp/${m}.jsonl" | \
            jq '.taxon,."(Intercept)",."PlantPart_WheatSpecies::Root_T.aestivum",."PlantPart_WheatSpecies::Root_T.durum",."PlantPart_WheatSpecies::Root_T.spelta",."PlantPart_WheatSpecies::Root_T.turgidum",."PlantPart_WheatSpecies::Seed_T.durum",."PlantPart_WheatSpecies::Seed_T.spelta",."PlantPart_WheatSpecies::Seed_T.turgidum",."PlantPart_WheatSpecies::Shoot_T.aestivum",."PlantPart_WheatSpecies::Shoot_T.durum",."PlantPart_WheatSpecies::Shoot_T.spelta",."PlantPart_WheatSpecies::Shoot_T.turgidum",."PlantPart_WheatSpecies::Soil_T.aestivum",."PlantPart_WheatSpecies::Soil_T.durum",."PlantPart_WheatSpecies::Soil_T.spelta",."PlantPart_WheatSpecies::Soil_T.turgidum",."PlantPart_WheatSpecies::Spike_T.aestivum",."PlantPart_WheatSpecies::Spike_T.durum",."PlantPart_WheatSpecies::Spike_T.spelta",."PlantPart_WheatSpecies::Spike_T.turgidum"' | \
            sed -e 's/^"//' -e 's/"$//' -e 's/null/NaN/' | paste - - - - - - - - - - - - - - - - - - - - - | \
            sed "1iTaxon\t(Intercept)|${m}\tPlantPart_WheatSpecies::Root_T.aestivum|${m}\tPlantPart_WheatSpecies::Root_T.durum|${m}\tPlantPart_WheatSpecies::Root_T.spelta|${m}\tPlantPart_WheatSpecies::Root_T.turgidum|${m}\tPlantPart_WheatSpecies::Seed_T.durum|${m}\tPlantPart_WheatSpecies::Seed_T.spelta|${m}\tPlantPart_WheatSpecies::Seed_T.turgidum|${m}\tPlantPart_WheatSpecies::Shoot_T.aestivum|${m}\tPlantPart_WheatSpecies::Shoot_T.durum|${m}\tPlantPart_WheatSpecies::Shoot_T.spelta|${m}\tPlantPart_WheatSpecies::Shoot_T.turgidum|${m}\tPlantPart_WheatSpecies::Soil_T.aestivum|${m}\tPlantPart_WheatSpecies::Soil_T.durum|${m}\tPlantPart_WheatSpecies::Soil_T.spelta|${m}\tPlantPart_WheatSpecies::Soil_T.turgidum|${m}\tPlantPart_WheatSpecies::Spike_T.aestivum|${m}\tPlantPart_WheatSpecies::Spike_T.durum|${m}\tPlantPart_WheatSpecies::Spike_T.spelta|${m}\tPlantPart_WheatSpecies::Spike_T.turgidum|${m}" > "${BIOM_FOLDER}/tmp/${m}.tsv"
    done
    csvjoin -t -I -y 0 -c 'Taxon' "${BIOM_FOLDER}/tmp/lfc.tsv" "${BIOM_FOLDER}/tmp/se.tsv" "${BIOM_FOLDER}/tmp/W.tsv" "${BIOM_FOLDER}/tmp/p.tsv" "${BIOM_FOLDER}/tmp/q.tsv" "${BIOM_FOLDER}/tmp/diff.tsv" "${BIOM_FOLDER}/tmp/passed_ss.tsv" | \
        csvformat -T > "${BIOM_FOLDER}/ancombc2-field-plantpart_wheatspecies@seed_aestivum-asv.tsv"
    rm -r "${BIOM_FOLDER}/tmp"

    micromamba run q2-ampl-2025_4 \
        qiime composition tabulate \
            --i-data "${DATA_FOLDER}/ancombc2-field-plantpart_wheatspecies@seed_aestivum-asv.qza" \
            --o-visualization "${DATA_FOLDER}/ancombc2-field-plantpart_wheatspecies@seed_aestivum-asv.qzv"

    micromamba run q2-ampl-2025_4 \
        qiime composition ancombc2-visualizer \
            --i-data "${DATA_FOLDER}/ancombc2-field-plantpart_wheatspecies@seed_aestivum-asv.qza" \
            --i-taxonomy "${DATA_FOLDER}/taxonomy.qza" \
            --o-visualization "${DATA_FOLDER}/ancombc2-barplot-field-plantpart_wheatspecies@seed_aestivum-asv.qzv"

    micromamba run q2-ampl-2025_4 \
        qiime composition ancombc2 \
            --i-table "${DATA_FOLDER}/genus-table-field.qza" \
            --m-metadata-file "${CONF_FOLDER}/${METADATA_FILE}" \
            --p-fixed-effects-formula 'PlantPart_WheatSpecies' \
            --p-reference-levels 'PlantPart_WheatSpecies::Seed_T.aestivum' \
            --o-ancombc2-output "${DATA_FOLDER}/ancombc2-field-plantpart_wheatspecies@seed_aestivum-genus.qza"

    mkdir -p "${BIOM_FOLDER}/tmp"
    micromamba run q2-ampl-2025_4 \
        qiime tools export \
            --input-path "${DATA_FOLDER}/ancombc2-field-plantpart_wheatspecies@seed_aestivum-genus.qza" \
            --output-path "${BIOM_FOLDER}/tmp"
    for m in ${METRICS[@]} ; do
        sed '1d' "${BIOM_FOLDER}/tmp/${m}.jsonl" | \
            jq '.taxon,."(Intercept)",."PlantPart_WheatSpecies::Root_T.aestivum",."PlantPart_WheatSpecies::Root_T.durum",."PlantPart_WheatSpecies::Root_T.spelta",."PlantPart_WheatSpecies::Root_T.turgidum",."PlantPart_WheatSpecies::Seed_T.durum",."PlantPart_WheatSpecies::Seed_T.spelta",."PlantPart_WheatSpecies::Seed_T.turgidum",."PlantPart_WheatSpecies::Shoot_T.aestivum",."PlantPart_WheatSpecies::Shoot_T.durum",."PlantPart_WheatSpecies::Shoot_T.spelta",."PlantPart_WheatSpecies::Shoot_T.turgidum",."PlantPart_WheatSpecies::Soil_T.aestivum",."PlantPart_WheatSpecies::Soil_T.durum",."PlantPart_WheatSpecies::Soil_T.spelta",."PlantPart_WheatSpecies::Soil_T.turgidum",."PlantPart_WheatSpecies::Spike_T.aestivum",."PlantPart_WheatSpecies::Spike_T.durum",."PlantPart_WheatSpecies::Spike_T.spelta",."PlantPart_WheatSpecies::Spike_T.turgidum"' | \
            sed -e 's/^"//' -e 's/"$//' -e 's/null/NaN/' | paste - - - - - - - - - - - - - - - - - - - - - | \
            sed "1iTaxon\t(Intercept)|${m}\tPlantPart_WheatSpecies::Root_T.aestivum|${m}\tPlantPart_WheatSpecies::Root_T.durum|${m}\tPlantPart_WheatSpecies::Root_T.spelta|${m}\tPlantPart_WheatSpecies::Root_T.turgidum|${m}\tPlantPart_WheatSpecies::Seed_T.durum|${m}\tPlantPart_WheatSpecies::Seed_T.spelta|${m}\tPlantPart_WheatSpecies::Seed_T.turgidum|${m}\tPlantPart_WheatSpecies::Shoot_T.aestivum|${m}\tPlantPart_WheatSpecies::Shoot_T.durum|${m}\tPlantPart_WheatSpecies::Shoot_T.spelta|${m}\tPlantPart_WheatSpecies::Shoot_T.turgidum|${m}\tPlantPart_WheatSpecies::Soil_T.aestivum|${m}\tPlantPart_WheatSpecies::Soil_T.durum|${m}\tPlantPart_WheatSpecies::Soil_T.spelta|${m}\tPlantPart_WheatSpecies::Soil_T.turgidum|${m}\tPlantPart_WheatSpecies::Spike_T.aestivum|${m}\tPlantPart_WheatSpecies::Spike_T.durum|${m}\tPlantPart_WheatSpecies::Spike_T.spelta|${m}\tPlantPart_WheatSpecies::Spike_T.turgidum|${m}" > "${BIOM_FOLDER}/tmp/${m}.tsv"
    done
    csvjoin -t -I -y 0 -c 'Taxon' "${BIOM_FOLDER}/tmp/lfc.tsv" "${BIOM_FOLDER}/tmp/se.tsv" "${BIOM_FOLDER}/tmp/W.tsv" "${BIOM_FOLDER}/tmp/p.tsv" "${BIOM_FOLDER}/tmp/q.tsv" "${BIOM_FOLDER}/tmp/diff.tsv" "${BIOM_FOLDER}/tmp/passed_ss.tsv" | \
        csvformat -T > "${BIOM_FOLDER}/ancombc2-field-plantpart_wheatspecies@seed_aestivum-genus.tsv"
    rm -r "${BIOM_FOLDER}/tmp"

    micromamba run q2-ampl-2025_4 \
        qiime composition tabulate \
            --i-data "${DATA_FOLDER}/ancombc2-field-plantpart_wheatspecies@seed_aestivum-genus.qza" \
            --o-visualization "${DATA_FOLDER}/ancombc2-field-plantpart_wheatspecies@seed_aestivum-genus.qzv"

    micromamba run q2-ampl-2025_4 \
        qiime composition ancombc2-visualizer \
            --i-data "${DATA_FOLDER}/ancombc2-field-plantpart_wheatspecies@seed_aestivum-genus.qza" \
            --o-visualization "${DATA_FOLDER}/ancombc2-barplot-field-plantpart_wheatspecies@seed_aestivum-genus.qzv"

    # Check effect of (formula) PlantPart + WheatSpecies using 'Seed' and 'T.aestivum' as references
    micromamba run q2-ampl-2025_4 \
        qiime composition ancombc2 \
            --i-table "${DATA_FOLDER}/asv-table-field.qza" \
            --m-metadata-file "${CONF_FOLDER}/${METADATA_FILE}" \
            --p-fixed-effects-formula 'PlantPart + WheatSpecies' \
            --p-reference-levels 'PlantPart::Seed' 'WheatSpecies::T.aestivum' \
            --o-ancombc2-output "${DATA_FOLDER}/ancombc2-field-plantpart~wheatspecies@seed~aestivum-asv.qza"

    mkdir -p "${BIOM_FOLDER}/tmp"
    micromamba run q2-ampl-2025_4 \
        qiime tools export \
            --input-path "${DATA_FOLDER}/ancombc2-field-plantpart~wheatspecies@seed~aestivum-asv.qza" \
            --output-path "${BIOM_FOLDER}/tmp"
    for m in ${METRICS[@]} ; do
        sed '1d' "${BIOM_FOLDER}/tmp/${m}.jsonl" | \
            jq '.taxon,."(Intercept)",."PlantPart::Root",."PlantPart::Shoot",."PlantPart::Soil",."PlantPart::Spike",."WheatSpecies::T.durum",."WheatSpecies::T.spelta",."WheatSpecies::T.turgidum"' | \
            sed -e 's/^"//' -e 's/"$//' -e 's/null/NaN/' | paste - - - - - - - - - | \
            sed "1iTaxon\t(Intercept)|${m}\tPlantPart::Root|${m}\tPlantPart::Shoot|${m}\tPlantPart::Soil|${m}\tPlantPart::Spike|${m}\tWheatSpecies::T.durum|${m}\tWheatSpecies::T.spelta|${m}\tWheatSpecies::T.turgidum|${m}" > "${BIOM_FOLDER}/tmp/${m}.tsv"
    done
    csvjoin -t -I -y 0 -c 'Taxon' "${BIOM_FOLDER}/tmp/lfc.tsv" "${BIOM_FOLDER}/tmp/se.tsv" "${BIOM_FOLDER}/tmp/W.tsv" "${BIOM_FOLDER}/tmp/p.tsv" "${BIOM_FOLDER}/tmp/q.tsv" "${BIOM_FOLDER}/tmp/diff.tsv" "${BIOM_FOLDER}/tmp/passed_ss.tsv" | \
        csvformat -T > "${BIOM_FOLDER}/ancombc2-field-plantpart~wheatspecies@seed~aestivum-asv.tsv"
    rm -r "${BIOM_FOLDER}/tmp"

    micromamba run q2-ampl-2025_4 \
        qiime composition tabulate \
            --i-data "${DATA_FOLDER}/ancombc2-field-plantpart~wheatspecies@seed~aestivum-asv.qza" \
            --o-visualization "${DATA_FOLDER}/ancombc2-field-plantpart~wheatspecies@seed~aestivum-asv.qzv"

    micromamba run q2-ampl-2025_4" \
        qiime composition ancombc2-visualizer \
            --i-data "${DATA_FOLDER}/ancombc2-field-plantpart~wheatspecies@seed~aestivum-asv.qza" \
            --i-taxonomy "${DATA_FOLDER}/taxonomy.qza" \
            --o-visualization "${DATA_FOLDER}/ancombc2-barplot-field-plantpart~wheatspecies@seed~aestivum-asv.qzv"

    micromamba run q2-ampl-2025_4 \
        qiime composition ancombc2 \
            --i-table "${DATA_FOLDER}/genus-table-field.qza" \
            --m-metadata-file "${CONF_FOLDER}/${METADATA_FILE}" \
            --p-fixed-effects-formula 'PlantPart + WheatSpecies' \
            --p-reference-levels 'PlantPart::Seed' 'WheatSpecies::T.aestivum' \
            --o-ancombc2-output "${DATA_FOLDER}/ancombc2-field-plantpart~wheatspecies@seed~aestivum-genus.qza"

    mkdir -p "${BIOM_FOLDER}/tmp"
    micromamba run q2-ampl-2025_4 \
        qiime tools export \
            --input-path "${DATA_FOLDER}/ancombc2-field-plantpart~wheatspecies@seed~aestivum-genus.qza" \
            --output-path "${BIOM_FOLDER}/tmp"
    for m in ${METRICS[@]} ; do
        sed '1d' "${BIOM_FOLDER}/tmp/${m}.jsonl" | \
            jq '.taxon,."(Intercept)",."PlantPart::Root",."PlantPart::Shoot",."PlantPart::Soil",."PlantPart::Spike",."WheatSpecies::T.durum",."WheatSpecies::T.spelta",."WheatSpecies::T.turgidum"' | \
            sed -e 's/^"//' -e 's/"$//' -e 's/null/NaN/' | paste - - - - - - - - - | \
            sed "1iTaxon\t(Intercept)|${m}\tPlantPart::Root|${m}\tPlantPart::Shoot|${m}\tPlantPart::Soil|${m}\tPlantPart::Spike|${m}\tWheatSpecies::T.durum|${m}\tWheatSpecies::T.spelta|${m}\tWheatSpecies::T.turgidum|${m}" > "${BIOM_FOLDER}/tmp/${m}.tsv"
    done
    csvjoin -t -I -y 0 -c 'Taxon' "${BIOM_FOLDER}/tmp/lfc.tsv" "${BIOM_FOLDER}/tmp/se.tsv" "${BIOM_FOLDER}/tmp/W.tsv" "${BIOM_FOLDER}/tmp/p.tsv" "${BIOM_FOLDER}/tmp/q.tsv" "${BIOM_FOLDER}/tmp/diff.tsv" "${BIOM_FOLDER}/tmp/passed_ss.tsv" | \
        csvformat -T > "${BIOM_FOLDER}/ancombc2-field-plantpart~wheatspecies@seed~aestivum-genus.tsv"
    rm -r "${BIOM_FOLDER}/tmp"

    micromamba run q2-ampl-2025_4 \
        qiime composition tabulate \
            --i-data "${DATA_FOLDER}/ancombc2-field-plantpart~wheatspecies@seed~aestivum-genus.qza" \
            --o-visualization "${DATA_FOLDER}/ancombc2-field-plantpart~wheatspecies@seed~aestivum-genus.qzv"

    micromamba run q2-ampl-2025_4 \
        qiime composition ancombc2-visualizer \
            --i-data "${DATA_FOLDER}/ancombc2-field-plantpart~wheatspecies@seed~aestivum-genus.qza" \
            --o-visualization "${DATA_FOLDER}/ancombc2-barplot-field-plantpart~wheatspecies@seed~aestivum-genus.qzv"

    echo -e '[X] Analyzing field differential abundance with ANCOM-BC2\n'
}

main
