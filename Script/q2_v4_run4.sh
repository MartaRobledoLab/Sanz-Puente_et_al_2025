#!/usr/bin/env bash

main(){
    TIME_START=`date +%s`
    echo -e "Starting ...\n  `date -d@${TIME_START}`\n"

    set_variables
    create_folders

    create_manifest_file
    create_metadata_file

    check_raw_reads_quality
    check_primers_position

    import_raw_reads
    trim_primers_from_reads
    check_trimmed_reads_quality
    dada2_denoise_reads

    TIME_END=`date +%s`
    echo -e "\nJob successfully processed.\n  `date -d@${TIME_END}`"
    echo "  Execution time was `expr ${TIME_END} - ${TIME_START}` s."
}

set_variables(){
    echo '[ ] Setting variables'

    THREADS=40

    BASE_FOLDER='PRJNA1282304'

    RUNDATE='2022-05-09'
    RAW_READS_FOLDER="${BASE_FOLDER}/run_220509"

    CONF_FOLDER="${BASE_FOLDER}/etc_v4_run4"
    MANIFEST_FILE='sample-manifest.csv'
    METADATA_FILE='sample-metadata.tsv'

    REGION='V4'
    PRIMER_NAME='U515F-Bakt_805R' # Ref.: ProbeBase (https://probebase.net)
    FORWARD_PRIMER='GTGCCAGCMGCCGCGGTAA' # Full name (Alm et al.): S-*-Univ-0515-a-S-19
    REVERSE_PRIMER='GACTACHVGGGTATCTAATCC' # Full name (Alm et al.): S-D-Bact-0785-a-A-21
    FORWARD_PRIMER_LEN=`echo ${FORWARD_PRIMER} | tr -d '\n' | wc -c`
    REVERSE_PRIMER_LEN=`echo ${REVERSE_PRIMER} | tr -d '\n' | wc -c`

    CUTADAPT_ERROR_RATE=0.1
    CUTADAPT_MIN_LENGTH=50

    # Expected length = (805 - 21) - (515 + 19) + 1 = 251 (E. coli)
    # Max theoretical = 300 - 19 + 300 - 21 - 12 (overlap) = 548
    # Conf length = 219 + 152 - 12 (overlap) = 359 >= 251 (min required)
    DADA_TRUNC_LEN_F=219 # Q1 >= 34
    DADA_TRUNC_LEN_R=152 # Q1 >= 34 (except at 3 spots)
    DADA_TRIM_LEFT_F=0
    DADA_TRIM_LEFT_R=0
    DADA_MAX_EE_F=2.0
    DADA_MAX_EE_R=2.0
    DADA_TRUNC_Q=2

    DATA_FOLDER="${BASE_FOLDER}/data_v4_run4"
    RAW_FASTQC_FOLDER="${DATA_FOLDER}/fastqc_raw"
    TRIMMED_FASTQC_FOLDER="${DATA_FOLDER}/fastqc_trimmed"

    echo -e '[X] Setting variables\n'
}

create_folders(){
    echo '[ ] Creating folders'

    mkdir -p "${CONF_FOLDER}"

    mkdir -p "${DATA_FOLDER}"
    mkdir -p "${RAW_FASTQC_FOLDER}"
    mkdir -p "${TRIMMED_FASTQC_FOLDER}"

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
EOF

    echo -e '[X] Creating sample-manifest.csv file\n'
}

create_metadata_file(){
    echo '[ ] Creating sample-metadata.tsv file'

    cat > "${CONF_FOLDER}/${METADATA_FILE}" <<EOF
sample-id	SampleName	ForwardPrimer	ReversePrimer	PrimerName	Region	RunDate	PNAs	Beads	Enzyme	SampleType	ControlType	PlantPart	WheatSpecies	PlantPart_WheatSpecies	PloidyLevel	Location	Domestication	Control	Seed	Gnotobiotic	Field
#q2:types	categorical	categorical	categorical	categorical	categorical	categorical	categorical	categorical	categorical	categorical	categorical	categorical	categorical	categorical	categorical	categorical	categorical	categorical	categorical	categorical	categorical
AE_246	Ata_RU	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Seed	NA	Seed	A.tauschii	Seed_A.tauschii	2n	Asia	Ancestral	No	Yes	No	No
AE_498	Ata_UZ	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Seed	NA	Seed	A.tauschii	Seed_A.tauschii	2n	Asia	Ancestral	No	Yes	No	No
BRAJua	Tae_SP5	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Seed	NA	Seed	T.aestivum	Seed_T.aestivum	6n	Europe	Commercial	No	Yes	No	No
BRNVil	Tae_SP4	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Seed	NA	Seed	T.aestivum	Seed_T.aestivum	6n	Europe	Commercial	No	Yes	No	No
BSC2021_27	BSC202127	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Seed	NA	Seed	T.aestivum	Seed_T.aestivum	6n	Europe	Commercial	No	No	Yes	Yes
BSCMatR3_SH	Tae_Sh_C	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Field	NA	Shoot	T.aestivum	Shoot_T.aestivum	6n	America	Commercial	No	No	No	Yes
BSCMatR3_SP	Tae_Sp_C	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Field	NA	Spike	T.aestivum	Spike_T.aestivum	6n	America	Commercial	No	No	No	Yes
CH4_RO	Tae_Ro_A	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Field	NA	Root	T.aestivum	Root_T.aestivum	6n	America	Commercial	No	No	No	Yes
CH4_SH	Tae_Sh_A	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Field	NA	Shoot	T.aestivum	Shoot_T.aestivum	6n	America	Commercial	No	No	No	Yes
CH4_SP	Tae_Sp_A	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Field	NA	Spike	T.aestivum	Spike_T.aestivum	6n	America	Commercial	No	No	No	Yes
CH_RO	Tae_Ro_B	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Field	NA	Root	T.aestivum	Root_T.aestivum	6n	America	Commercial	No	No	No	Yes
CH_SO	Tae_So_B	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Field	NA	Soil	T.aestivum	Soil_T.aestivum	6n	America	Commercial	No	No	No	Yes
CPCRW_77	CPCRW_77	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Control	H2O	NA	NA	NA	NA	NA	NA	Yes	No	No	No
C_ex_30	C_ex_30	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Control	Extraction	NA	NA	NA	NA	NA	NA	Yes	No	No	No
C_ex_81	C_ex_81	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	Yes	Yes	Kappa	Control	Extraction	NA	NA	NA	NA	NA	NA	Yes	No	No	No
SE16	Tpe_SP5	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Seed	NA	Seed	T.spelta	Seed_T.spelta	6n	Europe	Ancestral	No	Yes	No	No
SE1620	Tpe_SP2	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Seed	NA	Seed	T.spelta	Seed_T.spelta	6n	Europe	Ancestral	No	Yes	No	Yes
SE16_RO	Tpe_RO	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Field	NA	Root	T.spelta	Root_T.spelta	6n	Europe	Ancestral	No	No	No	Yes
SE16_SH	Tpe_SH	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Field	NA	Shoot	T.spelta	Shoot_T.spelta	6n	Europe	Ancestral	No	No	No	Yes
SE16_SO	Tpe_SO	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Field	NA	Soil	T.spelta	Soil_T.spelta	6n	Europe	Ancestral	No	No	No	Yes
SE16_SP	Tpe_SP	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Field	NA	Spike	T.spelta	Spike_T.spelta	6n	Europe	Ancestral	No	No	No	Yes
SE36	Tpe_SP3	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Seed	NA	Seed	T.spelta	Seed_T.spelta	6n	Europe	Ancestral	No	Yes	No	No
SE64_74	Tpe_SP4	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Seed	NA	Seed	T.spelta	Seed_T.spelta	6n	Europe	Ancestral	No	Yes	No	Yes
SE64_9	Tpe_SPC	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Seed	NA	Seed	T.spelta	Seed_T.spelta	6n	Europe	Ancestral	No	No	No	Yes
SPB	Tpe_SP1	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Seed	NA	Seed	T.spelta	Seed_T.spelta	6n	Europe	Ancestral	No	Yes	No	No
SPM20Q	Tmo_SP	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Seed	NA	Seed	T.monococcum	Seed_T.monococcum	2n	Europe	Ancestral	No	Yes	No	No
SPT	Ttu_SP	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Seed	NA	Seed	T.turgidum	Seed_T.turgidum	4n	Europe	Ancestral	No	Yes	No	Yes
SPT_RO	SPT_RO	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Field	NA	Root	T.turgidum	Root_T.turgidum	4n	Europe	Ancestral	No	No	No	Yes
SPT_SH	SPT_SH	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Field	NA	Shoot	T.turgidum	Shoot_T.turgidum	4n	Europe	Ancestral	No	No	No	Yes
SPT_SO	SPT_SO	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Field	NA	Soil	T.turgidum	Soil_T.turgidum	4n	Europe	Ancestral	No	No	No	Yes
SPT_SP	SPT_SP	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Field	NA	Spike	T.turgidum	Spike_T.turgidum	4n	Europe	Ancestral	No	No	No	Yes
TBR	Tae_SP1	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Seed	NA	Seed	T.aestivum	Seed_T.aestivum	6n	Europe	Commercial	No	Yes	No	No
TBS	Tae_SP2	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Seed	NA	Seed	T.aestivum	Seed_T.aestivum	6n	Europe	Commercial	No	Yes	Yes	No
TDR	Tdu_SP	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Seed	NA	Seed	T.durum	Seed_T.durum	4n	Europe	Commercial	No	Yes	No	Yes
TDR_RO	TDR_RO	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Field	NA	Root	T.durum	Root_T.durum	4n	Europe	Commercial	No	No	No	Yes
TDR_SH	TDR_SH	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Field	NA	Shoot	T.durum	Shoot_T.durum	4n	Europe	Commercial	No	No	No	Yes
TDR_SP	TDR_SP	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Field	NA	Spike	T.durum	Spike_T.durum	4n	Europe	Commercial	No	No	No	Yes
TKS	Tae_GE	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Seed	NA	Seed	T.aestivum	Seed_T.aestivum	6n	Europe	Commercial	No	Yes	No	No
TRI_10765	Tae_GR	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Seed	NA	Seed	T.aestivum	Seed_T.aestivum	6n	Europe	Commercial	No	Yes	No	No
TRI_1998	Tmo_BU	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Seed	NA	Seed	T.monococcum	Seed_T.monococcum	2n	Europe	Ancestral	No	Yes	No	No
TRI_4101	Tdu_AF	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Seed	NA	Seed	T.durum	Seed_T.durum	4n	Asia	Commercial	No	Yes	No	No
TRI_4168	Tdi_UN	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Seed	NA	Seed	T.dicoccum	Seed_T.dicoccum	4n	Unknown	Ancestral	No	Yes	No	No
TRI_4323	Tmo_UN	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Seed	NA	Seed	T.monococcum	Seed_T.monococcum	2n	Unknown	Ancestral	No	Yes	No	No
TRI_445	Tdi_US	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Seed	NA	Seed	T.dicoccum	Seed_T.dicoccum	4n	America	Ancestral	No	Yes	No	No
TRI_448	Tdu_AR	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Seed	NA	Seed	T.durum	Seed_T.durum	4n	America	Commercial	No	Yes	No	No
TRI_6263	Tdu_IR	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Seed	NA	Seed	T.durum	Seed_T.durum	4n	Asia	Commercial	No	Yes	No	No
TRI_7121	Tae_NZ	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Seed	NA	Seed	T.aestivum	Seed_T.aestivum	6n	Oceania	Commercial	No	Yes	No	No
EOF

    micromamba run q2-ampl-2025_4 \
        qiime metadata tabulate \
            --m-input-file "${CONF_FOLDER}/${METADATA_FILE}" \
            --o-visualization "${DATA_FOLDER}/sample-metadata.qzv"

    echo -e '[X] Creating sample-metadata.tsv file\n'
}

check_raw_reads_quality(){
    echo '[ ] Checking raw reads quality'

    RUN_FILES=`tail -n+2 "${CONF_FOLDER}/${MANIFEST_FILE}" | grep -v '^#' | cut -d',' -f2`

    micromamba run fastqc \
        fastqc -t ${THREADS} -q --noextract \
            -o "${RAW_FASTQC_FOLDER}" ${RUN_FILES}
    micromamba run fastqc \
        multiqc -o "${RAW_FASTQC_FOLDER}" "${RAW_FASTQC_FOLDER}"

    echo -e '[X] Checking raw reads quality\n'
}

check_primers_position(){
    echo '[ ] Checking primer position distribution'

    # Forward primers:
    RUN_FILES=`grep 'forward' "${CONF_FOLDER}/${MANIFEST_FILE}" | grep -v '^#' | cut -d',' -f2`

    PRIMER=`dna_to_regexp "${FORWARD_PRIMER}"`
    for f in ${RUN_FILES} ; do
        echo -e '\n###############'
        echo "File: ${f}"
        RUN_READS=`zgrep '^+' "${f}" | wc -l`
        echo " ${RUN_READS} -> Total reads of file"

        # Check primer position in reads using basic search
        zcat "${f}" | \
            awk -v primer="${PRIMER}" '(NR % 4 == 2) {idx = match($0, primer); if (idx != 0) {print idx}}' | \
            sort | uniq -c | sort -nk2

        # Print the first bases for reads with no primer found
        echo -e "\nPrimer: ${FORWARD_PRIMER}"
        zcat "${f}" | \
            awk -v primer="${PRIMER}" -v primer_len="${FORWARD_PRIMER_LEN}" '(NR % 4 == 2) {idx = match($0, primer); if (idx == 0) {print substr($0, 1, primer_len)}}' | \
            sort | uniq -c | sort -nrk1 | head -n10
    done

    # Reverse primers:
    RUN_FILES=`grep 'reverse' "${CONF_FOLDER}/${MANIFEST_FILE}" | grep -v '^#' | cut -d',' -f2`

    PRIMER=`dna_to_regexp "${REVERSE_PRIMER}"`
    for f in ${RUN_FILES} ; do
        echo -e '\n###############'
        echo "File: ${f}"
        RUN_READS=`zgrep '^+' "${f}" | wc -l`
        echo " ${RUN_READS} -> Total reads of file"

        # Check primer position in reads using basic search
        zcat "${f}" | \
            awk -v primer="${PRIMER}" '(NR % 4 == 2) {idx = match($0, primer); if (idx != 0) {print idx}}' | \
            sort | uniq -c | sort -nk2

        # Print the first bases for reads with no primer found
        echo -e "\nPrimer: ${REVERSE_PRIMER}"
        zcat "${f}" | \
            awk -v primer="${PRIMER}" -v primer_len="${REVERSE_PRIMER_LEN}" '(NR % 4 == 2) {idx = match($0, primer); if (idx == 0) {print substr($0, 1, primer_len)}}' | \
            sort | uniq -c | sort -nrk1 | head -n10
    done

    echo -e '[X] Checking primer position distribution\n'
}

import_raw_reads(){
    echo '[ ] Importing raw reads to Qiime2'

    micromamba run q2-ampl-2025_4 \
        qiime tools import \
            --input-format PairedEndFastqManifestPhred33 \
            --input-path "${CONF_FOLDER}/${MANIFEST_FILE}" \
            --type 'SampleData[PairedEndSequencesWithQuality]' \
            --output-path "${DATA_FOLDER}/paired-end-demux-seqs.qza"

    micromamba run q2-ampl-2025_4 \
        qiime demux summarize \
            --i-data "${DATA_FOLDER}/paired-end-demux-seqs.qza" \
            --o-visualization "${DATA_FOLDER}/paired-end-demux-seqs.qzv"

    echo -e '[X] Importing raw reads to Qiime2\n'
}

trim_primers_from_reads(){
    echo '[ ] Trimming primers from raw reads'

    micromamba run q2-ampl-2025_4 \
        qiime cutadapt trim-paired \
            --i-demultiplexed-sequences "${DATA_FOLDER}/paired-end-demux-seqs.qza" \
            --p-adapter-f ^"${FORWARD_PRIMER}"...`reverse_complement_dna "${REVERSE_PRIMER}"` \
            --p-adapter-r ^"${REVERSE_PRIMER}"...`reverse_complement_dna "${FORWARD_PRIMER}"` \
            --p-error-rate ${CUTADAPT_ERROR_RATE} \
            --p-overlap 3 \
            --p-minimum-length ${CUTADAPT_MIN_LENGTH} \
            --p-discard-untrimmed \
            --p-cores ${THREADS} \
            --o-trimmed-sequences "${DATA_FOLDER}/trimmed-seqs.qza" \
            --verbose > "${DATA_FOLDER}/trimmed-seqs.log"

    FORWARD_TRIMS=`grep "Sequence: ${FORWARD_PRIMER}" "${DATA_FOLDER}/trimmed-seqs.log" | \
        cut -d' ' -f9 | paste -d+ -s | bc`
    echo "Forward primer trimmed ${FORWARD_TRIMS} times"
    REVERSE_TRIMS=`grep "Sequence: ${REVERSE_PRIMER}" "${DATA_FOLDER}/trimmed-seqs.log" | \
        cut -d' ' -f9 | paste -d+ -s | bc`
    echo "Reverse primer trimmed ${REVERSE_TRIMS} times"

    micromamba run q2-ampl-2025_4 \
        qiime demux summarize \
            --i-data "${DATA_FOLDER}/trimmed-seqs.qza" \
            --o-visualization "${DATA_FOLDER}/trimmed-seqs.qzv"

    micromamba run q2-ampl-2025_4 \
        qiime tools export \
            --input-path "${DATA_FOLDER}/trimmed-seqs.qza" \
            --output-path "${DATA_FOLDER}/trimmed-seqs"

    echo -e '[X] Trimming primers from raw reads\n'
}

check_trimmed_reads_quality(){
    echo '[ ] Checking trimmed reads quality'

    RUN_FILES=`find "${DATA_FOLDER}/trimmed-seqs" -type f -name '*.fastq.gz'`

    micromamba run fastqc \
        fastqc -t ${THREADS} -q --noextract \
            -o "${TRIMMED_FASTQC_FOLDER}" ${RUN_FILES}
    micromamba run fastqc \
        multiqc -o "${TRIMMED_FASTQC_FOLDER}" "${TRIMMED_FASTQC_FOLDER}"

    echo -e '[X] Checking trimmed reads quality\n'
}

dada2_denoise_reads(){
    echo '[ ] Denoising reads with DADA2'

    micromamba run q2-ampl-2025_4 \
        qiime dada2 denoise-paired \
            --i-demultiplexed-seqs "${DATA_FOLDER}/trimmed-seqs.qza" \
            --p-trunc-len-f ${DADA_TRUNC_LEN_F} \
            --p-trunc-len-r ${DADA_TRUNC_LEN_R} \
            --p-trim-left-f ${DADA_TRIM_LEFT_F} \
            --p-trim-left-r ${DADA_TRIM_LEFT_R} \
            --p-max-ee-f ${DADA_MAX_EE_F} \
            --p-max-ee-r ${DADA_MAX_EE_R} \
            --p-trunc-q ${DADA_TRUNC_Q} \
            --p-chimera-method consensus \
            --p-n-threads ${THREADS} \
            --o-representative-sequences "${DATA_FOLDER}/asv-seqs.qza" \
            --o-table "${DATA_FOLDER}/asv-table.qza" \
            --o-denoising-stats "${DATA_FOLDER}/denoising-stats.qza"

    micromamba run q2-ampl-2025_4 \
        qiime metadata tabulate \
            --m-input-file "${DATA_FOLDER}/denoising-stats.qza" \
            --o-visualization "${DATA_FOLDER}/denoising-stats.qzv"

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

    echo -e '[X] Denoising reads with DADA2\n'
}

dna_to_regexp(){
    echo "$1" | \
        sed 's/W/[AT]/' | sed 's/w/[at]/' | \
        sed 's/S/[CG]/' | sed 's/s/[cg]/' | \
        sed 's/M/[AC]/' | sed 's/m/[ac]/' | \
        sed 's/K/[GT]/' | sed 's/k/[gt]/' | \
        sed 's/R/[AG]/' | sed 's/r/[ag]/' | \
        sed 's/Y/[CT]/' | sed 's/y/[ct]/' | \
        sed 's/B/[CGT]/' | sed 's/b/[cgt]/' | \
        sed 's/D/[AGT]/' | sed 's/d/[agt]/' | \
        sed 's/H/[ACT]/' | sed 's/h/[act]/' | \
        sed 's/V/[ACG]/' | sed 's/v/[acg]/' | \
        sed 's/N/[ACGT]/' | sed 's/n/[acgt]/'
}

reverse_complement_dna(){
    echo "$1" | rev | \
        tr "[ACGTWSMKRYBDHVNacgtwsmkrybdhvn]" "[TGCAWSKMYRVHDBNtgcawskmyrvhdbn]"
}

main
