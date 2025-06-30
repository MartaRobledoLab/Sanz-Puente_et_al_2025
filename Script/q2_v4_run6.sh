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

    RUNDATE='2025-04-07'
    RAW_READS_FOLDER="${BASE_FOLDER}/run_250407"

    CONF_FOLDER="${BASE_FOLDER}/etc_v4_run6"
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
    # Conf length = 235 + 179 - 12 (overlap) = 402 >= 251 (min required)
    DADA_TRUNC_LEN_F=235 # Q1 >= 34 (except at 2 spots)
    DADA_TRUNC_LEN_R=179 # Q1 >= 34
    DADA_TRIM_LEFT_F=0
    DADA_TRIM_LEFT_R=0
    DADA_MAX_EE_F=2.0
    DADA_MAX_EE_R=2.0
    DADA_TRUNC_Q=2

    DATA_FOLDER="${BASE_FOLDER}/data_v4_run6"
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
    echo '[ ] Creating sample-metadata.tsv file'

    cat > "${CONF_FOLDER}/${METADATA_FILE}" <<EOF
sample-id	SampleName	ForwardPrimer	ReversePrimer	PrimerName	Region	RunDate	PNAs	Beads	Enzyme	SampleType	ControlType	PlantPart	WheatSpecies	PlantPart_WheatSpecies	PloidyLevel	Location	Domestication	Control	Seed	Gnotobiotic	Field
#q2:types	categorical	categorical	categorical	categorical	categorical	categorical	categorical	categorical	categorical	categorical	categorical	categorical	categorical	categorical	categorical	categorical	categorical	categorical	categorical	categorical	categorical
Negative_Control	Negative_Control	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Control	H2O	NA	NA	NA	NA	NA	NA	Yes	No	No	No
RO_Ex	RO_Ex	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Control	Extraction	NA	NA	NA	NA	NA	NA	Yes	No	No	No
SP_Ex	SP_Ex	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Control	Extraction	NA	NA	NA	NA	NA	NA	Yes	No	No	No
Sto_Bag	Sto_Bag	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Control	Sto_Bag	NA	NA	NA	NA	NA	NA	Yes	No	No	No
Sto_Bag_2	Sto_Bag_2	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Control	Sto_Bag	NA	NA	NA	NA	NA	NA	Yes	No	No	No
Tae_RO_2	Tae_RO_2	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Field	NA	Root	T.aestivum	Root_T.aestivum	6n	America	Commercial	No	No	No	Yes
Tae_SE_1	Tae_SE_1	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Seed	NA	Seed	T.aestivum	Seed_T.aestivum	6n	Europe	Commercial	No	No	Yes	Yes
Tae_SE_5	Tae_SE_5	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Seed	NA	Seed	T.aestivum	Seed_T.aestivum	6n	Europe	Commercial	No	No	No	Yes
Tae_SH_2	Tae_SH_2	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Field	NA	Shoot	T.aestivum	Shoot_T.aestivum	6n	Europe	Commercial	No	No	No	Yes
Tae_SO_1	Tae_SO_1	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Field	NA	Soil	T.aestivum	Soil_T.aestivum	6n	Europe	Commercial	No	No	No	Yes
Tae_SO_2	Tae_SO_2	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Field	NA	Soil	T.aestivum	Soil_T.aestivum	6n	Europe	Commercial	No	No	No	Yes
Tae_SP_1	Tae_SP_1	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Field	NA	Spike	T.aestivum	Spike_T.aestivum	6n	Europe	Commercial	No	No	No	Yes
Tdu_RO_2	Tdu_RO_2	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Field	NA	Root	T.durum	Root_T.durum	4n	Europe	Commercial	No	No	No	Yes
Tdu_RO_3	Tdu_RO_3	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Field	NA	Root	T.durum	Root_T.durum	4n	Europe	Commercial	No	No	No	Yes
Tdu_SE_1	Tdu_SE_1	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Seed	NA	Seed	T.durum	Seed_T.durum	4n	Europe	Commercial	No	No	No	Yes
Tdu_SE_3	Tdu_SE_3	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Seed	NA	Seed	T.durum	Seed_T.durum	4n	Europe	Commercial	No	No	No	Yes
Tdu_SH_1	Tdu_SH_1	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Field	NA	Shoot	T.durum	Shoot_T.durum	4n	Europe	Commercial	No	No	No	Yes
Tdu_SH_2	Tdu_SH_2	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Field	NA	Shoot	T.durum	Shoot_T.durum	4n	Europe	Commercial	No	No	No	Yes
Tdu_SO_1	Tdu_SO_1	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Field	NA	Soil	T.durum	Soil_T.durum	4n	Europe	Commercial	No	No	No	Yes
Tdu_SO_2	Tdu_SO_2	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Field	NA	Soil	T.durum	Soil_T.durum	4n	Europe	Commercial	No	No	No	Yes
Tdu_SO_3	Tdu_SO_3	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Field	NA	Soil	T.durum	Soil_T.durum	4n	Europe	Commercial	No	No	No	Yes
Tdu_SP_1	Tdu_SP_1	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Field	NA	Spike	T.durum	Spike_T.durum	4n	Europe	Commercial	No	No	No	Yes
Tdu_SP_3	Tdu_SP_3	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Field	NA	Spike	T.durum	Spike_T.durum	4n	Europe	Commercial	No	No	No	Yes
Tpe_RO_1	Tpe_RO_1	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Field	NA	Root	T.spelta	Root_T.spelta	6n	Europe	Ancestral	No	No	No	Yes
Tpe_RO_2	Tpe_RO_2	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Field	NA	Root	T.spelta	Root_T.spelta	6n	Europe	Ancestral	No	No	No	Yes
Tpe_SH_1	Tpe_SH_1	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Field	NA	Shoot	T.spelta	Shoot_T.spelta	6n	Europe	Ancestral	No	No	No	Yes
Tpe_SH_2	Tpe_SH_2	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Field	NA	Shoot	T.spelta	Shoot_T.spelta	6n	Europe	Ancestral	No	No	No	Yes
Tpe_SO_1	Tpe_SO_1	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Field	NA	Soil	T.spelta	Soil_T.spelta	6n	Europe	Ancestral	No	No	No	Yes
Tpe_SO_2	Tpe_SO_2	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Field	NA	Soil	T.spelta	Soil_T.spelta	6n	Europe	Ancestral	No	No	No	Yes
Tpe_SP_1	Tpe_SP_1	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Field	NA	Spike	T.spelta	Spike_T.spelta	6n	Europe	Ancestral	No	No	No	Yes
Tpe_SP_2	Tpe_SP_2	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Field	NA	Spike	T.spelta	Spike_T.spelta	6n	Europe	Ancestral	No	No	No	Yes
Ttu_RO_1	Ttu_RO_1	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Field	NA	Root	T.turgidum	Root_T.turgidum	4n	Europe	Ancestral	No	No	No	Yes
Ttu_RO_2	Ttu_RO_2	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Field	NA	Root	T.turgidum	Root_T.turgidum	4n	Europe	Ancestral	No	No	No	Yes
Ttu_SE_1	Ttu_SE_1	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Seed	NA	Seed	T.turgidum	Seed_T.turgidum	4n	Europe	Ancestral	No	No	No	Yes
Ttu_SE_2	Ttu_SE_2	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Seed	NA	Seed	T.turgidum	Seed_T.turgidum	4n	Europe	Ancestral	No	No	No	Yes
Ttu_SH_1	Ttu_SH_1	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Field	NA	Shoot	T.turgidum	Shoot_T.turgidum	4n	Europe	Ancestral	No	No	No	Yes
Ttu_SH_2	Ttu_SH_2	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Field	NA	Shoot	T.turgidum	Shoot_T.turgidum	4n	Europe	Ancestral	No	No	No	Yes
Ttu_SO_1	Ttu_SO_1	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Field	NA	Soil	T.turgidum	Soil_T.turgidum	4n	Europe	Ancestral	No	No	No	Yes
Ttu_SO_2	Ttu_SO_2	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Field	NA	Soil	T.turgidum	Soil_T.turgidum	4n	Europe	Ancestral	No	No	No	Yes
Ttu_SP_1	Ttu_SP_1	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Field	NA	Spike	T.turgidum	Spike_T.turgidum	4n	Europe	Ancestral	No	No	No	Yes
Ttu_SP_2	Ttu_SP_2	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Field	NA	Spike	T.turgidum	Spike_T.turgidum	4n	Europe	Ancestral	No	No	No	Yes
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
