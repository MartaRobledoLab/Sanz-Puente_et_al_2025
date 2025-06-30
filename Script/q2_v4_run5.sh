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

    RUNDATE='2023-06-06'
    RAW_READS_FOLDER="${BASE_FOLDER}/run_230606"

    CONF_FOLDER="${BASE_FOLDER}/etc_v4_run5"
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
    # Conf length = 213 + 135 - 12 (overlap) = 336 >= 251 (min required)
    DADA_TRUNC_LEN_F=213 # Q1 >= 34
    DADA_TRUNC_LEN_R=135 # Q1 >= 34 (except at 6 spots)
    DADA_TRIM_LEFT_F=0
    DADA_TRIM_LEFT_R=0
    DADA_MAX_EE_F=2.0
    DADA_MAX_EE_R=2.0
    DADA_TRUNC_Q=2

    DATA_FOLDER="${BASE_FOLDER}/data_v4_run5"
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
EOF

    echo -e '[X] Creating sample-manifest.csv file\n'
}

create_metadata_file(){
    echo '[ ] Creating sample-metadata.tsv file'

    cat > "${CONF_FOLDER}/${METADATA_FILE}" <<EOF
sample-id	SampleName	ForwardPrimer	ReversePrimer	PrimerName	Region	RunDate	PNAs	Beads	Enzyme	SampleType	ControlType	PlantPart	WheatSpecies	PlantPart_WheatSpecies	PloidyLevel	Location	Domestication	Control	Seed	Gnotobiotic	Field
#q2:types	categorical	categorical	categorical	categorical	categorical	categorical	categorical	categorical	categorical	categorical	categorical	categorical	categorical	categorical	categorical	categorical	categorical	categorical	categorical	categorical	categorical
A1	Tae_Sh1	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Gnotobiotic	NA	Shoot	T.aestivum	Shoot_T.aestivum	6n	Europe	Commercial	No	No	Yes	No
A2	Tae_Sh2	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Gnotobiotic	NA	Shoot	T.aestivum	Shoot_T.aestivum	6n	Europe	Commercial	No	No	Yes	No
A3	Tae_Sh3	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Gnotobiotic	NA	Shoot	T.aestivum	Shoot_T.aestivum	6n	Europe	Commercial	No	No	Yes	No
R1	Tae_Ro1	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Gnotobiotic	NA	Root	T.aestivum	Root_T.aestivum	6n	Europe	Commercial	No	No	Yes	No
R2	Tae_Ro2	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Gnotobiotic	NA	Root	T.aestivum	Root_T.aestivum	6n	Europe	Commercial	No	No	Yes	No
R3	Tae_Ro3	${FORWARD_PRIMER}	${REVERSE_PRIMER}	${PRIMER_NAME}	${REGION}	${RUNDATE}	No	Yes	Kappa	Gnotobiotic	NA	Root	T.aestivum	Root_T.aestivum	6n	Europe	Commercial	No	No	Yes	No
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
