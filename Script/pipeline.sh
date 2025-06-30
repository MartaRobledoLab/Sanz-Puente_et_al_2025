#!/usr/bin/env bash

if true ; then

./q2_v4_run4.sh 2>&1 | tee q2_v4.log
./q2_v4_run5.sh 2>&1 | tee -a q2_v4.log
./q2_v4_run6.sh 2>&1 | tee -a q2_v4.log

fi

if true ; then

./q2_v4_merg.sh 2>&1 | tee -a q2_v4.log

fi

if true ; then

./Metabarcoding_analysis.R
./qPCR.R

fi
