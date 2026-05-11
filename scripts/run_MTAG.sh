#!/usr/bin/env bash

MTAG_SCRIPT=$1
SUMSTATS=$2
OUT=$3
SNP=$4
A1=$5
A2=$6
Z=$7
N=$8
MAF_MIN=$9

python $MTAG_SCRIPT \
  --sumstats $SUMSTATS \
  --out $OUT \
  --snp_name $SNP \
  --a1_name $A1 \
  --a2_name $A2 \
  --z_name $Z \
  --n_name $N \
  --maf_min $MAF_MIN \
  --stream_stdout \
  --force
