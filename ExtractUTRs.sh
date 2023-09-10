#!/usr/bin/env bash

VERSION="0.1"

## Pre-run setup
set -euo pipefail
IFS=$'\n\t'

# Absolute path to this script, e.g. /home/user/bin/foo.sh
SCRIPT=$(readlink -f "$0")
# Absolute path this script is in, thus /home/user/bin
SCRIPTPATH=$(dirname "$SCRIPT")
# Script name
SCRIPT_NAME=$(basename ${SCRIPT})


## Helper functions
function run_cmd(){
local CMD="$@"
  echo -t "[`date`]\tCMD: $CMD"
  eval $CMD
}

function log(){
  echo -e "[`date`]\tLOG: $@"
}

function err(){
  echo -e "[`date`]\tERROR: $@" >&2
}

function check_empty(){
  if [[ ! -s "${1}" ]];
  then
    err "${1} is empty!"
    exit 2
  fi
}


## Set option envs
TRANSCRIPTS=""
DIAMOND_DB=""
NCPUS=6
CLEAN_START=1
MIN_ORF_LENGTH=20
MIN_UTR_LENGTH=10
GENETIC_CODE="Universal"
SUFFIX="ExUTRs"
DIAMOND_OPTS="--max-target-seqs 1 --ultra-sensitive --iterate --outfmt 6 --evalue 1e-5"


## Useage information
usage() {
echo -e "##
## ${SCRIPT_NAME} v${VERSION}
##

Predict ORFs in transcripts and extract 5- and 3-prime UTR regions

Usage: 
./${SCRIPT_NAME} -t transcripts.fasta -d /path/to/rn.dmnd

Options (all optional):
-t, --transcripts          Transcripts fasta file (required)
-d, --database             File with diamond database to use to help guild ORF prediction (optional)
-n, --ncpus                Number of threads to use for DIAMOND (default: ${NCPUS})

-v, --version              Script version (v${VERSION})
-h, --help                 This help message
--debug                    Run debug mode
" 1>&2
exit 1
}


# See https://stackoverflow.com/questions/192249/how-do-i-parse-command-line-arguments-in-bash

POSITIONAL=()
while [[ $# -gt 0 ]]; do
  key="$1"
  case $key in
    -t|--transcripts)
      TRANSCRIPTS="$2"
      shift # past argument
      shift # past value
      ;;
    -d|--database)
      DIAMOND_DB="$2"
      shift # past argument
      shift # past value
      ;;
    -n|--ncpus)
      NCPUS="$2"
      shift # past argument
      shift # past value
      ;;
    -h|--help)
      usage
      exit 1;
      ;;
    -v|--version)
      echo "v${VERSION}"
      exit 0;
      ;;
    --debug)
      set -x
      shift # past argument
      ;;
    *) # unknown option
      POSITIONAL+=("$1") # save it in an array for later
      shift # past argument
      ;;
  esac
done

set -- "${POSITIONAL[@]}" # restore positional parameters


## Check required args are set
set +eu
if [ -z "${TRANSCRIPTS}" ]; then
  err "\n\n\t*** ERROR: Not all required values provided ***\n\n"
  usage
fi
set -eu


## Check options set
if [[ ! -z $TRANSCRIPTS && ! -s $TRANSCRIPTS ]]; then log "$TRANSCRIPTS does not exist!"; exit 1; fi
if [[ ! -z $DIAMOND_DB  && ! -s $DIAMOND_DB  ]]; then log "$DIAMOND_DB does not exist!"; exit 1; fi


## Check execuatble
FAIL_COUNT=0

TRANSDECODER_LONGESTORFS=$(which TransDecoder.LongOrfs || echo "")
if [[ -z ${TRANSDECODER_LONGESTORFS} ]]; then
  err "TransDecoder.LongOrfs missing from PATH."
  FAIL_COUNT=$((FAIL_COUNT+1))
fi

TRANSDECODER_PREDICT=$(which TransDecoder.Predict || echo "")
if [[ -z ${TRANSDECODER_LONGESTORFS} ]]; then
  err "TransDecoder.Predict missing from PATH."
  FAIL_COUNT=$((FAIL_COUNT+1))
fi

DIAMOND=$(which diamond || echo "")
if [[ -z ${DIAMOND} ]]; then
  err "diamond missing from PATH."
  FAIL_COUNT=$((FAIL_COUNT+1))
fi

SEQKIT=$(which seqkit || echo "")
if [[ -z ${SEQKIT} ]]; then
  err "seqkit missing from PATH."
  FAIL_COUNT=$((FAIL_COUNT+1))
fi

if [ $FAIL_COUNT -gt 0 ];
then
  err "...Aborting!"
  exit 1
fi


## Set envs
OUT="${TRANSCRIPTS}.${SUFFIX}"


## Clean up existing output 
if [[ $CLEAN_START -eq 1 ]];
then
  log "Removing any old analysis before we start"
  rm -fr "${OUT}"* "${TRANSCRIPTS}.transdecoder"* "${TRANSCRIPTS}.seqkit.fai"
fi


################################################
################                ################
################ START ANALYSIS ################
################                ################
################################################

QUERY="${TRANSCRIPTS}.transdecoder_dir/longest_orfs.pep"
DIAMOND_OUT="${QUERY}.diamondp.outfmt6"

## Run TransDecoder.LongOrfs
log "Running TransDecoder.LongOrfs to predict initial ORFs in transcripts"

CMD="TransDecoder.LongOrfs --genetic_code ${GENETIC_CODE} -m ${MIN_ORF_LENGTH} -t ${TRANSCRIPTS}"
run_cmd "${CMD}"


## Run DIAMOND
log "Running DIAMOND BLASTP against the provided database to assist TransDecoder ORF prediction (if a database was given!)"

CMD="diamond blastp --query ${QUERY} --db ${DIAMOND_DB} --out ${DIAMOND_OUT} --threads ${NCPUS} ${DIAMOND_OPTS}"
if [[ ! -z $DIAMOND_DB ]];
then
  run_cmd "${CMD}"
else
  log "  - No database given, skipping this step."
fi


## Run TransDecoder.Predict
log "Running TransDecoder.Predict to produce the final ORFs in transcripts"

CMD="TransDecoder.Predict --genetic_code ${GENETIC_CODE} -t ${TRANSCRIPTS} --single_best_only"
if [[ ! -z $DIAMOND_DB ]];
then
  CMD="${CMD} --retain_blastp_hits ${DIAMOND_OUT}"
fi
run_cmd "${CMD}"


## Extract UTRs
GFF="${TRANSCRIPTS}.transdecoder.gff3"

log "Extracting predicted 3-prime UTRs"
cat "$GFF" \
  | awk -F'\t' 'NF==9 && $3=="three_prime_UTR" {split($9,a,";"); NAME=""; for(i=1; i<=length(a); i++){split(a[i],b,"="); if(b[1]=="ID"){NAME=b[2]}}; print $1"\t"$4-1"\t"$5"\t"NAME"\t"$6"\t"$7}' \
  | awk -F'\t' -vMIN_UTR_LENGTH="${MIN_UTR_LENGTH}" '($3-$2)>=MIN_UTR_LENGTH' \
  > "${OUT}.3UTR.bed"
seqkit subseq -U --line-width 0 --bed "${OUT}.3UTR.bed" "${TRANSCRIPTS}" | awk -F' ' '{ if($1~"^>"){sub(">","",$1); print ">"$2" "$1}else{print} }' > "${OUT}.3UTR.fasta"


log "Extracting predicted 5-prime UTRs"
cat "$GFF" \
  | awk -F'\t' 'NF==9 && $3=="five_prime_UTR" {split($9,a,";"); NAME=""; for(i=1; i<=length(a); i++){split(a[i],b,"="); if(b[1]=="ID"){NAME=b[2]}}; print $1"\t"$4-1"\t"$5"\t"NAME"\t"$6"\t"$7}' \
  | awk -F'\t' -vMIN_UTR_LENGTH="${MIN_UTR_LENGTH}" '($3-$2)>=MIN_UTR_LENGTH' \
  > "${OUT}.5UTR.bed"
seqkit subseq -U --line-width 0 --bed "${OUT}.5UTR.bed" "${TRANSCRIPTS}" | awk -F' ' '{ if($1~"^>"){sub(">","",$1); print ">"$2" "$1}else{print} }' > "${OUT}.5UTR.fasta"


## Print stats
N_3UTR=$(cat "${OUT}.3UTR.bed" | wc -l)
N_5UTR=$(cat "${OUT}.5UTR.bed" | wc -l)

log "Workflow has finished running!\nFound:\n${N_3UTR} 3-prime UTR regions\n${N_5UTR} 5-prime UTR regions"


###################################################
################                   ################
################ FINISHED ANALYSIS ################
################                   ################
###################################################
