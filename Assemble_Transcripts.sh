#!/usr/bin/env bash

VERSION="0.1.1"

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
  echo -e "[`date`]\tCMD: $CMD"
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
CLEAN_START=0
SRA_IDS=""
SS="rf"
OUT=""
NCPUS=24
MEM=120

## Useage information
usage() {
echo -e "##
## ${SCRIPT_NAME} v${VERSION}
##

Takes a list of SRA Run IDs and downloads, trimms, and assembles them into a single combined transcriptome.

Usage: 
# Single SRA Run
./${SCRIPT_NAME} -s DRR330246 -o test/single
# Multiple SRA Runs
./${SCRIPT_NAME} -s DRR330244,DRR330245,DRR330246 -o test/multiple
# Set number of cpus and memory to use (mainly affects assembly)
./${SCRIPT_NAME} -s DRR330244,DRR330245,DRR330246 -o test/multiple -m 60 -n 48


Options (all optional):
-s, --sra                  SRA ID to download and assemble (required).
                             - Multiple (upto 9) SRA IDs can be provided as a comma separated list.
-o, --out                  Output name (required)
--ss                       Strandedness of samples (default: ${SS})
-n, --ncpus                Max number of threads to use for sample download, trimming, and assembly (default: ${NCPUS})
-m, --mem                  Max memory (in gigabytes) for spades to use (default: ${MEM})

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
    -s|--sra)
      SRA_IDS="$2"
      shift # past argument
      shift # past value
      ;;
    -o|--out)
      OUT="$2"
      shift # past argument
      shift # past value
      ;;
    --ss)
      SS="$2"
      shift # past argument
      shift # past value
      ;;
    -n|--ncpus)
      NCPUS="$2"
      shift # past argument
      shift # past value
      ;;
    -m|--mem)
      MEM="$2"
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
if [[ -z "${SRA_IDS}" || -z "${OUT}" ]]; then
  err "\n\n\t*** ERROR: Not all required values provided ***\n\n"
  usage
fi
set -eu


## Check execuatble
FAIL_COUNT=0

KINGFISHER=$(which kingfisher || echo "")
if [[ -z ${KINGFISHER} ]]; then
  err "kingfisher missing from PATH."
  FAIL_COUNT=$((FAIL_COUNT+1))
fi

FASTP=$(which fastp || echo "")
if [[ -z ${FASTP} ]]; then
  err "fastp missing from PATH."
  FAIL_COUNT=$((FAIL_COUNT+1))
fi

RNASPADES=$(which rnaspades.py || echo "")
if [[ -z ${RNASPADES} ]]; then
  err "rnaspades.py missing from PATH."
  FAIL_COUNT=$((FAIL_COUNT+1))
fi

STATS=$(which stats.sh || echo "")
if [[ -z ${STATS} ]]; then
  err "stats.sh from the bbmap package missing from PATH."
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


## Clean up existing output and create new working dir
if [[ $CLEAN_START -eq 1 ]];
then
  log "Removing any old analysis before we start"
  rm -fr "${OUT}" "${OUT}.transcripts.fasta" "${OUT}.transcripts.fasta.stats"
fi
mkdir -p "${OUT}"


################################################
################                ################
################ START ANALYSIS ################
################                ################
################################################

SPADES_READS="" # Store all the reads that we want to run spades on
SRA_COUNT=0 # Count how many samples we have processed

# Set max download NCPUS to 16
#   - 'aria2c' throws an error if NCPUS > 16
#   - 'fastp' only uses a max of 16 threads
if [[ $NCPUS -gt 16 ]]; then tNCPUS=16; else tNCPUS=$NCPUS; fi

while read SRA;
do
  SRA_COUNT=$((SRA_COUNT+1))
  log "Processing $SRA"
  
  ##
  ## Download
  ##
  CHECKPOINT="${OUT}/${SRA}.download.done"
  
  # Set max download NCPUS to 16 - 'aria2c' throws an error if NCPUS > 16
  if [[ $NCPUS -gt 16 ]]; then tNCPUS=16; else tNCPUS=$NCPUS; fi
  
  CMD="(cd ${OUT}; ${KINGFISHER} get -r ${SRA} -f fastq.gz --check-md5sums -m ena-ftp aws-http prefetch aws-cp --download-threads ${tNCPUS} --extraction-threads ${tNCPUS})"
  CMD="${CMD} && touch ${CHECKPOINT}"
  
  if [[ ! -f $CHECKPOINT ]];
  then
    log "  - Downloading"
    run_cmd "${CMD}"
    log "    - Done."
  else
    log "    - Already done."
  fi
  
  
  ##
  ## Trim
  ##
  CHECKPOINT="${OUT}/${SRA}.trim.done"
  
  if [[ -f "${OUT}/${SRA}_2.fastq.gz" || -f "${OUT}/${SRA}_2.trimmed.fastq.gz" ]];
  then
    I1="${OUT}/${SRA}_1.fastq.gz"
    I2="${OUT}/${SRA}_2.fastq.gz"
    O1="${OUT}/${SRA}_1.trimmed.fastq.gz"
    O2="${OUT}/${SRA}_2.trimmed.fastq.gz"
    CMD="fastp --in1 ${I1} --in2 ${I2} --out1 ${O1} --out2 ${O2} --json ${OUT}/${SRA}_fastp.json --html ${OUT}/${SRA}_fastp.html --thread ${tNCPUS}"
    CMD="${CMD} && rm -fr ${I1} ${I2}"
    CMD="${CMD} && touch ${CHECKPOINT}"
    SPADES_READS="${SPADES_READS} --pe${SRA_COUNT}-1 ${O1} --pe${SRA_COUNT}-2 ${O2}"
  else
    I1="${OUT}/${SRA}.fastq.gz"
    O1="${OUT}/${SRA}_1.trimmed.fastq.gz"
    CMD="fastp --in1 ${I1} --out1 ${O1} --json ${OUT}/${SRA}_fastp.json --html ${OUT}/${SRA}_fastp.html --thread ${tNCPUS}"
    CMD="${CMD} && rm -fr ${I1}"
    CMD="${CMD} && touch ${CHECKPOINT}"
    SPADES_READS="${SPADES_READS} --s ${SRA_COUNT} ${O1}"
  fi
  
  if [[ ! -f $CHECKPOINT ]];
  then
    log "  - Trimming"
    run_cmd "${CMD}"
    log "    - Done."
  else
    log "    - Already done."
  fi
  
done < <(echo "${SRA_IDS}" | sed -e 's/,/\n/g')


##
## Assemble transcripts
##
CHECKPOINT="${OUT}/spades.done"

if [[ ! -d "${OUT}/spades" ]];
then
  CMD="rnaspades.py --threads ${NCPUS} --memory ${MEM} -o ${OUT}/spades --ss ${SS} ${SPADES_READS}"
else
  CMD="rnaspades.py -o ${OUT}/spades --continue"
fi
CMD="${CMD} && touch ${CHECKPOINT}"

if [[ ! -f $CHECKPOINT ]];
then
  log "  - Asssembling"
  run_cmd "${CMD}"
  log "    - Done."
else
  log "    - Already done."
fi


##
## Copy and cleanup assembly results
##
CHECKPOINT="${OUT}/rename.done"

CMD="cat ${OUT}/spades/transcripts.fasta | sed -e 's@>@>${OUT##*/}-@' > ${OUT}.transcripts.fasta"
CMD="${CMD} && touch ${CHECKPOINT}"

if [[ ! -f $CHECKPOINT ]];
then
  log "  - Renaming"
  run_cmd "${CMD}"
  log "    - Done."
else
  log "    - Already done."
fi


##
## Assembly stats
##
CHECKPOINT="${OUT}/stats.done"

CMD="stats.sh ${OUT}.transcripts.fasta > ${OUT}.transcripts.fasta.stats"
CMD="${CMD} && touch ${CHECKPOINT}"

if [[ ! -f $CHECKPOINT ]];
then
  log "  - Stats"
  run_cmd "${CMD}"
  log "    - Done."
  cat "${OUT}.transcripts.fasta.stats"
else
  log "    - Already done."
fi


##
## Finsihed
##
log "Workflow has finished running!"


###################################################
################                   ################
################ FINISHED ANALYSIS ################
################                   ################
###################################################
