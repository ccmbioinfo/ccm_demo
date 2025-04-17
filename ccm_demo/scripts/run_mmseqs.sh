#! /bin/bash

set -oe pipefail

MMSEQS=$CONDA_PREFIX/bin/mmseqs
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

while getopts 'i:d:o:c' flag; do
  # shellcheck disable=SC2220
  case "${flag}" in
    i) input="${OPTARG}" ;;
    d) database="${OPTARG}" ;;
    o) output="${OPTARG}" ;;
    c) cleanup='true' ;;
  esac
done

$MMSEQS createdb $input query -v 0
$MMSEQS search query $database search tmp --num-iterations 2 -v 0
$MMSEQS align query $database search align -a -v 0
$MMSEQS convertalis query $database align query.tab --format-output target,qlen,qstart,qend,tstart,tend,tseq,cigar,taln -v 0
python $SCRIPT_DIR/process_mmseqs.py -i query.tab -o $output



