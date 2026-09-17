#!/usr/bin/env bash
# Usage: sbatch combine_vmr_dmr_input.sbatch THRESHOLD 10000|30000
#SBATCH --job-name=mvi_vmr_dmr_join
#SBATCH --account=lijia
#SBATCH --partition=cpu
#SBATCH --nodelist=cu03
#SBATCH --cpus-per-task=4
#SBATCH --mem=40G
#SBATCH --time=1-00:00:00
#SBATCH --output=/home/lijia/luozhixiong/IPF_tissue/Scripts/Methylvi/logs/vmr_dmr_join_%j.out
#SBATCH --error=/home/lijia/luozhixiong/IPF_tissue/Scripts/Methylvi/logs/vmr_dmr_join_%j.err

set -Eeuo pipefail
if (( $# != 2 )); then echo "Usage: $0 0.01|0.02|0.05 10000|30000" >&2; exit 2; fi
case "$1" in 0.01|0.02|0.05) ;; *) exit 3;; esac
case "$2" in 10000|30000) ;; *) exit 4;; esac
project=/home/lijia/luozhixiong/IPF_tissue
old="$project/Results/MethylVI_clean6264_10k30k/vmr/var_${1}/features_${2}/input.h5mu"
root="$project/Results/MethylVI_vmr_plus_all_unique_pooled_DMR"
dmr="$root/shared/dmr/all_unique_hypo_DMRs.input.h5mu"
out="$root/vmr_dmr/var_${1}/features_${2}/input.h5mu"
[[ -f "$root/shared/dmr/counts.COMPLETE" ]] || exit 5
cd "$project"
"/home/lijia/luozhixiong/miniconda3/envs/methylvi/bin/python" \
  Scripts/Methylvi/vmr_dmr/02_combine_vmr_and_dmr_inputs.py \
  --vmr-input "$old" --dmr-input "$dmr" --output "$out"
touch "$(dirname "$out")/input.COMPLETE"
