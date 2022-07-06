#!/usr/bin/env bash


L_CONSENSUS_DIR="/Users/jjuma/Documents/PhD_RVF2019/pipelines/rvfvampliconseq/outdir/L-segment/ivar/consensus"
#L_CONSENSUS_DIR="/Users/jjuma/Documents/PhD_RVF2019/pipelines/rvfvampliconseq/data/genomics/rvfv/genome/"
#M_CONSENSUS_DIR="/Users/jjuma/Documents/PhD_RVF2019/pipelines/rvfvampliconseq/outdir/M-segment/ivar/consensus"
#M_CONSENSUS_DIR="/Users/jjuma/Documents/PhD_RVF2019/pipelines/rvfvampliconseq/data/genomics/rvfv/genome/"
#S_CONSENSUS_DIR="/Users/jjuma/Documents/PhD_RVF2019/pipelines/rvfvampliconseq/outdir/S-segment/ivar/consensus"
#S_CONSENSUS_DIR="/Users/jjuma/Documents/PhD_RVF2019/pipelines/rvfvampliconseq/data/genomics/rvfv/genome/"


# ncbi sequences
L_BED="/Users/jjuma/Documents/PhD_RVF2019/projects/rvfv-primers/RVFV_L_400.primer.tsv"
#M_BED="/Users/jjuma/Documents/PhD_RVF2019/projects/rvfv-primers/RVFV_M_400.primer.tsv"
#S_BED="/Users/jjuma/Documents/PhD_RVF2019/projects/rvfv-primers/RVFV_S_400.primer.tsv"

# vipr sequences
# L_BED="/Users/jjuma/Documents/PhD_RVF2019/projects/primalscheme_rvfv/primalscheme_output/vipr/rvfv-L.primer.tsv"
# M_BED="/Users/jjuma/Documents/PhD_RVF2019/projects/primalscheme_rvfv/primalscheme_output/vipr/rvfv-M.primer.tsv"
# S_BED="/Users/jjuma/Documents/PhD_RVF2019/projects/primalscheme_rvfv/primalscheme_output/vipr/rvfv-S.primer.tsv"

# sinngle reference
# S_BED="/Users/jjuma/Documents/PhD_RVF2019/projects/primalscheme_rvfv/primalscheme_output/web_One_reference/rvfv-S-400.primer.tsv"
# M_BED="/Users/jjuma/Documents/PhD_RVF2019/projects/primalscheme_rvfv/primalscheme_output/web_One_reference/rvfv-M-300.primer.tsv"
# L_BED="/Users/jjuma/Documents/PhD_RVF2019/projects/primalscheme_rvfv/primalscheme_output/web_One_reference/rvfv-L-350.primer.tsv"

BIN_DIR="/Users/jjuma/Documents/PhD_RVF2019/pipelines/rvfvampliconseq/bin"



cd "${L_CONSENSUS_DIR}"
# find primer positions in each consensus fasta file
for file in *.fa;
 do
   echo -e "getting primer positions in: " "$file";
   sample_id="${file%.fa}"
   [[ -e "$file" ]] || break
   python ${BIN_DIR}/find_amplicons.py --fasta $file --bed $L_BED
done


# cd "${M_CONSENSUS_DIR}"
# # find primer positions in each consensus fasta file
# for file in *.fasta;
#  do
#    echo -e "getting primer positions in: " "$file";
#    sample_id="${file%.fasta}"
#    [[ -e "$file" ]] || break
#    python ${BIN_DIR}/find_amplicons.py --fasta $file --bed $M_BED
# done

# cd "${S_CONSENSUS_DIR}"
# # find primer positions in each consensus fasta file
# for file in *.fasta;
#   do
#     echo -e "getting primer positions in: " "$file";
#     sample_id="${file%.fasta}"
#     [[ -e "$file" ]] || break
#     python ${BIN_DIR}/find_amplicons.py --fasta $file --bed $S_BED
# done