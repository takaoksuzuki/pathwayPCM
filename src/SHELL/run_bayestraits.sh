#PBS -S /bin/bash
#PBS -l ncpus=1
#PBS -M suzukitakaotakao@gmail.com
#PBS -m abe

cd /data/takaoksuzuki/BayesTraits/20230504_GTDB207_cellshape_motility_sporulation_ALL/

./BayesTraitsV4 Proteobacteria_cellshape_motility_sporulation_bac120_r207.tree ch_GTDB207_cellshape_motility_sporulation_Proteobacteria.tsv < bac_1_r207_Proteobacteria_cellshape_motility_sporulation_inputcommands_Res_g0101_a.txt > fout_a_cellshape_motility_sporulation_Proteobacteria_rj_g0101.txt
