#!/bin/bash
#SBATCH --cpus-per-task=8

# cflo vs sinv
#################
# proteinortho5 -step=1 -project=cflo_v_sinv_1 -p=blastp+ ./sinv ./cflo_v7dot5
## done
# proteinortho5 -step=2 -project=cflo_v_sinv_1 -p=blastp+ ./sinv ./cflo_v7dot5
## done

# cflo vs amel
#################
#proteinortho5 -step=1 -project=cflo_v_amel_1 -p=blastp+ ./protein_seqs/amel.faa ./protein_seqs/cflo_v7dot5
#proteinortho5 -step=2 -project=cflo_v_amel_1 -p=blastp+ ./protein_seqs/amel.faa ./protein_seqs/cflo_v7dot5

# cflo vs dmel
#################
#proteinortho5 -step=1 -project=cflo_v_dmel_1 -p=blastp+ ./protein_seqs/dmel.faa ./protein_seqs/cflo_v7dot5
#proteinortho5 -step=2 -project=cflo_v_dmel_1 -p=blastp+ ./protein_seqs/dmel.faa ./protein_seqs/cflo_v7dot5

# cflo vs mmus
#################
#proteinortho5 -step=1 -project=cflo_v_mmus_1 -p=blastp+ ./protein_seqs/mmus ./protein_seqs/cflo_v7dot5
#proteinortho5 -step=2 -project=cflo_v_mmus_1 -p=blastp+ ./protein_seqs/mmus ./protein_seqs/cflo_v7dot5

# cflo vs human
#################
#proteinortho5 -step=1 -project=cflo_v_human_1 -p=blastp+ ./protein_seqs/human ./protein_seqs/cflo_v7dot5
#proteinortho5 -step=2 -project=cflo_v_human_1 -p=blastp+ ./protein_seqs/human ./protein_seqs/cflo_v7dot5


# cflo vs agam
#################
#proteinortho5 -step=1 -project=cflo_v_agam_1 -p=blastp+ ./protein_seqs/agam.faa ./protein_seqs/cflo_v7dot5
#proteinortho5 -step=2 -project=cflo_v_agam_1 -p=blastp+ ./protein_seqs/agam.faa ./protein_seqs/cflo_v7dot5

# cflo vs mono
#################
#proteinortho5 -step=1 -project=cflo_v_mono_1 -p=blastp+ ./protein_seqs/mono.faa ./protein_seqs/cflo_v7dot5
#proteinortho5 -step=2 -project=cflo_v_mono_1 -p=blastp+ ./protein_seqs/mono.faa ./protein_seqs/cflo_v7dot5

# cflo vs pogo
#################
#proteinortho5 -step=1 -project=cflo_v_pogo_1 -p=blastp+ ./protein_seqs/pogo.faa ./protein_seqs/cflo_v7dot5
#proteinortho5 -step=2 -project=cflo_v_pogo_1 -p=blastp+ ./protein_seqs/pogo.faa ./protein_seqs/cflo_v7dot5

# cflo vs temno
#################
#proteinortho5 -step=1 -project=cflo_v_temno_1 -p=blastp+ ./protein_seqs/temno.faa ./protein_seqs/cflo_v7dot5
#proteinortho5 -step=2 -project=cflo_v_temno_1 -p=blastp+ ./protein_seqs/temno.faa ./protein_seqs/cflo_v7dot5

# cflo vs clonal
#################
proteinortho5 -step=1 -project=cflo_v_clonal_1 -p=blastp+ ./protein_seqs/clonal.faa ./protein_seqs/cflo_v7dot5
proteinortho5 -step=2 -project=cflo_v_clonal_1 -p=blastp+ ./protein_seqs/clonal.faa ./protein_seqs/cflo_v7dot5
