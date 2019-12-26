#!/bin/csh

## Pairwise FP comparison can be reduced to O((N^2)/2-N) instead of O(N^2)
## because the matrix is symmetrical
## 'Total' is the sum of weighted Daylight (2), ECFP_4 (2), MACCS (1) FPs
## in 2:2:1 ratio. 

../../1_fp_2_compare.symmetric-mpi.py          \
  ../recpt.zfg19.sch_top500.pick.sdf.bz2      \
  total                                       \
  recpt.zfg19.pairwise.fp_all


## recpt.zfg19.pairwise.fp_all.txt
#  pairwise mol-mol comparison with Daylight, ECFP_4, MACCS, Total FPs

# 19.12.19
