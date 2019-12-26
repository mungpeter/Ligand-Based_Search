#!/bin/csh

## Asymmetrical FP comparison is inefficient with order O(N^2)
## Use MPI to accelerate the calculation
##
## 'Total' is the sum of weighted Daylight (2), ECFP_4 (2), MACCS (1) FPs
## in 2:2:1 ratio. 


../../1_fp_2_compare.mpi.py                \
  ../query_mol.smi                         \
  ../recpt.zfg19.sch_top500.pick.sdf.bz2   \
  total                                    \
  4                                        \
  recpt.zfg19.query_mol.fp_all

## recpt.zfg19.query_mol.fp_all.txt
#  result of all FP comparison, numbers are Tanimoto coefficient

#  19.12.19
