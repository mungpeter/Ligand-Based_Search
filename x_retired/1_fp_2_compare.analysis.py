#!/usr/bin/env python3

##########################################################################
#
#   Peter M.U. Ung @ MSSM
#
#   v1.0    -   15.05.18
#
#   Purpose:    Take in a set of fingerprint datafiles generated by RDKIT - 
#               1_fp_2_compare.py, parse the data and print out the ligands
#               with Tanimoto coefficient higher than the threshold value.
#
#               Redundent SMILES are not removed since RDKit fp comparison
#               does not distinguish chiral centers so IDs can be different.
#   
#               1_fp_2_compare.py
#               2_fp_2_compare.run.csh
#
##########################################################################
import sys,glob,csv
from CommonUtility import *

msg = '''
    Usage: {0}
             [Input File(s)] 
	     [Fingerprint for Sorting: str] ('ecfp4', 'dl', 'maccs', 'total')
	     [Threshold to Select:   float]
	     [Output Smiles Filename]\n
    e.g.: x.py '*.txt.bz2' ecfp4 0.5 output.smi
      '''.format(sys.argv[0])
if len(sys.argv) != 5: sys.exit(msg)

inp = sys.argv[1]
fp  = sys.argv[2]
thr = float(sys.argv[3])
out = sys.argv[4]

Cols  = {'ecfp4':2, 'dl':3, 'maccs':4, 'total':5}
Title = ['#id','ref', 'ecfp4','dl','maccs','total','smiles']

##########################################################################
fo   = open(out, 'w')
fc   = csv.writer(open(out+'.csv','wb'), delimiter='\t')
In   = glob.glob(inp)
fc.writerow(Title)

count = 0
for inp in In:
  Data = []
  print('# Reading {0} ...'.format(inp))
  with file_handle(inp) as fi:
    for line in fi:
      Itm = line.split()
      try: float(Itm[2])
      except ValueError: continue
      Data.append(line.split())
  Pass = [Itm for Itm in Data if float(Itm[Cols[fp]]) > thr]
  Pass.sort(key=lambda x: float(x[Cols[fp]]), reverse=True)

  print('  Found items: '+str(len(Data)))
  print('  Passed:      '+str(len(Pass)))
  count += len(Pass)
  for Itm in Pass:
    fo.write('{0} {1}\n'.format(Itm[-1], Itm[0]))
    fc.writerow(Itm)
fo.close()

print('Total: '+str(count))
