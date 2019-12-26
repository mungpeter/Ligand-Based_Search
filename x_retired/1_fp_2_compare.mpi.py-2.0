#!/usr/bin/python

#########################################################################
##
##	Peter M.U. Ung @MSSM
##
##	v1.0	- 15.05.14  based on serial 0_fp_2_compare.py
##      v2.0    - 15.08.24  correct map_async error due to Python Update
##
##	Calculate chemical similarity coefficient of 2 sets of input molecules
##	Takes gzip or bzip2 files of both SDF (.sdf) and SMILES (.smi) formats.
##	Files can contain multiple molecules.
##	This version performs asychronized parallelization to calculate fp.
##
##	For smiles, title is turned off.
##      Return output file:
##	  [prefix].fp.txt   	all ligands and their coeff to template
##
#########################################################################

import sys
msg = '''\n  ## Usage: {0}
             [Template MOL] [Compare MOL]  (.smi, .sdf, multi-mol file ok)
	     [Rank mol with Fingerprint: str] (ecfp4, dl, maccs, total)
	     [No. Processor for Parallel: int]
             [Output Prefix]\n
             e.g.>   x.py a.smi b.sdf.bz2 ecfp4 4 output\n
	     return: output.fp.txt (all scores)\n
'''.format(sys.argv[0])
if len(sys.argv) != 6: sys.exit(msg)

from rdkit import Chem
from rdkit import DataStructs
from rdkit_open import *

from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem import MACCSkeys
from rdkit.Chem import AllChem
from CommonUtility import *
import multiprocessing
import pandas as pd

Cols  = {'ecfp4':2, 'dl':3, 'maccs':4, 'total':5}
Title = ['#id','ref','ECFP4','DL','MACCS','Sum','SMILES']

##########################################################################
def main(template, library, rank_fp, mpi_np, output):
#  mpi = multiprocessing.Pool(processes = mpi_np)

  for idx, ref in enumerate(rdkit_open([template])):
    mpi = multiprocessing.Pool(processes = mpi_np)

    print ' # Working on reference molecule: '+ref.GetProp('_Name')
    print '  > {0}'.format(Chem.MolToSmiles(ref))
    pRef = CompareMolFingerprints()     # Calculate reference fp only once
    pLib = CompareMolFingerprints(
             ref_fp = pRef.generate_fingerprints(ref),
	     ref_id = ref.GetProp('_Name'))
    print '   Calculating fingerprint Tanimoto Coefficient'

    ## Load Lib after creating the CompareMolFingerprints object to avoid 
    ## duplication of memory use, otherwise for N processors to run MPI, 
    ## there will be MEM*N% memory to be used   
    #  Somehow when the list of mol object (Lib) is imported into any MPI
    #  process, the basic info in each ligand is lost. For any required 
    #  ligand info to be used in the MPI process, a list of the info will
    #  need to be extracted and supplied to the MPI process alongside with
    #  the list of mol objects thru ZIP function.
    Lib  = rdkit_open([library])
    LbNm = [m.GetProp('_Name') for m in Lib]
    Data = mpi.map(pLib, zip(Lib, LbNm))
    mpi.close()
    mpi.join()

    ## Use Pandas to print result file
    Rank = sorted(Data, key=lambda x: x[Cols[rank_fp]], reverse=True)
    df = pd.DataFrame(Rank, columns=Title)
    if idx == 0:
      df.to_csv(output+'.fp.txt', sep='\t', float_format='%.3f', index=False)
    else:
      df.to_csv(output+'.fp.txt', sep='\t', float_format='%.3f',
                index=False, mode='a')
 

##########################################################################
class CompareMolFingerprints(object):

  def __init__(self, ref_fp=[], ref_id=''):
    self.ref_fp = ref_fp
    self.ref_id = ref_id

  def __call__(self, Item):
    return self.compare_fingerprints(Item)

##########################################################################
  def compare_fingerprints(self, Item):
    mol, mol_id = Item

    mol_sm = Chem.MolToSmiles(mol)

    FP     = self.generate_fingerprints(mol)
    ref_FP = self.ref_fp
    Coeff  = [ DataStructs.FingerprintSimilarity(ref_FP[x], FP[x])
                 for x in range(3) ]

    Data   = [mol_id, self.ref_id, Coeff[0], Coeff[1], Coeff[2], 
              sum(Coeff), mol_sm]
    return Data


  def generate_fingerprints(self, mol):
    return [ AllChem.GetMorganFingerprintAsBitVect(mol, 2, 2048),
             FingerprintMols.FingerprintMol(mol),
             MACCSkeys.GenMACCSKeys(mol) ]


##########################################################################

if __name__ == '__main__':
  main(sys.argv[1], sys.argv[2], sys.argv[3],
       int(sys.argv[4]), sys.argv[5])
