#!/usr/bin/env python3

#########################################################################
##
##	Peter M.U. Ung @MSSM
##
##	v1.0	- 15.05.14  based on serial 0_fp_2_compare.py
##      v2.0    - 15.08.24  correct map_async error due to Python Update
##      
##      v3.0    - 18.02.14  special case: SELF-comparing pairwise distance
##                          reduce calculation to O(N^2/2-N) cases
##
##	Calculate chemical similarity coefficient of 2 sets of input molecules
##	Takes gzip or bzip2 files of both SDF (.sdf) and SMILES (.smi) formats.
##	Files can contain multiple molecules.
##
##	Since pairwise matrix is diagonally mirroring, i.e. 2-1 = 1-2, can
##	reduce the needed non-redundant calculation to half, while self-pairs
##	are always 1, i.e. 1-1 = 2-2 = 999-999 == 1
##
##	This version performs asychronized parallelization to calculate fp.
##      Use maximum number of CPU automatically.
##
##	For smiles, title is turned off.
##      Return output file:
##	  [prefix].fp.txt   	all ligands and their coeff to template
##
#########################################################################

import sys
msg = '''\n  ## Usage: {0}
             [self-comparing MOL]  (.smi, .sdf, multi-mol file ok)
	     [Rank mol with Fingerprint: str] (ecfp4, dl, maccs, total)
             [Output Prefix]\n
             e.g.>   x.py a.smi total output\n
	     return: output.fp.txt (all scores)\n
'''.format(sys.argv[0])
if len(sys.argv) != 4: sys.exit(msg)

import gzip,bz2,re,glob
import itertools
import pandas as pd

from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import MACCSkeys
from rdkit.Chem import PandasTools as rdpd
from rdkit.Chem.Fingerprints import FingerprintMols

from tqdm import tqdm
from pathos import multiprocessing

Cols  = {'ecfp4':2, 'dl':3, 'maccs':4, 'total':5}
Title = ['#Ref','Que','ECFP4','DL','MACCS' ,'Sum'] #,'Ref_smi', 'Que_smi']
dl_r  = 2
ec_r  = 2
mc_r  = 1

##########################################################################
def main(template, rank_fp, output):

  Lib  = RDkitRead(template, keep_Hs=False)
  LbNm = Lib.ID.tolist()
  
  print(' # Working on Generating FPs for molecules: '+str(len(LbNm)))
  mpi = multiprocessing.Pool(processes = multiprocessing.cpu_count())
  FP  = [x for x in tqdm(mpi.imap_unordered(
                            GenerateFPs, list(zip(Lib.mol.tolist(), LbNm)) ),
                            total=len(LbNm) ) ]
#  FP   = mpi.map(GenerateFPs, list(zip(Lib, LbNm)))
  mpi.close()
  mpi.join()

#########################
  # Calculate all combinations of non-redundant pairs, also ignore self-pairs
  for i in range(2,3):
    Comb = list(itertools.combinations(FP,i))
 
  print(' # Working on pairwise FP Taninmoto calculation #')
#  Data = [CompareFPs(comb) for comb in Comb]
  mpi  = multiprocessing.Pool(processes = multiprocessing.cpu_count())
  Data = [x for x in tqdm(mpi.imap_unordered(
                                   CompareFPs, Comb),
                                   total=len(Comb) ) ]
#  Data = mpi.map(CompareFPs, Comb)
  mpi.close()
  mpi.join()

  ## Use Pandas to print result file
  Rank = sorted(Data, key=lambda x: x[Cols[rank_fp]], reverse=True)
  df = pd.DataFrame(Rank, columns=Title)

  df.to_csv(output+'.txt', sep='\t', float_format='%.3f',
              index=False )


##########################################################################
def GenerateFPs(inp):
  mol, name = inp
  mol_sm = Chem.MolToSmiles(mol)

  return [ name, mol_sm,
           AllChem.GetMorganFingerprintAsBitVect(mol, 2, 2048),
           FingerprintMols.FingerprintMol(mol),
           MACCSkeys.GenMACCSKeys(mol) ]


##########################################################################
def CompareFPs(comb):
  ref_FP, que_FP = comb
  ref_nm = ref_FP[0]
  que_nm = que_FP[0]
  ref_sm = ref_FP[1]
  que_sm = que_FP[1]

  Coeff  = [ DataStructs.FingerprintSimilarity(ref_FP[x], que_FP[x])
               for x in range(2,5) ]
  Info   = [ref_nm, que_nm, Coeff[0], Coeff[1], Coeff[2],
             CombineFPs(Coeff) ] #, ref_sm, que_sm]
  return Info


def CombineFPs(Coeff):
  return (Coeff[0]*dl_r + Coeff[1]*ec_r + Coeff[2]*mc_r)/(dl_r+ec_r+mc_r)

##########################################################################
## Read in SMILES or SDF input and add Hs to it
def RDkitRead( in_file, keep_Hs=True, add_Hs=False ):
  ## Read in SDF file; can choose to add hydrogens or not
  if re.search(r'.sdf', in_file):
    print(' # Reading SDF')
    df = rdpd.LoadSDF(  file_handle(in_file), removeHs=keep_Hs,
                        smilesName='smiles', molColName='mol' )
    if add_Hs:
      df['mol'] = df.mol.apply(Chem.AddHs)

  ## Read in SMILES file, check if there is a header "smiles"
  if re.search(r'.smi', in_file):
    print('# Reading SMI')
    with file_handle(in_file) as fi:
      if re.search('smi', str(fi.readline()), re.IGNORECASE):
        print('\n # Smiles input has Header #\n')
        df = pd.read_csv(in_file, sep='\s+', comment='#').dropna()
        df.columns = ['smiles','ID']
      else:
        print('\n # Smiles input has NO Header #\n')
        df = pd.read_csv(in_file, header=None, sep='\s+', comment='#').dropna()
        df.columns = ['smiles','ID']
    df['mol'] = df.smiles.apply(Chem.MolFromSmiles)

  print('## Number of MOL read from {}: {}\n'.format(in_file,len(df.smiles)))
  return df

#########################################################################
## Handle gzip and bzip2 file if the extension is right. otherwise, just open
## outuput: file handle
def file_handle(file_name):
  if re.search(r'.gz$', file_name):
    handle = gzip.open(file_name, 'rb')
  elif re.search(r'.bz2$', file_name):
    handle = bz2.BZ2File(file_name, 'rb')
  else:
    if re.search(r'.smi', file_name):
      handle = open(file_name, 'r')
    else:
      handle = file_name
  return handle


##########################################################################
if __name__ == '__main__':
  main(sys.argv[1], sys.argv[2], sys.argv[3])
