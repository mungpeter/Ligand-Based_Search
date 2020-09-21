#!/usr/bin/env python3

#########################################################################
##
##	Peter M.U. Ung @MSSM
##
##	v1.0	- 15.05.14  based on serial 0_fp_2_compare.py
##  v2.0  - 15.08.24  correct map_async error due to Python Update
##	v3.0	- 18.02.15  optimize the fingerprint calculation step
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

import re,glob
import gzip,bz2
import pandas as pd

from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import MACCSkeys
from rdkit.Chem.AtomPairs import Pairs
from rdkit.Chem import PandasTools as rdpd
from rdkit.Chem.Fingerprints import FingerprintMols

from tqdm import tqdm
from pathos import multiprocessing

Cols  = {'ecfp4':2, 'dl':3, 'apair':4, 'maccs':5, 'total':6}
Title = ['#id','ref','ECFP4','DL','APair','MACCS','Sum','SMILES']
dl_r  = 2
ec_r  = 2
ap_r  = 2
ms_r  = 1

##########################################################################
def main( reference, queue, rank_fp, mpi_np, output ):

  Ref  = RDkitRead(reference, removeHs=False)
  RfNm = Ref['ID'].tolist()
  Que  = RDkitRead(queue, removeHs=False)
  QuNm = Que['ID'].tolist()
  
  mpi = multiprocessing.Pool(processes = mpi_np)
  ## Generate fingerprints for all inputs
#  Rfp = [x for x in tqdm( mpi.imap( GenerateFPs, list(zip(Ref.mol.tolist(), RfNm))),
#                              total=len(list(RfNm)) ) ]
#  Qfp = [x for x in tqdm( mpi.imap( GenerateFPs, list(zip(Que.mol.tolist(), QuNm))), 
#                              total=len(list(QuNm)) ) ]
  Rfp = mpi.map( GenerateFPs, list(zip(Ref.mol.tolist(), RfNm)) )
  Qfp = mpi.map( GenerateFPs, list(zip(Que.mol.tolist(), QuNm)) )

  idx = len(Rfp)
  for ref in tqdm( Rfp ):

#    print(' # Working on reference molecule: '+ref[1])
#    print('  > {0}'.format(Chem.MolToSmiles(ref[0])))
#    pRef = FPFunctions()     # Calculate reference fp only once
    pLib = FPFunctions( ref_id = ref[1], ref_fp = ref[2:] )

    ## Load Lib after creating the CompareMolFingerprints object to avoid 
    ## duplication of memory use, otherwise for N processors to run MPI, 
    ## there will be MEM*N% memory to be used   
    #  Somehow when the list of mol object (Lib) is imported into any MPI
    #  process, the basic info in each ligand is lost. For any required 
    #  ligand info to be used in the MPI process, a list of the info will
    #  need to be extracted and supplied to the MPI process alongside with
    #  the list of mol objects thru ZIP function.

#    Data = [x for x in tqdm(mpi.imap(pLib, Qfp),total=len(Qfp))]
    Data = mpi.map(pLib, Qfp)

    ## Use Pandas to print result file
    Rank = sorted(Data, key=lambda x: x[Cols[rank_fp]], reverse=True)
    df = pd.DataFrame(Rank, columns=Title)
    if idx == 0:
      df.to_csv(output+'.txt', sep='\t', float_format='%.3f', index=False)
    else:
      df.to_csv(output+'.txt', sep='\t', float_format='%.3f',
                index=False, mode='a')
 
  mpi.close()
  mpi.join()


##########################################################################
def GenerateFPs( inp ):
  mol, name = inp
  return [ mol, name,
           AllChem.GetMorganFingerprintAsBitVect(mol, 2, 2048),
           FingerprintMols.FingerprintMol(mol),
           Pairs.GetAtomPairFingerprintAsBitVect(mol),
           MACCSkeys.GenMACCSKeys(mol) ]

##########################################################################
class FPFunctions(object):

  def __init__(self, ref_fp=[], ref_id=''):
    self.ref_fp = ref_fp
    self.ref_id = ref_id

  def __call__(self, Item):
    return self.CompareFPs(Item)

############
  def CompareFPs(self, Item):
    que, que_id = Item[:2]
    que_fp = Item[2:]
    que_sm = Chem.MolToSmiles(que)

    Tc = [ DataStructs.FingerprintSimilarity(self.ref_fp[x], que_fp[x])
              for x in range(4) ]

    Data  = [ que_id, self.ref_id, Tc[0], Tc[1], Tc[2], Tc[3],
              ( Tc[0]*dl_r + Tc[1]*ec_r + Tc[2]*ap_r + Tc[3]*ms_r )/(dl_r + ec_r + ap_r + ms_r),
              que_sm ]
    return Data


##########################################################################
## Read in SMILES or SDF input and add Hs to it
def RDkitRead( in_file, removeHs=True, add_Hs=False ):
  ## Read in SDF file; can choose to add hydrogens or not
  if re.search(r'.sdf', in_file):
    print(' # Reading SDF')
    df = rdpd.LoadSDF(  file_handle(in_file), removeHs=removeHs,
                        idName='ID', molColName='mol' )
    df['smiles'] = df.mol.apply(lambda m:Chem.MolToSmiles(Chem.RemoveHs(m)))
    if add_Hs:
      df['mol'] = df.mol.apply(Chem.AddHs)

  ## Read in SMILES file, check if there is a header "smiles"
  if re.search(r'.smi', in_file):
    print('# Reading SMI')
    with file_handle(in_file) as fi:
      if re.search('smi', str(fi.readline()), re.IGNORECASE):
        print(' # Smiles input has Header #\n')
        df = pd.read_csv(in_file, sep='\s+', comment='#').dropna()
        df.columns = ['smiles','ID']
      else:
        print(' # Smiles input has NO Header #\n')
        df = pd.read_csv(in_file, header=None, sep='\s+', comment='#').dropna()
        df.columns = ['smiles','ID']
    rdpd.AddMoleculeColumnToFrame(df, smilesCol='smiles', molCol='mol')
    df['smiles'] = df.mol.apply(Chem.MolToSmiles)

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
  main(sys.argv[1], sys.argv[2], sys.argv[3],
       int(sys.argv[4]), sys.argv[5])
