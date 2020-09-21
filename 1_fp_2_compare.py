#!/usr/bin/env python3

import sys
msg = '''\n  {0}
      -ref  < >  [ Reference MOL ]   (.smi, .sdf ok; gzip,bzip2 ok)
      -que  < >  [ Query MOL ]  
      -rank < >  [ Rank mol with Fingerprint: str ] (ecfp4, dl, apair, maccs, wavg)
      -pref < >  [ Output Prefix ]\n
Optional:
      -rid  < >  [ SDF tag of Reference MOL (Def: ID) ]
      -qid  < >  [ SDF tag of Query MOL (Def: ID) ]
      -coff < >  [ Cutoff to Similarity Coefficient (-rank) to save (Def: 0.25) ]
      -dice      [ Using *Dice* instead of *Tanimoto* for Similarity (Def: Tanimoto) ]
      -ratio <+> [ Ratio of ECFP4:Daylight:AtomPair:MACCS Weighted Average (Def: 2 2 2 1) ]\n
e.g.> x.py -ref a.smi -que b.sdf.bz2 -rank ecfp4 -pref output -ratio 2 2 1 0 -dice\n
      return: output.txt (all scores); output.sdf.gz (new tag to -ref)\n'''.format(sys.argv[0])
if len(sys.argv) < 5: sys.exit(msg)

import re
import gzip,bz2
import pandas as pd

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import MACCSkeys
from rdkit.Chem.AtomPairs import Pairs
from rdkit.Chem import PandasTools as rdpd
from rdkit.DataStructs import DiceSimilarity
from rdkit.DataStructs import TanimotoSimilarity
from rdkit.DataStructs import FingerprintSimilarity
from rdkit.Chem.Fingerprints import FingerprintMols

from tqdm import tqdm
from argparse import ArgumentParser

Title = ['Query','Ref','ECFP4','DL','APair','MACCS','W_Avg','Query_smi','Ref_smi']
fp_list = ['ecfp4','dl','apair','maccs']

##########################################################################
def main( ):

  args = UserInput()
  if args.cutoff:
    cutoff = float(args.cutoff)
  else:
    cutoff = 0.25
  if args.rid:
    r_id = args.rid
  else:
    r_id = 'ID'
  if args.qid:
    q_id = args.qid
  else:
    q_id = 'ID'
  ## fp_ratio is ordered: [ ecfp4, daylight, atom-pairs, maccs ] 
  if args.ratio:
    fp_ratio = [float(args.ratio[0]),float(args.ratio[1]),float(args.ratio[2]),float(args.ratio[3])]
  else:
    fp_ratio = [2., 2., 2., 1.]
  if args.dice:
    coeff = DiceSimilarity
  else:
    coeff = TanimotoSimilarity

##########################

  ## Import data and calculate FPs
  Ref = RDkitRead(args.reference, removeHs=False)
  Ref['ecfp4'] = Ref.mol.apply(lambda m:AllChem.GetMorganFingerprintAsBitVect(m, 2, 2048))
  Ref['apair'] = Ref.mol.apply(lambda m:Pairs.GetAtomPairFingerprintAsBitVect(m))
  Ref['dl']    = Ref.mol.apply(lambda m:FingerprintMols.FingerprintMol(m))
  Ref['maccs'] = Ref.mol.apply(lambda m:MACCSkeys.GenMACCSKeys(m))

  Que = RDkitRead(args.query, removeHs=False)
  Que['ecfp4'] = Que.mol.apply(lambda m:AllChem.GetMorganFingerprintAsBitVect(m, 2, 2048))
  Que['apair'] = Que.mol.apply(lambda m:Pairs.GetAtomPairFingerprintAsBitVect(m))
  Que['dl']    = Que.mol.apply(lambda m:FingerprintMols.FingerprintMol(m))
  Que['maccs'] = Que.mol.apply(lambda m:MACCSkeys.GenMACCSKeys(m))

  ## Calculate Similarity Coefficient
  df_list, df_tops = [], []
  for idx, row in tqdm(Ref.iterrows(), total=len(Ref)):
    FPCalc = FPDistCalculator(ref_ec=row.ecfp4, ref_dl=row.dl, ref_ap=row.apair, 
                              ref_ms=row.maccs, ref_sm=row.smiles, ref_id=row[r_id],
                              que_id=q_id, fp_r=fp_ratio, coeff=coeff )
    Ranked = FPCalc(Que).sort_values(by=[args.rank_fp], ascending=False)
    df_list.append(Ranked)
    df_tops.append(Ranked[0:1])


  ## summarize and output similarity between Ref and Query mol
  all_df = pd.concat(df_list).reset_index(drop=True)
  all_df = all_df[all_df[args.rank_fp] >= cutoff]
  all_df.columns = Title
  all_df.to_csv(args.output+'.txt', sep='\t', float_format='%.3f', index=False)

  ## put most similar Ref smiles into Query mol
  top_df = pd.concat(df_tops).reset_index(drop=True)
  Ref['check_ref_id']= top_df['ref_id']
  Ref['top_que_id']  = top_df['que_id']
  Ref['top_que_ec4'] = top_df['ecfp4']
  Ref['top_que_dl']  = top_df['dl']
  Ref['top_que_ap']  = top_df['apair']
  Ref['top_que_ms']  = top_df['maccs']
  Ref['top_que_wavg']= top_df['wavg']
  Ref['top_que_smiles'] = top_df['que_smi']
  Ref.drop(columns=['ecfp4','dl','apair','maccs'], inplace=True)
  rdpd.WriteSDF(Ref, args.output+'.sdf.gz', molColName='mol', properties=list(Ref.columns))


##########################################################################
## Calculate Similarity Coefficient for different FPs and their weighted average
class FPDistCalculator(object):
  def __init__( self, ref_ec='', ref_dl='', ref_ap='', ref_ms='', ref_sm='',
                      ref_id='', que_id='', fp_r=[],   coeff=''  ):
    self.ref_ec = ref_ec  # reference fp: ecfp4
    self.ref_dl = ref_dl  # reference fp: daylight
    self.ref_ap = ref_ap  # reference fp: atom pairs
    self.ref_ms = ref_ms  # reference fp: maccs
    self.ref_sm = ref_sm  # reference smiles
    self.ref_id = ref_id  # reference ID column name
    self.que_id = que_id  # query ID column name
    self.fp_r   = fp_r    # ratio of fps
    self.coeff  = coeff   # similarity coefficient function

  def __call__( self, Query ):
    return self.CompareFPs(Query)

############
  def _weighted_avg( self, ecfp4, dl, apair, maccs ):
    ec_r, dl_r, ap_r, ms_r = self.fp_r
    return ( ecfp4*ec_r + dl*dl_r + apair*ap_r + maccs*ms_r )/(ec_r + dl_r + ap_r + ms_r)

  ## Default similarity algorithm is Tanimoto. To control simi algorithm, add
  ## "metric=DataStructs.DiceSimilarity" into FingerprintSimilarity()
  def CompareFPs( self, Query ):
    df = pd.DataFrame(columns=['_x'], index=list(range(len(Query))))
    df['que_id'] = Query[self.que_id]
    df['ref_id'] = self.ref_id
    df['ecfp4']  = Query.ecfp4.apply(lambda m: FingerprintSimilarity(self.ref_ec, m, metric=self.coeff))
    df['dl']     = Query.dl.apply(   lambda m: FingerprintSimilarity(self.ref_dl, m, metric=self.coeff))
    df['apair']  = Query.apair.apply(lambda m: FingerprintSimilarity(self.ref_ap, m, metric=self.coeff))
    df['maccs']  = Query.maccs.apply(lambda m: FingerprintSimilarity(self.ref_ms, m, metric=self.coeff))
    df['wavg']   = self._weighted_avg(df.ecfp4, df.dl, df.apair, df.maccs)
    df['que_smi']= Query.smiles
    df['ref_smi']= self.ref_sm
    df.drop(columns=['_x'], inplace=True)
    return df


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
        print('# Smiles input has Header #\n')
        df = pd.read_csv(in_file, sep='\s+', comment='#').dropna()
        df.columns = ['smiles','ID']
      else:
        print('# Smiles input has NO Header #\n')
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
def UserInput():
  p = ArgumentParser(description='Command Line Arguments')

  p.add_argument('-ref', dest='reference', required=True,
                  help='Reference MOL (sdf|smi) * gzip,bzip2 accepted')
  p.add_argument('-que', dest='query', required=True,
                  help='Query MOL (sdf|smi) * gzip,bzip2 accepted')
  p.add_argument('-rank', dest='rank_fp', required=True,
                  help='Rank mol with Fingerprint ("ecfp4", "dl", "apair", "maccs", "wavg")')
  p.add_argument('-pref', dest='output', required=True,
                  help='Output prefix')

## Optional
  p.add_argument('-rid', dest='rid', required=False,
                  help='Optional: SDF tag of Reference MOL identifier (Def: ID)')
  p.add_argument('-qid', dest='qid', required=False,
                  help='Optional: SDF tag of Query MOL identifier (Def: ID)')
  p.add_argument('-coff', dest='cutoff', required=False,
                  help='Optional: Cutoff to Similarity Coefficient (use -rank) to save (Def: 0.25)')
  p.add_argument('-dice', dest='dice', required=False, action='store_true',
                  help='Optional: Using *Dice* instead of *Tanimoto* for Similarity (Def: Tanimoto)')
  p.add_argument('-coeff', dest='coeff', required=False,
                  help='Optional: Cutoff to Tanimoto Coefficient to save (use -rank)')
  p.add_argument('-ratio', dest='ratio', required=False, nargs='+',
                  help='Optional: Ratio of ecfp4:dl:apair:maccs Tc average (Def: 2 2 2 1)')

  return p.parse_args()


##########################################################################
if __name__ == '__main__':
  main()


#########################################################################
##
##	Peter M.U. Ung @ gRED
##
##  *OLD* version using MPI: 1_fp_2_compare.mpi.py

##  v1.0  - 20.09.19  modernize with rdkit PandasTools

##	Calculate chemical similarity coefficient of 2 sets of input molecules
##	Takes gzip files of both SDF (.sdf) and SMILES (.smi) formats.
##	Files can contain multiple molecules.

##	This version uses single-cpu + Pandas to calculate fp similarity.
##  For comparison, a 767x850 pairwise set:
##  2-cpu mpi-only version: ~ 170s (4 FPs)
##  8-cpu mpi-only version: ~ 145s (4 FPs)
##  1-cpu Pandas version:   ~  40s (4 FPs)

##  MPI-only is simply inefficient for this type of many-small-but-fast job.
##  Overhead to create the many processes is too much.

##  Dask is mpi-cpu + Pandas but also not good for this job. partitions=2 is ~6x
##  slower than single-cpu; partitions=8 is ~10x slower! Probably has to do
##  with number of small and fast jobs to perform where mpi overhead is too
##  costly for this type of job.
##  For dask, need to do: x.apply(x, meta=pd.Series(dtype='int',column='ecfp4'))
##
##  Return output files:
##	  [prefix].txt       pairwise query-reference similarity coeff
##    [prefix].sdf.gz    1 closest query molecule for each reference molecule
##
#########################################################################