
import numpy as np
import pandas as pd

from rdkit import Chem
from rdkit import Geometry
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem import MolStandardize
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem.rdMolAlign import AlignMol

from rdkit.Chem import PandasTools as rdpd
import sys

acamide = pd.DataFrame(columns=['Name','Smiles'])
acamide = acamide.append(
  [ {'Smiles':'CC1(C(N2C(S1)C(C2=O)NC(=O)CC3=CC=CC=C3)C(=O)O)C','Name':'Penicilline_G'},
    {'Smiles':'CC1(C2CC3C(C(=O)C(=C(C3(C(=O)C2=C(C4=C1C=CC=C4O)O)O)O)C(=O)N)N(C)C)O','Name':'Tetracycline'},
    {'Smiles':'CC1(C(N2C(S1)C(C2=O)NC(=O)C(C3=CC=CC=C3)N)C(=O)O)C','Name':'Ampicilline'}], ignore_index=True)
#rdpd.AddMoleculeColumnToFrame(acamide,'Smiles','Molecule')
acamide.Molecule[0]

acamide['Molecule'] = acamide.Smiles.apply(CleanSmilesForMol)
acamide
def CleanSmilesForMol( smi ):
  std_smi = MolStandardize.standardize_smiles(Chem.CanonSmiles(smi))
  tau_smi = MolStandardize.canonicalize_tautomer_smiles(std_smi)
  tau_mol = Chem.MolFromSmiles(tau_smi)
  mol_h   = Chem.AddHs(tau_mol)
  return mol_h


confgen = ETKDG_ConfGenerator(numConfs=0, run_uff=1)
Data = confgen(acamide.loc[0])
Data
Data[0].GetNumConformers()
mol = Data[0]
confs = mol.GetConformers()
AllChem.UFFGetMoleculeForceField(mol, confid=1)

x = [AllChem.UFFGetMoleculeForceField(mol, confid=conf) for conf in mol.GetConformers()]


new = prune_conformers(Data[0],50)

conversion = ConvertConformerResult()
xxx = conversion(Data)
xxx[3].GetProp('_Name')
xxx
for data in [Data]:
  m, name, confs = data
  outStream = Chem.SDWriter('test.sdf')
  for _id in ids:
    outStream.write(m, confId=_id)
    outStream.flush()
  outStream.close()
  sys.stderr.write('\n')


class ETKDG_ConfGenerator(object):
  def __init__( self, numConfs=1, run_uff=False ):
    self.numConfs  = numConfs
    self.run_uff   = run_uff

    self.conf_parm = AllChem.ETKDGv2()
    self.conf_parm.pruneRmsThresh = 0.5
    self.conf_parm.randomSeed     = -1

  def __call__( self, mol ):
    return self.conf_generator( mol )

#######################################
  def adaptive_conf( self, mol ):
    '''
        default conformer number based on no. rotatable bond
          rb >= 13       max_conformers = 300
          10 < rb <= 13  max_conformers = 250
          8  < rb <= 10  max_conformers = 200
          5  < rb <= 8   max_conformers = 150
          3  < rb <= 5   max_conformers = 100
          rb <= 3        max_conformers = 50
        override these settings if numConfs != 0
    '''
    rb_num = Chem.rdMolDescriptors.CalcNumRotatableBonds(mol)
    if  rb_num <= 3:
      max_conf = 50
    elif rb_num > 3  and rb_num <= 5:
      max_conf = 100
    elif rb_num > 5  and rb_num <= 8:
      max_conf = 150
    elif rb_num > 8  and rb_num <= 10:
      max_conf = 200
    elif rb_num > 10 and rb_num <= 13:
      max_conf = 250
    else:
      max_conf = 300

    # enforce input conformer number
    if self.numConfs != 0:
      max_conf = self.numConfs
    return max_conf

################################
  def conf_generator( self, inp ):
    mol  = inp.Molecule
    name = inp.Name

    ## get rotatable bond-dependent adaptive conformation number
    numConfs = self.adaptive_conf(mol)

    ## Generate 3D conformers, map atom 3D vectors in 'ids' to 'mol'
    mol = Chem.AddHs(mol)
    ids = AllChem.EmbedMultipleConfs(mol, numConfs, self.conf_parm)

    ## align all conformers to 1st frame
    rmslist = []
    AllChem.AlignMolConformers(mol, RMSlist=rmslist)

    ## Minimize conformers with UFF, 2x slower than without
    if self.run_uff:
      for _id in ids: 
        AllChem.UFFOptimizeMolecule(mol, confId = _id)

    return [mol, name, ids]

##########################################
def prune_conformers( mol, num_max_conf, rmsd_threshold=0.5):

  if rmsd_threshold < 0 or mol.GetNumConformers() <= 1:
    return mol

  rmsd = get_conformer_rmsd(mol)

  keep = []  # always keep lowest-rmsd conformer
  discard = []
  ## sort by increasing RMSD
  for i in np.argsort(rmsd):

    # always keep lowest-energy conformer
    if len(keep) == 0:
      keep.append(i)
      continue

    # discard conformers after max_conformers is reached
#    if len(keep) >= self.max_conformers:
    if len(keep) >= num_max_conf:
      discard.append(i)
      continue

    # get RMSD to selected conformers
    this_rmsd = rmsd[i][np.asarray(keep, dtype=int)]

    # discard conformers within the RMSD threshold
    if np.all(this_rmsd >= rmsd_threshold):
      keep.append(i)
    else:
      discard.append(i)

  # create a new molecule to hold the chosen conformers
  # this ensures proper conformer IDs and energy-based ordering
  new = Chem.Mol(mol)
  new.RemoveAllConformers()
  conf_ids = [conf.GetId() for conf in mol.GetConformers()]
  print(len(keep))
  for i in keep:
    conf = mol.GetConformer(conf_ids[i])
    new.AddConformer(conf, assignId=True)
  return new

###################################################################
def get_conformer_rmsd(mol):
  rmsd = np.zeros((mol.GetNumConformers(), mol.GetNumConformers()),
                        dtype=float)
  for i, ref_conf in enumerate(mol.GetConformers()):
    for j, fit_conf in enumerate(mol.GetConformers()):
      if i >= j:
        continue
#                rmsd[i, j] = AllChem.GetBestRMS(mol, mol, ref_conf.GetId(),
#                                                fit_conf.GetId())
#      GetBestRMS may go 'permutation explosion' for large molecules
      rmsd[i, j] = AlignMol(mol, mol,
                            ref_conf.GetId(), fit_conf.GetId(),
                            maxIters=200)
#                print( 'xxx {0} {1} {2}'.format(i, j, rmsd[i, j]))
      rmsd[j, i] = rmsd[i, j]
  return rmsd



###################################################################
class ConvertConformerResult(object):
  def __init__( self, cg=None, log=None ):
    self.cg  = cg
    self.log = log
  
  def __call__( self, inp ):
    return self._convert_conformer( inp )

  def _convert_conformer( self, inp ):
    m, name, ids = inp
    Mols = []

    for id in range(0, m.GetNumConformers()):
      mol = Chem.MolFromMolBlock(Chem.MolToMolBlock(m,confId=id),removeHs=False)
      mol.SetProp('_Name', name+'_'+str(id+1))
      mol.SetProp('Group', name)
      mol.SetProp('Conformation', str(id+1))
      Mols.append(mol)

    ## write how many conformer generated for each molecule
#    self.log.write('{0}\t{1}\n'.format(name, len(Mols)))
    return Mols
