#!/usr/bin/env python3

##########################################################################
#
#   Peter M.U. Ung @ MSSM
#
#   v1.0    17.12.12    
#
#   Calculate molecular properties via RDKit
#
###########################################################################

import os,sys

msg = '''\n\t{0}\n\t\t[library of ligands smi|sdf]
\t\t[Output .csv/.xlsx prefix]
\n\te.g.>\t x.py mol.smi mol_result\n'''.format(sys.argv[0])
if len(sys.argv) != 3: sys.exit(msg)

import gzip,bz2,re
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors as rd


##########################################################################

def main(in_file, output):

  Cmpds  = {}
  InMols = rdkit_open([in_file])
  print('\n # Number of input molecule: {0}'.format(len(InMols)))
  for mol in InMols:
    m = {}

    name = mol.GetProp('_Name').split()[0]
    
    m['Name'] = name
    m['Formula'] = rd.CalcMolFormula(mol)
    m['SMILES'] = Chem.MolToSmiles(mol)

    m['MW']   = rd._CalcMolWt(mol)               # Molecular Weight
    m['logP'] = rd.CalcCrippenDescriptors(mol)[0]  # Partition coefficient
    m['HDon'] = rd.CalcNumLipinskiHBD(mol)      # Lipinski Hbond donor
    m['HAcc'] = rd.CalcNumLipinskiHBA(mol)      # Lipinski Hbond acceptor
    m['TPSA'] = rd.CalcTPSA(mol)                # Topological polar surface area

    m['Rotat'] = rd.CalcNumRotatableBonds(mol, strict=True) # Rotatable bond
    m['MolRef'] = rd.CalcCrippenDescriptors(mol)[1]         # Molar refractivity
    m['AliRing'] = rd.CalcNumAliphaticRings(mol)        # Aliphatic ring number
    m['AroRing'] = rd.CalcNumAromaticRings(mol)         # Aromatic ring number
#    m['Stereo'] = rd.CalcNumAtomStereoCenters(mol)      # Stereo center number
#    m['UnspStereo'] = rd.CalcNumUnspecifiedAtomStereoCenters(mol)  # unspecified stereo

    m['SMILES'] = Chem.MolToSmiles(mol, 
                    isomericSmiles=True, allHsExplicit=False)
    Cmpds[name] = m

  ####################################

  df = pd.DataFrame.from_dict(Cmpds, orient='index')
  df.index.name = 'Name'

  # Columns of data to print out
  Columns = [ 'Formula',
              'MW',    'logP',   'HDon',    'HAcc',    'TPSA',
              'Rotat', 'MolRef', 'AliRing', 'AroRing', 
              #'Stereo', 'UnspStereo', 
              'SMILES', ]
  reorder = df[Columns]

  # Output to CSV
  reorder.to_csv( output+'.csv', sep=',', na_rep='NA', encoding='utf-8',
                  float_format='%.5f', header=True )

  # Output to Excel
  reorder.to_excel( output+'.xlsx', header=True, na_rep='NA' )


##########################################################################
def file_handle(file_name):
  if re.search(r'.gz$', file_name):
    handle = gzip.open(file_name, 'r')
  elif re.search(r'.bz2$', file_name):
    handle = bz2.BZ2File(file_name, 'r')
  else:
    handle = open(file_name)

#  print "## Opening "+file_name
  return handle

#######################################################################
## new version of rdkit distinguish the input source of the file, treating
## regular utf-8 file as str input and bytes file (zipped) as object input.
## Forward_supplier only takes bytes files and Regular_supplier takes regulars.
## To get around this, use file('x.sdf') to make utf-8 file as an object.

def rdkit_open(File_Tuple):

  List = []

  for f in (File_Tuple):
    handle = file_handle(f)

    if re.search(r'.sdf', f):
      Mol = [x for x in Chem.ForwardSDMolSupplier(handle, removeHs=True)
             if x is not None]

    if re.search(r'.smi', f):
      with handle as fi:
        first_line = fi.readline()

      if re.search(r'smiles', first_line, re.IGNORECASE):
        Mol = [x for x in Chem.SmilesMolSupplier(f, titleLine=True,
                 delimiter=' |\t|,') if x is not None]
      else:
        Mol = [x for x in Chem.SmilesMolSupplier(f, titleLine=False,
                 delimiter=' |\t|,') if x is not None]

    print( "# Found mol in {0}: {1}".format(f, len(Mol)))
    for mol in Mol: List.append(mol)

  return List

#######################################################################
##########################################################################
if __name__ == '__main__':
  main( sys.argv[1], sys.argv[2] )
