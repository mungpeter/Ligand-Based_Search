# Ligand-Based_Search
Ligand-based searching for molecules using FingerPrint similarity

- Molecules can be described by the chemical Fingerprints and these fingerprints can be compared with each other to calculate Tanimoto coefficient (standard) or other newer coefficients (e.g. Dice). Three commonly used Fingerprints availble in [RDKit](https://www.rdkit.org/UGM/2012/Landrum_RDKit_UGM.Fingerprints.Final.pptx.pdf): DayLight (dl), ECFP_4 (ecfp4), and 166-bit MACCS key (maccs). However, each fingerprint has it strengths and weaknesses:
. ECFP in principle has unlimited byte space to describe any molecule, but this unlimited space can be to vast and the Tanimoto coefficient will be lower than other FPs
. MACCS is limited to the 166-bit public keys to describe the molecule
- To get a good balance between the different FPs, I combine the Tanimoto coefficients of the FPs DL:ECFP_4:MACCS in a _**2:2:1**_ ratio to get a **Total** Tanimoto coefficient. In practice, it does better than any FP alone. An example would be to pick out analogs of **Sorafenib** from a set of kinase inhibitors. Any single FP will have a harder time to pick out **regorafenib** and **LS1-15**, which differ by only 1 additional Fluorine atom. **Total** coefficient on the other hand can pick them out.
 
> Note: RDKit uses Morgan FP, same principle as ECFP.
#######################################################################################
```
|---- /Examples
          |---- /1_query_compare_fp    # compare 2 list of molecules using FP
          |---- /2_self_compare_fp     # self compare of a list of molecules
          |---- /3_lig_property        # generate molecular properties
```
#####################################
```
> 1_fp_2_compare.mpi.py
      [Template MOL] [Compare MOL]  (.smi, .sdf, multi-mol file ok)
	    [Rank mol with Fingerprint: str] (ecfp4, dl, maccs, total)
	    [No. Processor for Parallel: int]
      [Output Prefix]

   e.g.> *.py ../query_mol.smi   ../recpt.zfg19.sch_top500.pick.sdf.bz2
              total  4  recpt.zfg19.query_mol.fp_all

return >  recpt.zfg19.query_mol.fp_all.txt
```
- Compare one list of molecules to another list of molecules. The number of calculations is proportional to the number of item, _O_(MxN). If the 2 lists are the same, use the script **1_fp_2_compare.symmetric-mpi.py** below to reduce the number of calculations.

```
> 1_fp_2_compare.symmetric-mpi.py
      [ Self-comparing Mol        ]  # smi, sdf
      [ Rank Mol with Fingerprint ]  # ecfp_4, dl, maccs, total
      [ Output prefix ]

   e.g.> *.py  ../recpt.zfg19.sch_top500.pick.sdf.bz2  total                                       \
               recpt.zfg19.pairwise.fp_all

result > recpt.zfg19.pairwise.fp_all.txt
```
- Compare a list of molecules to itself. Since it is a symmetric pairwise comparison, it will be a (NxN) matrix and allows a trick to reduce the number of calculations to the order of _O_((N^2)/2 - N) instead of _O_(N^2).

#####################################
```
> 3_lig_property.py
      [ Library of ligand:  smi|sdf ]
      [ Outprefix ]

   e.g.> *.py  ../recpt.zfg19.sch_top500.pick.sdf.bz2
               recpt.zfg19.property

result > recpt.zfg19.property.csv   recpt.zfg19.property.xlsx
```
- Generate a set of chemical properties relevant to medicinal chemistry. Name, Formula, MW, logP, HDon, HAcc, TPSA, Rotat, MolRef, AliRing, AroRing, SMILES

#######################################################################################
# Required software / packages
```
LIGSHFT
```

```
csh/tcsh          # shell
python            # 3.6.8+
  numpy           # 1.16.2+
  pandas          # 0.24.2+
  matplotlib      # 3.0.3+
  rdkit           # 2019.09.02
  pathos          # 0.2.3+
  tqdm            # 4.31.1+
  argparse        # 1.1
  html            #
  gzip            #
  bz2             #
```
