# Ligand-Based_Search
Ligand-based searching for molecules using FingerPrint similarity

- Molecules can be described by the chemical Fingerprints and these fingerprints can be compared with each other to calculate Tanimoto coefficient (standard) or other newer coefficients (e.g. Dice). Three commonly used Fingerprints availble in [RDKit](https://www.rdkit.org/UGM/2012/Landrum_RDKit_UGM.Fingerprints.Final.pptx.pdf): DayLight (dl), ECFP_4 (ecfp4), and 166-bit MACCS key (maccs). However, each fingerprint has it strengths and weaknesses:
1. ECFP in principle has unlimited byte space to describe any molecule, but this unlimited space can be to vast and the Tanimoto coefficient will be lower than other FPs
2. MACCS is limited to the 166-bit public keys to describe the molecule

- To get a good balance between the different FPs, I combine the Tanimoto coefficients of the FPs DL:ECFP_4:MACCS in a _**2:2:1**_ ratio to get a **Total** Tanimoto coefficient. In practice, it does better than any FP alone. An example would be to pick out analogs of **Sorafenib** from a set of kinase inhibitors. Any single FP will have a harder time to pick out **regorafenib** and **LS1-15**, which differ by only 1 additional Fluorine atom. **Total** coefficient on the other hand can pick them out.

> Note: RDKit uses Morgan FP, same principle as ECFP.<

#######################################################################################
- Folder structure
```
|---- /Examples
          |---- /1_query_compare_fp    # compare 2 list of molecules using FP
          |---- /2_self_compare_fp     # self compare of a list of molecules
          |---- /3_lig_property        # generate molecular properties
```
#####################################
```
> 1_fp_2_compare.py
      -ref  < >  [ Reference MOL ]   (.smi, .sdf ok; gzip,bzip2 ok)
      -que  < >  [ Query MOL ]  
      -rank < >  [ Rank mol with Fingerprint: str ] (ecfp4, dl, apair, maccs, wavg)
      -pref < >  [ Output Prefix ]

Optional:
      -rid  < >  [ SDF tag of Reference MOL (Def: ID) ]
      -qid  < >  [ SDF tag of Query MOL (Def: ID) ]
      -coff < >  [ Cutoff to Similarity Coefficient (-rank) to save (Def: 0.25) ]
      -dice      [ Using *Dice* instead of *Tanimoto* for Similarity (Def: Tanimoto) ]
      -ratio <+> [ Ratio of ECFP4:Daylight:AtomPair:MACCS Weighted Average (Def: 2 2 2 1) ]

e.g.> x.py -ref a.smi -que b.sdf.bz2 -rank ecfp4 -pref output -ratio 2 2 1 0 -dice
      return: output.txt (all scores); output.sdf.gz (new tag to -ref)
```
- Calculate chemical similarity coefficient of 2 sets of input molecules. Takes gzip files of both SDF (.sdf) and SMILES (.smi) formats.
- This version uses **single-cpu + Pandas** to calculate fp similarity. For comparison, a 767x850 pairwise set:
```
  2-cpu mpi-only version: ~ 170s (4 FPs)
  8-cpu mpi-only version: ~ 145s (4 FPs)
  1-cpu Pandas version:   ~  40s (4 FPs)
```
- **MPI-only** is simply inefficient for this type of *many-small-but-fast* job. *Overhead* to create the many processes is too much.

-  *Dask* is **mpi-cpu + Pandas** but also not good for this job. partitions=2 is ~6x slower than single-cpu; partitions=8 is ~10x slower! Probably has to do with number of small and fast jobs to perform where mpi overhead is too costly for this type of job. For dask, need to do: x.apply(x, meta=pd.Series(dtype='int',column='ecfp4'))

- Return output files: **prefix.txt** = pairwise query-reference similarity coeff ; **prefix.sdf.gz** = 1 closest query molecule for each reference molecule

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
- **Depricated** -- superseded by **1_fp_2_compare.py** that uses *Pandas* to calculate and manage.
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
> 2_LIGSIFT_parse.py
      -q     3D Query/Reference molecule     (sdf, mol2)
      -db    List of 3D Molecule databases name, no extension
             (txt of SDF, SDF only)
      -cpu   Number of CPU to run            [default: 1]
      -o     Output prefix
      -rank  Output top ranking conformer    [default: 'ShapeSimPval']
          'ShapeTanimoto', 'ChemTanimoto', 'ShapeSim', 'ChemSim',
          'ShapeSimPval', 'ChemSimPval', 'TverskyShape', 'TverskyChem', 
          'TverskyChem'\n
   e.g.>  x.py  -q=ref.mol2  -db=db.list        -cpu=4 
                -o=result    -rank=ShapeSimPval

result >
```
- This uses [LIGSIFT](https://doi.org/10.1093/bioinformatics/btu692) program to perform a ligand-based 3D conformer screening.
- To a 3D structure of a query/reference molecule, screen a library of molecules' 3D conformer to find those with similar shape and chemical features. Similar idea to OpenEye's [ROCS](https://www.eyesopen.com/rocs)+[EON](https://www.eyesopen.com/eon). Will need a pre-generated 3D conformer library from either OpenEye [OMEGA](https://www.eyesopen.com/omega), Schrödinger [ConfGen](https://www.schrodinger.com/confgen) or a RDKit [ETKDG](https://doi.org/10.1021/acs.jcim.5b00654) and my [script](https://github.com/mungpeter/Structure-Based_docking/tree/master/A_docking_scripts/3_conformer_gen), and perhaps others (Conformator [A](https://doi.org/10.1021/acs.jcim.8b00704) [B](https://www.zbh.uni-hamburg.de/forschung/amd/software/conformator.html) , [CSD](https://doi.org/10.1021/acs.jcim.7b00697) ).


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
LIGSIFT           # an open-source tool for ligand structural alignment and virtual screening
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
