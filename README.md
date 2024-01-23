# Mutation utilities for protein design

A very raw version, work in progress. To install, run `pip install git+https://github.com/anton-bushuiev/mutils.git`.

## `Mutation` class

```python
from mutils.mutations import Mutation

# Parse a double-point mutation
mutation = Mutation('TC13A,GC13aA')
mutation
> Mutation(TC13A,GC13aA)

# Revert
mutation.revert()
> Mutation(AC13T,AC13aG)

# Check for the presence of insertions
mutation.insertion()
> [False, True]

# Convert wild types to Graphein format
mutation.wt_to_graphein()
> ['THR:C:13', 'GLY:C:13:a']
```

## `MutationSpace` class

```python
from mutils.mutations import MutationSpace

# Define possible substitutions
space = MutationSpace({'A10': 'AG', 'A11': 'A', 'A12': ''})
space
> MutationSpace({'A10': 'AG', 'A11': 'A', 'A12': ''})

# Construct all single-point mutations in the space
space.construct(d=1)
> ['A10A', 'A10G', 'A11A']

# Construct all double-point mutations in the space
space.construct(d=2)
> ['A10A,A11A', 'A10G,A11A']

# Construct all mutations in the space
space.construct()
> ['A10A,A11A', 'A10A', 'A10G,A11A', 'A10G', 'A11A']

# Get the size of the space without constructing
space.size(d=2)
> 2
```

## Utilities for `.pdb` and `.fasta` files

```python
from mutils.pdb import download_pdb, pdb2fasta
from mutils.proteins import parse_fasta

# Download structure from PDB
download_pdb('1bui')

# Convert sequences to FASTA
fasta = pdb2fasta('1bui.pdb', verbose=False)
fasta
> '>1BUI:A\nPSFDCGKPQVEPKKCXPGXVVGGCVAHPHSWPWQVSLRTRFGMHFCGGTLISPEWVLTAAHCLEKSPRPSSYKVILGAHQEVNLEPHVQEIEVSRLFLEPTRXXXXXXKDIALLKLSSPAVITDKVIPACLPSPNYVVADRTECFITGWGETQGXXTFGAGLLKEAQLPVIENKVCNRYEFLNGRVQSTELCAGHLAGGTDSCQGDSGGPLVCFEKDKYILQGVTSWGLXGCARPNKPGVYVRVSRFVTWIEGVMRNN\n>1BUI:B\nAPSFDCGKPQVEPKKCXPGXVVGGCVAHPHSWPWQVSLRTRFGMHFCGGTLISPEWVLTAAHCLEKSPRPSSYKVILGAHQEVNLEPHVQEIEVSRLFLEPTRXXXXXXKDIALLKLSSPAVITDKVIPACLPSPNYVVADRTECFITGWGETQGXXTFGAGLLKEAQLPVIENKVCNRYEFLNGRVQSTELCAGHLAGGTDSCQGDSGGPLVCFEKDKYILQGVTSWGLXGCARPNKPGVYVRVSRFVTWIEGVMRNN\n>1BUI:C\nASYFEPTGPYLMVNVTGVDSKGNELLSPHYVEFPIKPGTTLTKEKIEYYVEWALDATAYKEFRVVELDPSAKIEVTYYDKNKKKEETKSFPITEKGFVVPDLSEHIKNPGFNLITKVVIEKK\n'

# Get sequences as a dict
parse_fasta(fasta)
> {'A': 'PSFDCGKPQVEPKKCPGVVGGCVAHPHSWPWQVSLRTRFGMHFCGGTLISPEWVLTAAHCLEKSPRPSSYKVILGAHQEVNLEPHVQEIEVSRLFLEPTRKDIALLKLSSPAVITDKVIPACLPSPNYVVADRTECFITGWGETQGTFGAGLLKEAQLPVIENKVCNRYEFLNGRVQSTELCAGHLAGGTDSCQGDSGGPLVCFEKDKYILQGVTSWGLGCARPNKPGVYVRVSRFVTWIEGVMRNN',
>  'B': 'APSFDCGKPQVEPKKCPGVVGGCVAHPHSWPWQVSLRTRFGMHFCGGTLISPEWVLTAAHCLEKSPRPSSYKVILGAHQEVNLEPHVQEIEVSRLFLEPTRKDIALLKLSSPAVITDKVIPACLPSPNYVVADRTECFITGWGETQGTFGAGLLKEAQLPVIENKVCNRYEFLNGRVQSTELCAGHLAGGTDSCQGDSGGPLVCFEKDKYILQGVTSWGLGCARPNKPGVYVRVSRFVTWIEGVMRNN',
>  'C': 'ASYFEPTGPYLMVNVTGVDSKGNELLSPHYVEFPIKPGTTLTKEKIEYYVEWALDATAYKEFRVVELDPSAKIEVTYYDKNKKKEETKSFPITEKGFVVPDLSEHIKNPGFNLITKVVIEKK'}
```

## Reading preprocessed datasets

Please cite the corresponding papers if you find datasets useful (see `data/README.md`).

```python
from mutils.data import load_SKEMPI2

# Load and preprocess SKEMPI2 dataset
df = load_SKEMPI2()[0]
df
             #Pdb Mutation(s)_PDB Mutation(s)_cleaned iMutation_Location(s) Hold_out_type  Hold_out_proteins  ...     dG_mut      dG_wt       ddG  PDB Id Partner 1 Partner 2
0        1CSE_E_I           LI45G               LI38G                   COR         Pr/PI              Pr/PI  ... -14.022334 -16.302911  2.280577    1CSE         E         I
1        1CSE_E_I           LI45S               LI38S                   COR         Pr/PI              Pr/PI  ... -15.114136 -16.302911  1.188776    1CSE         E         I
2        1CSE_E_I           LI45P               LI38P                   COR         Pr/PI              Pr/PI  ...  -9.537466 -16.302911  6.765446    1CSE         E         I
3        1CSE_E_I           LI45I               LI38I                   COR         Pr/PI              Pr/PI  ... -13.320410 -16.302911  2.982502    1CSE         E         I
4        1CSE_E_I           LI45D               LI38D                   COR         Pr/PI              Pr/PI  ... -11.891069 -16.302911  4.411843    1CSE         E         I
...           ...             ...                 ...                   ...           ...                ...  ...        ...        ...       ...     ...       ...       ...
7080  3QIB_ABP_CD            KP9R                KP8R                   COR      TCR/pMHC  TCR/pMHC,1JCK_A_B  ...  -4.938011  -7.175045  2.237034    3QIB       ABP        CD
7081  3QIB_ABP_CD           TP12A               TP11A                   COR      TCR/pMHC  TCR/pMHC,1JCK_A_B  ...  -4.036047  -7.175045  3.138999    3QIB       ABP        CD
7082  3QIB_ABP_CD           TP12S               TP11S                   COR      TCR/pMHC  TCR/pMHC,1JCK_A_B  ...  -6.099323  -7.175045  1.075723    3QIB       ABP        CD
7083  3QIB_ABP_CD           TP12N               TP11N                   COR      TCR/pMHC  TCR/pMHC,1JCK_A_B  ...  -5.951210  -7.175045  1.223835    3QIB       ABP        CD
7084  3QIB_ABP_CD      YP7F,TP12S          YP6F,TP11S               COR,COR      TCR/pMHC  TCR/pMHC,1JCK_A_B  ...  -5.958076  -7.175045  1.216970    3QIB       ABP        CD

[6706 rows x 35 columns]
```
