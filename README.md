## SLCFold
The source code is the implementation of the fitness score, which was designed for new fold. We successfully found a new fold SLC44 according the the fitness score.

### Data
#### Predicted information:
1. 3D structures from multiple prediction methods. RaprorX-3DModeling and trRosetta were used in the paper. 
2. Contact maps. TripleRes was adopted. 
3. Secondary structures. DeepCNF was adopted. 

Files in `input` directory:
```
├── SLC44A1                           # the name of target protein
│   ├── SLC44A1.fold.229863005.pdb    # predicted structure by RaptorX-3DModeling
│   ├── SLC44A1.partial.contact       # predicted contact by TripletRes
│   ├── SLC44A1.partial.fasta         # sequence
│   ├── SLC44A1.partial.pdb           # predicted structure by trRosetta
│   └── SLC44A1.ss3                   # secondary structure by DeepCNF
└── SLC6A12
    ├── SLC6A12.fold.361676800.pdb
    ├── SLC6A12.partial.contact
    ├── SLC6A12.partial.fasta
    ├── SLC6A12.partial.pdb
    └── SLC6A12.ss3
```
Note: The suffixes of the three kinds of predicted information need to be set as "pdb", "contact", and "ss3".

#### Experimental structures
Structure from RCSB PDB in PDB format. Here, only a few pdb files are included in *PBDs*. All membrane protein structures considered in the paper are listed in *memPDBs*.

#### TMalign and TMscore
The tool TMalign and TMscore are required.
And the corresponding paths of executables should be set in `scripts/ranking_score.py`.

#### Set the number of cpu
The variable *NUMPROCESS* in `scripts/ranking_score.py` can be set according to the number of cpu availiable in your machine.

### Environment
The code is implemented in python3.6.
The required packages are as follows:
```
glob
pandas
numpy
biopython
```

### Running
```
cd scripts/
python ranking_score.py --input ../input/ --term   # get each term of the score for all proteins
python ranking_score.py --input ../input/ --score  # get the score values for all proteins
```
The results are in `output`.

### Reference
Tengyu Xie, Ximin Chi, Bangdong Huang, Fangfei Ye, Qiang Zhou, Jing Huang. Rational Exploration of Fold Atlas for Human Solute Carrier Proteins. bioRxiv, 2021.08.21.457230; doi: https://doi.org/10.1101/2021.08.21.457230.

