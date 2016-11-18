# antaRNAdp - ant assembled RNA: Designing bistable RNA

## Description

antaRNA is extended towards the usage of the ViennaRNA Tools python wrappers to 
get started with base pair probability matrices in order to calculate and model 
bistable RNA molecules and mimicking ligand interaction by applying 
constraints to simulate the induced conformation.

## Usage

antaRNA can be used in several ways. This includes python object management or a direct command line call of the program.
The program differentiates between two modi: MFE (dot-bracket structure optimization) and DP (dot plot structure optimization). Each 
of the modi has specific options and individual means of controlling its input. However, both modi also share options and input. Therefore
a commandline call of the programm constitutes as:

```
$python antaRNA.py [COMMON INPUT AND OPTIONS] [MODUS{MFE,DP}] [MODUS INPUT AND OPTIONS]
```
See in Sources for detailed description of the explicit usage...

## VERSION HISTORY OF antaRNAdp

 - v2.0.1
  - Bioconda Version

 - v2.0.0 
  - new Terrain Graph and associated functions
  - new way of defining structural constraints


## SOURCES

### Algrithmic description and exemplary documentation

 - R. Kleinkauf: [Ant colony optimization based inverse folding of mono- and bistable RNA macromolecules](https://freidok.uni-freiburg.de/data/11095) DOI: 10.6094/UNIFR/11095, 2016

### References
When using antaRNA in your work, please cite one of the following, suiting your work.

 - R. Kleinkauf, M. Mann and R. Backofen; [antaRNA: ant colony-based RNA sequence
design](http://bioinformatics.oxfordjournals.org/content/31/19/3114.long),doi:10.1093/bioinformatics/btv319, Bioinformatics, 2015 Oct 1;31(19):3114-21

 - R. Kleinkauf, T. Houwaart, R. Backofen and M. Mann; antaRNA: [antaRNA: Multi-Objective Inverse Folding of Pseudoknot RNA using Ant-Colony Optimization](http://www.biomedcentral.com/content/pdf/s12859-015-0815-6.pdf);doi:10.1186/s12859-015-0815-6, 16:389, BMC Bioinformatics(2015)

