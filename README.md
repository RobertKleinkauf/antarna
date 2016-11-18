# antaRNAdp - ant assembled RNA: Designing bistable RNA

## Description

antaRNA is extended towards the usage of the ViennaRNA Tools python wrappers to 
get started with base pair probability matrices in order to calculate and model 
bistable RNA molecules and mimicking ligand interaction by applying 
constraints to simulate the induced conformation.

## Usage

antaRNA can be used in several ways. 

#### Regular command line example:
```
$python antaRNA.py -tGC 0.5 -Cstr "...(((...)))..." -ov
>antaRNA#0|Cstr:...((((...))))...|Cseq:NNNNNNNNNNNNNNNNN|Alpha:1.0|Beta:1.0|tGC:1-17>0.5|ER:0.2|Struct_CT:0.5|GC_CT:5.0|Seq_CT:1.0|UsedProgram:RNAfold|Ants:2|Resets:0/5|AntsTC:50|CC:130|IP:s|BSS:0|LP:0|ds:0.0|dGC:0.0|GC:52.9411764706|dseq:0.0|L:17|Time:0.0177099704742
Rseq:GAAACACGAGGUGUAGG
Rstr:...((((...))))...
```

#### Python object
```python
import antaRNA # importing antaRNA

anthill = antaRNA.AntHill() # creating an anthill object

anthill.params.tGC.append(0.5) # setting the generat target GC value to 0.5
anthill.params.Cstr = "...(((...)))..." # setting structural constraint
anthill.params.output_verbose = True # activating verbose output 
anthill.params.check() # checking the input

anthill.swarm() if anthill.params.error == "0" else anthill.params.error # execution if no error was detected
anthill.result # retrieve result
[('>antaRNA#0', 'Cstr:...(((...)))...', 'Cseq:NNNNNNNNNNNNNNN', 'Alpha:1.0', 'Beta:1.0', 'tGC:1-15>0.5', 'ER:0.2', 'Struct_CT:0.5', 'GC_CT:5.0', 'Seq_CT:1.0', 'UsedProgram:RNAfold', 'Ants:7', 'Resets:0/5', 'AntsTC:50', 'CC:130', 'IP:s', 'BSS:0', 'LP:0', 'ds:0.0', 'dGC:0.0', 'GC:46.6666666667', 'dseq:0.0', 'L:15', 'Time:0.0780849456787', 'Rseq:UCACGUCUUACGAAG', 'Rstr:...(((...)))...')]
```





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

