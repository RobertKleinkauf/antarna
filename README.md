# antaRNA---ant assembled RNA

antaRNA solves the RNA inverse folding problem for RNA on the secondary structure level by applying 
ant-colony optimization meta-heuristics. The tool mimics the behavior of ants when they were foraging.
By doing so, a solution sequence is found, which satisfies the user defined constraints.

Starting with version 2.0.0 bistable RNA can be modelled. This can either be done within one dotplot for intrinsically oszillating RNA or in two conformations dotplots simulating the binding event of an ligand, which should induce structural rearrangement.

Please consult repsective release description for more info.

### VERSION HISTORY
 - v220(with Version Tag):
  - Introducing bistable inverse folding based on the partition function and their dotplots. General restructuring of the 	terrain graph due to the new structural obstacles.
  - new updating and bonification
  - new sequence emission routine
 - v118(with Version Tag): 
  - Complete Error Handling via the Variables.error. If there is an error, the program will not do the exit(1),
  but report to that variable. The user has to check manually. In command line mode, this is done automatically.
 - v117(with Version Tag):
  - Error handling is moved into the check() function of Variables class. Variables.error can be used as checker.
 - v116(with Version Tag):
  - adding the AntHill class. A Variables object is instanciated by the call of that class
 - v115(with Version Tag):
  - adding the Variables class
 - v114(with Version Tag):
  - additionally adding HotKnots and IPKnot into structure predictors such that the calculations can be done with them too
 - v113:
  - multiple target GC values in different areas of the RNA
 - v112:
  - including soft sequence constraint: lower case letters invoke the terrain at this position to be possible for all nucleotides. The position is then allocated according to the seq distance penalty within the quality measure.
  - improved reachableGC measure: now determines a minimum and a maximum rachable GC interval. Exits the program if	the entered tGC is outside this interval.
 - v111:
  - including LP management [1/0]
  - including strict and soft fuzzy structure constraints
   -Lower case: soft fuzzy constraint: all structure is true, no structure is explicitly enforced
   -Upper case: hard fuzzy constraint: at least one BP is required within one definition block in order to not rise a penalty
 - v110:
  - parametrized pseudoknot webserver related paper version
 - v109:
  - has now some extended RNAfold functionality (parameterfile selection, Temperature)
	



###### References
When using antaRNA in your work, please cite antaRNA.
 - R. Kleinkauf, M. Mann and R. Backofen; antaRNA: ant colony-based RNA sequence design,doi:10.1093/bioinformatics/btv319, Bioinformatics, 2015 Oct 1;31(19):3114-21
  - [antaRNA: ant colony-based RNA sequence
design](http://bioinformatics.oxfordjournals.org/content/31/19/3114.long).
 - R. Kleinkauf, T. Houwaart, R. Backofen and M. Mann; antaRNA: Multi-Objective Inverse Folding of Pseudoknot RNA using Ant-Colony Optimization;doi:10.1186/s12859-015-0815-6, 16:389, BMC Bioinformatics(2015)
  - [antaRNA: Multi-Objective Inverse Folding of Pseudoknot RNA using Ant-Colony Optimization](http://www.biomedcentral.com/content/pdf/s12859-015-0815-6.pdf)
