# antaRNA---ant assembled RNA

antaRNA solves the RNA inverse folding problem for RNA on the secondary structure level by applying 
ant-colony optimization meta-heuristics. The tool mimics the behavior of ants when they were foraging.
By doing so, a solution sequence is found, which satisfies the user defined constraints.

## Description

antaRNA comes in a single python file so far.
A suitabe packaging and maybe as well introduction of class wrappers in sighted.


### VERSION HISTORY (before checking)
 - v118: 
  - Complete Error Handling via the Variables.error. If there is an error, the program will not do the exit(1),
  but report to that variable. The user has to check manually. In command line mode, this is done automatically.
 - v117:
  - Error handling is moved into the check() function of Variables class. Variables.error can be used as checker.
 - v116:
  - adding the AntHill class. A Variables object is instanciated by the call of that class
 - v115:
  - adding the Variables class
 - v114:
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
	


## Constraints
### Structure Constraint
Structure constraints are standard input to antaRNA. They should be provided in regular dot-bracket notation.
 - regular
  - nested structures, e.g. `...(((...)))...`
  - pseudo-knot structures, e.g. `((((.[[[[...)))).]]]]`
 - *fuzzy*
  - e.g. `AAAA(((...)))AAAA`, s.t. each base pair within the `A` block is counted legal. At least one base pair has to occure in the solution sequence.
  - e.g. `aaaa(((...)))aaaa`, s.t. each base pair within the `a` block is counted legal. No base pair has to be there, but of course can occure in the solution.

### Sequence Constraint
Sequence constraint is hold in IUPAC RNA nucleotide definitions and can be set to 'soft mode', by writing the constraint in lower case

### GC contraint
GC contraint can be concepted in various ways:
 - singular or multiple target GC definitions wihtin one design for non-overlapping stretches of the sequence.
 - GC sampling, such that either normal or gaussian distributions are populated

###### References
When using antaRNA in your work, please cite antaRNA.
 - R. Kleinkauf, M. Mann and R. Backofen, antaRNA: ant colony-based RNA sequence design,doi:10.1093/bioinformatics/btv319, Bioinformatics, 2015
  - [antaRNA: ant colony-based RNA sequence
design](http://bioinformatics.oxfordjournals.org/content/early/2015/06/24/bioinformatics.btv319.full.pdf+html).
 - R. Kleinkauf, T. Houwaart, R. Backofen and M. Mann, antaRNA – Multi-Objective Inverse Folding of Pseudoknot RNA using Ant-Colony Optimization 
  - Submitted: Jul 2015