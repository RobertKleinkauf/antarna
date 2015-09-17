# antaRNA --- ant assembled RNA

antaRNA solves the RNA inverse folding problem for RNA on the secondary structure level by applying 
ant-colony optimization meta-heuristics. The tool mimics the behavior of ants when they were foraging.
By doing so, a solution sequence is found, which satisfies the user defined constraints.

#### Constraints
### Structure Constraint
Structure constraints are standard input to antaRNA. They should be provided in regular dot-bracket notation.
 - regular
  - nested structures, e.g. "...(((...)))..."
  - pseudo-knot structures, e.g. "((((.[[[[...)))).]]]]"
 - *fuzzy*
  - e.g. "AAAA(((...)))AAAA", s.t. each base pair within the `A` block is counted legal. At least one base pair has to occure in the solution sequence.
  - e.g. "aaaa(((...)))aaaa", s.t. each base pair within the `a` block is counted legal. No base pair has to be there, but of course can occure in the solution.

### Sequence Constraint
Sequence constraint is hold in IUPAC RNA nucleotide definitions and can be set to 'soft mode', by writing the constraint in lower case

### GC contraint
GC contraint can be concepted in various ways:
 - singular or multiple target GC definitions wihtin one design for non-overlapping stretches of the sequence.
 - GC sampling, such that either normal or gaussian distributions are populated

#### References
When using antaRNA in your work, please cite antaRNA.
 - R. Kleinkauf, M. Mann and R. Backofen, antaRNA: ant colony-based RNA sequence design,doi:10.1093/bioinformatics/btv319, Bioinformatics, 2015
  - [antaRNA: ant colony-based RNA sequence
design](http://bioinformatics.oxfordjournals.org/content/early/2015/06/24/bioinformatics.btv319.full.pdf+html).
 - R. Kleinkauf, T. Houwaart, R. Backofen and M. Mann, antaRNA – Multi-Objective Inverse Folding of Pseudoknot RNA using Ant-Colony Optimization 
  - Submitted: Jul 2015