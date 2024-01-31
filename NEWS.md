# MIRit 0.99.8

This patch further reduces R CMD check time by limiting unnecessary examples
and by using smaller datasets for unit tests.

# MIRit 0.99.7

This version introduces parallel computing capabilities for the
`mirnaIntegration()` function. This is particularly useful for Boschloo's
exact test, whose execution is now faster. Moreover, the test suite has been
redefined to reduce running times.

# MIRit 0.99.6

The testing suite has been redefined to allow different results for different
versions of packages employed in differential expression analysis.

# MIRit 0.99.5

This update fixes a bug in the `batchCorrection()` function that prevented the
correct use of this function with newer versions of the `MultiAssayExperiment`
package.

# MIRit 0.99.4

After the implementation of the IS_BIOC_BUILD_MACHINE variable to the Single
Package Builder (SPB), this version bump drives a new build to fix errors 
during R CMD check on SPB.

# MIRit 0.99.3

This version fixes a bug in the `topologicalAnalysis()` function that
prevented the use of a functional progress bar during permutation testing.
Moreover, the example for the `addDifferentialExpression()` function has been updated to reduce its running time.

# MIRit 0.99.2

MIRit now allows to filter the pathways used for topological analysis based
on the number of nodes.

# MIRit 0.99.1

Functional enrichment analyses and TAIPA now use cached databases to reduce
running times.

# MIRit 0.99.0

Initial version for Bioconductor submission.
