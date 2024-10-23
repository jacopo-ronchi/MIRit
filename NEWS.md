# MIRit 1.1.1

Changes were made to ensure that R CMD check runs without errors or warnings.
Moreover, MIRit has been updated to use the new miRTarBase version 10 database.

# MIRit 0.99.13

A minor fix was made to fix the undefined variables note during R CMD check.

# MIRit 0.99.12

This version includes several improvements, including a completely revised
vignette where all chunks are evaluated, minor tweaks to default values for
differential expression analysis, and some bug fixes to the error bars in the
`plotDE()` function. Other issues, such as artifacts in show methods, lacks in
documentation, and dependence in DESCRIPTION, have been addressed too.

# MIRit 0.99.11

This new version introduces the possibility of limiting validated targets
retrieval from miRTarBase to only those interactions supported by extensive
experimental evidence. Moreover, minor fixes were made to the vignette,
documentation, and internal data. Finally, some examples have been redefined
to reduce checking time.

# MIRit 0.99.10

This patch introduces a minor fix for one unit test.

# MIRit 0.99.9

With this version, significant p-values originating from functional enrichment
analyses now include extreme values. Further, examples have been shortened and
unit tests now use smaller datasets.

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
