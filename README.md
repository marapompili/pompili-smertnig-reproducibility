# Reproducibility files for Pompili–Smertnig computations

This repository contains SageMath scripts for independently checking
computations appearing in:

M. Pompili and D. Smertnig,
*Factoriality and Class Groups of Upper Cluster Algebras and Finite
Laurent Intersection Rings: A Computational Approach*,
arXiv:2601.07520.

The scripts verify:

1. the performance-comparison table;
2. two factorization examples.

The computations use the `FLIR` and `BanffClusterAlgebra`
implementations proposed in SageMath pull request #42538.

## Reproducibility status

The code in this repository is intended to make every input matrix,
mutation sequence, element, timeout and reported invariant explicit.
Benchmark timings are machine-dependent; algebraic outputs should be
independent of the machine and are checked against committed expected
results.
