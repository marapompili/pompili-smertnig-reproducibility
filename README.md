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

The implementation relies on the SageMath development version containing
the Banff cluster algebra and finite Laurent intersection ring (FLIR)
implementation introduced in SageMath pull request #42538.

---

## Repository contents

- `performance_comparison.py`
  
  Reproduces the performance comparison between the FLIR algorithm and the
  presentation/primary decomposition approach.

- `factorization_examples.py`
  
  Reproduces the factorization examples appearing in the paper.

---

## Requirements

This repository **does not work with the current official SageMath release**.

It requires a development version of SageMath containing the implementation
introduced in SageMath pull request **#42538**.

In particular, the scripts use

```python
from sage.algebras.banff_cluster_algebra import BanffClusterAlgebra
```

which is not yet available in an official SageMath release.

---

## Installation

Clone SageMath and check out the branch containing pull request #42538:

```bash
git clone https://github.com/sagemath/sage.git
cd sage
git fetch origin pull/42538/head:pr-42538
git checkout pr-42538
```

Build SageMath following the official developer installation instructions.

---

## Running the examples

After building SageMath, execute

```bash
/path/to/sage performance_comparison.py
```

to reproduce the benchmark computations.

To reproduce the factorization examples, execute

```bash
/path/to/sage factorization_examples.py
```

