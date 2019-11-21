# putils: A package for personal utility functions.
## Description
The putils package provides personal utility functions.

## Installation
```r
library(devtools)
install_github('wangtengyao/putils')
```

## List of all functions 
### Mathematical functsions
* expit
* logit
* lambertW

### Matrices and vectors
* vector.norm
* vector.normalise
* vector.clip
* vector.soft.thresh
* vector.hard.thresh
* vector.scramble
* vector.cleanup
* matrix.GramSchmidt
* matrix.trace
* matrix.rank
* matrix.standardise
* matrix.power
* sinThetaLoss
* powerMethod
* ql
* offdiag

### Random element generation
* random.rademacher
* random.bernoulli
* random.UnitVector
* random.OrthogonalMatrix
* random.WishartMatrix
* random.WignerMatrix
* random.psdMatrix
* random.SymmetricMatrix

### Auxiliary functions
* printPercentage
* visualise
* snippet
* find.first
* find.last
* sf
* sf_exp
* sim.params
* show.params
* '%=%'

### Statistical functions
* CvM.test

### String operations
* strlen
* strstr
* prefix
* suffix

### NA handling
* setNA


