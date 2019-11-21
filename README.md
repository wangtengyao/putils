# putils: A package for personal utility functions.
## Description
The putil package provides personal utility functions.

### Mathematical functsions
expit
logit
lambertW

### Matrices and vectors
vector.norm
vector.normalise
vector.clip
vector.soft.thresh
vector.hard.thresh
vector.scramble
vector.cleanup
matrix.GramSchmidt
matrix.trace
matrix.rank
matrix.standardise
matrix.power
sinThetaLoss
powerMethod
ql
offdiag

### random element generation
random.rademacher
random.bernoulli
random.UnitVector
random.OrthogonalMatrix
random.WishartMatrix
random.WignerMatrix
random.psdMatrix
random.SymmetricMatrix

### auxiliary functions
printPercentage
visualise
snippet
find.first
find.last
sf
sf_exp
sim.params
show.params
'%=%'

### statistical functions
CvM.test

### string operations
strlen
strstr
prefix
suffix

### NA handling
setNA










# To create a new package:
library(devtools)
setwd('parent_dir_of_package')
package.skeleton('putils')
setwd('putils')
# copy the source file(s) to ./putils/R, delete NAMESPACE
# annotate source file with roxygen2 comments
document()
package <- as.package('.')
load_all(package)
# open terminal to the parent folder of package and execute
# R CMD build putils



## To recompile this package:
# Open RStudio
setwd('~/Dropbox/Programming/putils')
library(devtools)
document()
package <- as.package('.')
load_all(package)

# open terminal to ~/Dropbox/Programming
R CMD build putils