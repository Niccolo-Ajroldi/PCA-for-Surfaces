
#' Functional Principal component estimation for functional data defined on a 
#' bivariate domain.

rm(list=ls())

source("D:/Poli/TESI/Code/Time-Series-CP/FAR_2D/Simulation/simulate_FAR.R")

# SIMULATION PARAMETERS
nbasis.1.sim = 5
nbasis.2.sim = 5
sample_size = 50
x1.grid = seq(from=0, to=1, length=50)
x2.grid = seq(from=0, to=1, length=50)

# seed
simulation_seed = 1

# Data -------------------------------------------------------------------------

# total number of basis functions for the simulation (tensor product basis)
nbasis.sim = nbasis.1.sim*nbasis.2.sim

# Parameters
sample_size = sample_size

## FAR(1) non-concurrent parameters ----

# d = dimension of the underlying VAR(1)
d=nbasis.sim 
Psi1=matrix(0.3,d,d)
diag(Psi1)=0.8
Psi1=Psi1/norm(Psi1, type = "F")/2
Psi=array(0,c(d,d,1))
Psi[,,1]=Psi1
my_Sigma=matrix(0.6,d,d)
diag(my_Sigma)=1
my_Sigma=my_Sigma/2

# set seed
set.seed(simulation_seed)

## simulate data ----
out <- simulate_FAR(n = sample_size, 
                    Psi = Psi, 
                    x1.grid = x1.grid, 
                    x2.grid = x2.grid, 
                    nbasis.1 = nbasis.1.sim, 
                    nbasis.2 = nbasis.2.sim, 
                    sigma = my_Sigma,
                    basis.type = "bspline")

rm(nbasis.1.sim, nbasis.2.sim, nbasis.sim)

# Basis definition -------------------------------------------------------------

# grid: [0,1]x[0,1]
x1.grid <- out$grid.1
x2.grid <- out$grid.2
length(x2.grid)

# number of basis functions in each dimension
nbasis.1 = 5
nbasis.2 = 5
nbasis = nbasis.1*nbasis.2

# build 1D basis functions on each of the 1D domains
basis.1 <- fda::create.bspline.basis(rangeval=range(x1.grid), nbasis=nbasis.1)
basis.2 <- fda::create.bspline.basis(rangeval=range(x2.grid), nbasis=nbasis.2)

# evaluate 1D basis functions on the grid
basis.1.eval <- fda::eval.basis(evalarg=x1.grid, basisobj=basis.1)
basis.2.eval <- fda::eval.basis(evalarg=x2.grid, basisobj=basis.2)

dim(basis.1.eval)
dim(basis.2.eval)

# 2D basis functions evaluated on the grid
basis.grid.eval <- vector("list", length = nbasis)
grid.2D <- expand.grid(x1.grid, x2.grid)
counter = 1

# tensor product basis construction
for(i in 1:nbasis.1)
{
  base_i_j = matrix(nrow=length(x1.grid),ncol=length(x2.grid))
  for(j in 1:nbasis.2)
  {
    # outer product sui punti della griglia
    base_i_j = outer(basis.1.eval[,i],basis.2.eval[,j])
    
    basis.grid.eval[[counter]] = base_i_j
    counter = counter+1
  }
}
rm(counter)

# Project data on the basis ----------------------------------------------------

Xt = out$Xt
is(Xt) # list of matrices
length(Xt)
n_points = length(out$grid.1) * length(out$grid.2)
n_points # total number of points in the grid

## Vectorize Xt ----

Xt_mat = vapply(Xt, as.numeric, numeric(length=n_points))
dim(Xt_mat) # n_points x n

### centering ----

mean.functional = rowMeans(Xt_mat)
Xt_mat = sweep(Xt_mat, 1, mean.functional)

## Vectorize basis -----

basis_eval = vapply(basis.grid.eval, as.numeric, numeric(length=n_points))
dim(basis_eval) # n_points x nbasis

## Project on the basis ----

# basis_eval %*% C = Xt_vector 
# crossprod(basis_eval) %*% C = crossprod(basis_eval, Xt_mat)
# B %*% C = X

B = crossprod(basis_eval)
X = crossprod(basis_eval, Xt_mat)
C = solve(B,X)

dim(C) # nbasis x n

## check ----
Xt_reconstructed = basis_eval %*% C
all.equal(Xt_reconstructed, Xt_mat)

## COMPUTE W ! -----------------------------------------------------------------

# Let's define a dicretization grid for numerical integration
x1.discretized = seq(0,1,by=0.001)
x2.discretized = seq(0,1,by=0.001)

# n_points of dicretization grid
n_points.discretized = length(x1.discretized) * length(x2.discretized)

# evaluate 1D basis on such grid
basis.1.eval.discretized <- fda::eval.basis(evalarg=x2.discretized, basisobj=basis.1)
basis.2.eval.discretized <- fda::eval.basis(evalarg=x2.discretized, basisobj=basis.2)

# 2D basis functions evaluated on the bidimensional grid
basis.grid.eval.discretized <- vector("list", length = nbasis)
grid.2D.discretized <- expand.grid(x1.discretized, x2.discretized)
counter = 1

# considero ogni combinazione di basi
for(i in 1:nbasis.1){
  base_i_j = matrix(nrow=length(x1.discretized),ncol=length(x2.discretized))
  for(j in 1:nbasis.2){
    # outer product sui punti della griglia
    base_i_j = outer(basis.1.eval.discretized[,i],basis.2.eval.discretized[,j])
    basis.grid.eval.discretized[[counter]] = base_i_j
    counter = counter+1
  }
}
rm(counter)

# Vectorize basis
basis_eval.discretized = vapply(basis.grid.eval.discretized, as.numeric, numeric(length=n_points.discretized))

# check orthogonality
w_discretized = 1/( (length(x1.discretized)-1) * (length(x2.discretized)-1))
W_discretized = crossprod(basis_eval.discretized, basis_eval.discretized) * w_discretized
W = W_discretized

# W is diagonal only if the basis system is orthonormal
# this happens  for example when we use as basis system the tensor
# product of two fourier basis systems. If we instead use the product
# of bspline basis, W will not be diagonal.
all.equal(W, diag(nbasis))

# Estimate FPCS ----------------------------------------------------------------

# library for matrix square root
library(expm) 

eigen_dec = eigen(expm::sqrtm(W) %*% tcrossprod(C)%*% expm::sqrtm(W)/sample_size, symmetric=TRUE)
lambda = eigen_dec$values
efpc = eigen_dec$vectors # vectors, efpc[i,j] = <x_i, x_j>

dim(efpc) # nbasis x nbasis
efpc = solve(expm::sqrtm(W)) %*% efpc

# cumulative proportion of explained variance
cumsum(lambda)/sum(lambda)

# select a number of harmonics (M in the paper)
nharm = nbasis

# subset eigenvalues and eigenvectors
lambda = lambda[1:nharm]
efpc = efpc[,1:nharm]
dim(efpc) # nbasis x nharm 

# Evaluate FPCS on the 2D grid -------------------------------------------------

efpc_eval = basis_eval %*% efpc
dim(efpc_eval) #  npoints x nharm

# Project data onto the FPCS ---------------------------------------------------

# C_efpc: matrix of coefficients of basis projection of X onto efpc, nharm x n
#
# efpc_eval %*% C_efpc = Xt_mat
# crossprod(efpc_eval) %*% C_efpc = crossprod(efpc_eval, Xt_mat)
# B %*% C_efpc = X

B = crossprod(efpc_eval)
X = crossprod(efpc_eval, Xt_mat)
C_efpc = solve(B,X)
dim(C_efpc) # nharm x n

## check ----
Xt_reconstructed = efpc_eval %*% C_efpc
all.equal(Xt_reconstructed, Xt_mat) # TRUE solo se metto nbasis>=nbasis.sim e pure le stesse basi!









