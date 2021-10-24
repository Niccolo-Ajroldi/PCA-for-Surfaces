
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

# Multivariate PCA -------------------------------------------------------------

Xt = out$Xt
is(Xt) # list of matrices
length(Xt)
n_points = length(out$grid.1) * length(out$grid.2)
n_points # total number of points in the grid

## Vectorize Xt ----

Xt_mat = vapply(Xt, as.numeric, numeric(length=n_points))
dim(Xt_mat) # n_points x n

# discretization step
w = 1/(n_points-1)

# pca
pca = prcomp(t(Xt_mat), center=FALSE, scale.=FALSE)

# the number of FPC's is found by cumulative proportion of variance
nharm = 20
print(cumsum(pca$sd^2)/sum(pca$sd^2))

# eigenvalues of FPCA
lambda = (pca$sdev^2)[1:nharm] * w # w needed to get back to the functional form

# projection of training data on FPC's
scores.train = pca$x[,1:nharm] * sqrt(w) # n x nharm

# FPC's evaluated on the time grid
efpc_eval = pca$rotation[,1:nharm] / sqrt(w) # length(grid) x nharm
dim(efpc_eval)

rm(pca)

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








