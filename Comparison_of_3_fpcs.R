
#' Performo i due tipi di FPCA e poi faccio dei plot delle prime 3 armoniche nei due casi
#' TODO: fare una funzione che faccia FPCA 2D e snellire questo script!

rm(list=ls())
source("D:/Poli/TESI/Code/Time-Series-CP/FAR_2D/Simulation/simulate_FAR.R")

# SIMULATION PARAMETERS
nbasis.1.sim = 5
nbasis.2.sim = 5
sample_size = 50
x1.grid = seq(from=0, to=1, length=50)
x2.grid = seq(from=0, to=1, length=50)

# number of FPC's
nharm = 3

# seed
simulation_seed = 1 # 2,3,7 for best pics

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

# FPCA-basis-fourier ~~~~~~~~~ -------------------------------------------------

## Basis definition -------------------------------------------------------------

# grid: [0,1]x[0,1]
x1.grid <- x1.grid
x2.grid <- x2.grid
length(x2.grid)

# number of basis functions in each dimension
nbasis.1 = nbasis.1.sim
nbasis.2 = nbasis.2.sim
nbasis = nbasis.1*nbasis.2

# build 1D basis functions on each of the 1D domains
basis.1 <- fda::create.fourier.basis(rangeval=range(x1.grid), nbasis=nbasis.1)
basis.2 <- fda::create.fourier.basis(rangeval=range(x2.grid), nbasis=nbasis.2)

# evaluate 1D basis functions on the grid
basis.1.eval <- fda::eval.basis(evalarg=x1.grid, basisobj=basis.1)
basis.2.eval <- fda::eval.basis(evalarg=x2.grid, basisobj=basis.2)

dim(basis.1.eval)
dim(basis.2.eval)

# 2D basis functions evaluated on the grid
basis.grid.eval <- vector("list", length = nbasis)
grid.2D <- expand.grid(x1.grid, x2.grid)
counter = 1

# considero ogni combinazione di basi
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
rm(counter, base_i_j, grid.2D, i, j)

## Project data on the basis ----------------------------------------------------

Xt = out$Xt
is(Xt) # list of matrices
length(Xt)
n_points = length(out$grid.1) * length(out$grid.2)
n_points # total number of points in the grid

### Vectorize Xt ----

Xt_mat = vapply(Xt, as.numeric, numeric(length=n_points))
dim(Xt_mat) # n_points x sample_size

### centering ----
mean.functional = rowMeans(Xt_mat)
Xt_mat = sweep(Xt_mat, 1, mean.functional)

### Vectorize basis -----

basis_eval = vapply(basis.grid.eval, as.numeric, numeric(length=n_points))
dim(basis_eval) # n_points x nbasis

### Project on the basis ----

# basis_eval %*% C = Xt_vector 
# crossprod(basis_eval) %*% C = crossprod(basis_eval, Xt_mat)
# B %*% C = X

B = crossprod(basis_eval)
X = crossprod(basis_eval, Xt_mat)
C = solve(B,X)

dim(C) # nbasis x sample_size

### check ----
Xt_reconstructed = basis_eval %*% C
all.equal(Xt_reconstructed, Xt_mat) # TRUE solo se metto nbasis>=nbasis.sim e pure le stesse basi!

## Estimate FPCS ----------------------------------------------------------------

eigen_dec = eigen(tcrossprod(C)/sample_size, symmetric=TRUE)
lambda = eigen_dec$values
efpc = eigen_dec$vectors # vectors, efpc[i,j] = <x_i, x_j

dim(efpc) # nbasis x nbasis

lambda = lambda[1:nharm]

efpc = efpc[,1:nharm]
dim(efpc) # nbasis x nharm 

## Evaluate FPCS on the time grid -----------------------------------------------

efpc_eval = basis_eval %*% efpc
dim(efpc_eval) #  npoints x nharm

#matplot(efpc_eval, type='l')

## Project data onto the FPCS ---------------------------------------------------

# C_efpc: matrix of coefficients of basis projection of X onto efpc, nharm x sample_size
#
# efpc_eval %*% C_efpc = Xt_mat
# crossprod(efpc_eval) %*% C_efpc = crossprod(efpc_eval, Xt_mat)
# B %*% C_efpc = X

B = crossprod(efpc_eval)
X = crossprod(efpc_eval, Xt_mat)
C_efpc = solve(B,X)
dim(C_efpc) # nharm x sample_size

### check ----
Xt_reconstructed = efpc_eval %*% C_efpc
all.equal(Xt_reconstructed, Xt_mat) # TRUE solo se metto nharm=nbasis

## SAVE and REMOVE -------------------------------------------------------------

# vectorize efpc_eval
efpc_eval_list = list()
efpc_eval_list[[1]] = matrix(efpc_eval[,1], nrow=length(x1.grid), ncol=length(x2.grid))
efpc_eval_list[[2]] = matrix(efpc_eval[,2], nrow=length(x1.grid), ncol=length(x2.grid))
efpc_eval_list[[3]] = matrix(efpc_eval[,3], nrow=length(x1.grid), ncol=length(x2.grid))

# save efpc_eval_list for comparison
efpc_eval_list_BASIS = efpc_eval_list
lambda_BASIS = lambda
Xt_reconstructed_BASIS = Xt_reconstructed

rm(efpc_eval_list, Xt_reconstructed, C_efpc, X, B, efpc_eval, efpc, lambda, nharm, 
   eigen_dec, C, basis_eval, Xt_mat, n_points, Xt, basis.1, basis.2, basis.1.eval, basis.2.eval, nbasis.1, nbasis.2)


# FPCA-discretization ~~~~~~~~~ ------------------------------------------------

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

# discretization step
w = 1/(n_points-1)

# pca
pca = prcomp(t(Xt_mat), center=FALSE, scale.=FALSE)

# the number of FPC's is found by cumulative proportion of variance
nharm = 3
#nharm = as.integer(which.max(cumsum(pca$sd^2)/sum(pca$sd^2)>cum_prop_var))
#print(paste0("nharm: ", nharm))
print(cumsum(pca$sd^2)/sum(pca$sd^2))

# eigenvalues of FPCA
lambda = (pca$sdev^2)[1:nharm] * w # w needed to get back to the functional form

# projection of training data on FPC's
scores.train = pca$x[,1:nharm] * sqrt(w) # n x nharm

# FPC's evaluated on the time grid
efpc_eval = pca$rotation[,1:nharm] / sqrt(w) # length(grid) x nharm
dim(efpc_eval)

#rm(pca)

#matplot(efpc_eval, type='l')


## Project data onto the FPCS ---------------------------------------------------

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
all.equal(Xt_reconstructed, Xt_mat) # TRUE solo se metto nharm=nbasis

## SAVE and REMOVE -------------------------------------------------------------

# vectorize efpc_eval
efpc_eval_list = list()
efpc_eval_list[[1]] = matrix(efpc_eval[,1], nrow=length(x1.grid), ncol=length(x2.grid))
efpc_eval_list[[2]] = matrix(efpc_eval[,2], nrow=length(x1.grid), ncol=length(x2.grid))
efpc_eval_list[[3]] = matrix(efpc_eval[,3], nrow=length(x1.grid), ncol=length(x2.grid))

# save efpc_eval_list for comparison
efpc_eval_list_GRID = efpc_eval_list
lambda_GRID = lambda
Xt_reconstructed_GRID = Xt_reconstructed

rm(efpc_eval_list, Xt_reconstructed, C_efpc, X, B, efpc_eval, lambda, pca, w, Xt_mat, n_points, Xt)

# FPCA-basis-Bspline ~~~~~~~~~ -------------------------------------------------

## Basis definition ------------------------------------------------------------

# grid: [0,1]x[0,1]
x1.grid <- x1.grid
x2.grid <- x2.grid
length(x2.grid)

# number of basis functions in each dimension
nbasis.1 = 10
nbasis.2 = 10
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
rm(counter, base_i_j, grid.2D, i, j)

## Project data on the basis ----------------------------------------------------

Xt = out$Xt
is(Xt) # list of matrices
length(Xt)
n_points = length(out$grid.1) * length(out$grid.2)
n_points # total number of points in the grid

### Vectorize Xt ----

Xt_mat = vapply(Xt, as.numeric, numeric(length=n_points))
dim(Xt_mat) # n_points x sample_size

### centering ----

mean.functional = rowMeans(Xt_mat)
Xt_mat = sweep(Xt_mat, 1, mean.functional)

### Vectorize basis -----

basis_eval = vapply(basis.grid.eval, as.numeric, numeric(length=n_points))
dim(basis_eval) # n_points x nbasis

### Project on the basis ----

# basis_eval %*% C = Xt_vector 
# crossprod(basis_eval) %*% C = crossprod(basis_eval, Xt_mat)
# B %*% C = X

B = crossprod(basis_eval)
X = crossprod(basis_eval, Xt_mat)
C = solve(B,X)

dim(C) # nbasis x sample_size

### check ----
Xt_reconstructed = basis_eval %*% C
all.equal(Xt_reconstructed, Xt_mat) # TRUE solo se metto nbasis>=nbasis.sim e pure le stesse basi!

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

# W is not diagonal!
all.equal(W, diag(nbasis))

## Estimate FPCS ---------------------------------------------------------------

library(expm) # for matrix square root
eigen_dec = eigen(sqrtm(W) %*% tcrossprod(C)%*% sqrtm(W) /sample_size, symmetric=TRUE)
lambda = eigen_dec$values
efpc = eigen_dec$vectors # vectors, efpc[i,j] = <x_i, x_j>
dim(efpc)
efpc = solve(sqrtm(W)) %*% efpc

dim(efpc) # nbasis x nbasis

lambda = lambda[1:nharm]

efpc = efpc[,1:nharm]
dim(efpc) # nbasis x nharm 

## Evaluate FPCS on the time grid ----------------------------------------------

efpc_eval = basis_eval %*% efpc
dim(efpc_eval) #  npoints x nharm

#matplot(efpc_eval, type='l')

## Project data onto the FPCS ---------------------------------------------------

# C_efpc: matrix of coefficients of basis projection of X onto efpc, nharm x sample_size
#
# efpc_eval %*% C_efpc = Xt_mat
# crossprod(efpc_eval) %*% C_efpc = crossprod(efpc_eval, Xt_mat)
# B %*% C_efpc = X

B = crossprod(efpc_eval)
X = crossprod(efpc_eval, Xt_mat)
C_efpc = solve(B,X)
dim(C_efpc) # nharm x sample_size

### check ----
Xt_reconstructed = efpc_eval %*% C_efpc
all.equal(Xt_reconstructed, Xt_mat) # TRUE solo se metto nharm=nbasis

## SAVE and REMOVE -------------------------------------------------------------

# listerize efpc_eval
efpc_eval_list = list()
efpc_eval_list[[1]] = matrix(efpc_eval[,1], nrow=length(x1.grid), ncol=length(x2.grid))
efpc_eval_list[[2]] = matrix(efpc_eval[,2], nrow=length(x1.grid), ncol=length(x2.grid))
efpc_eval_list[[3]] = matrix(efpc_eval[,3], nrow=length(x1.grid), ncol=length(x2.grid))

# save efpc_eval_list for comparison
efpc_eval_list_BSPLINE = efpc_eval_list
lambda_BSPLINE = lambda
Xt_reconstructed_BSPLINE = Xt_reconstructed

rm(efpc_eval_list, Xt_reconstructed, C_efpc, X, B, efpc_eval, efpc, lambda, 
   eigen_dec, C, basis_eval, Xt_mat, n_points, Xt, basis.1, basis.2, basis.1.eval, basis.2.eval, nbasis.1, nbasis.2)


# RAPPORTO ---------------------------------------------------------------------

## GRID vs FOURIER ----

k_GRID_1 = mean(efpc_eval_list_BASIS[[1]]/efpc_eval_list_GRID[[1]])
k_GRID_2 = mean(efpc_eval_list_BASIS[[2]]/efpc_eval_list_GRID[[2]])
k_GRID_3 = mean(efpc_eval_list_BASIS[[3]]/efpc_eval_list_GRID[[3]])

k_GRID_1; k_GRID_2; k_GRID_3;

var(as.numeric(efpc_eval_list_BASIS[[1]]/efpc_eval_list_GRID[[1]]))
var(as.numeric(efpc_eval_list_BASIS[[2]]/efpc_eval_list_GRID[[2]]))
var(as.numeric(efpc_eval_list_BASIS[[3]]/efpc_eval_list_GRID[[3]]))

efpc_eval_list_GRID[[1]] = efpc_eval_list_GRID[[1]]*k_GRID_1
efpc_eval_list_GRID[[2]] = efpc_eval_list_GRID[[2]]*k_GRID_2
efpc_eval_list_GRID[[3]] = efpc_eval_list_GRID[[3]]*k_GRID_3

## GRID vs BSPLINE ----

k_BSPLINE_1 = mean(efpc_eval_list_BASIS[[1]]/efpc_eval_list_BSPLINE[[1]])
k_BSPLINE_2 = mean(efpc_eval_list_BASIS[[2]]/efpc_eval_list_BSPLINE[[2]])
k_BSPLINE_3 = mean(efpc_eval_list_BASIS[[3]]/efpc_eval_list_BSPLINE[[3]])

k_BSPLINE_1; k_BSPLINE_2; k_BSPLINE_3;

var(as.numeric(efpc_eval_list_BASIS[[1]]/efpc_eval_list_BSPLINE[[1]]))
var(as.numeric(efpc_eval_list_BASIS[[2]]/efpc_eval_list_BSPLINE[[2]]))
var(as.numeric(efpc_eval_list_BASIS[[3]]/efpc_eval_list_BSPLINE[[3]]))

efpc_eval_list_BSPLINE[[1]] = efpc_eval_list_BSPLINE[[1]]*k_BSPLINE_1
efpc_eval_list_BSPLINE[[2]] = efpc_eval_list_BSPLINE[[2]]*k_BSPLINE_2
efpc_eval_list_BSPLINE[[3]] = efpc_eval_list_BSPLINE[[3]]*k_BSPLINE_3

# MSE --------------------------------------------------------------------------

MLmetrics::MSE(efpc_eval_list_BASIS[[1]],efpc_eval_list_BSPLINE[[1]])
MLmetrics::MSE(efpc_eval_list_BASIS[[1]],efpc_eval_list_GRID[[1]])

MLmetrics::MSE(efpc_eval_list_BASIS[[2]],efpc_eval_list_BSPLINE[[2]])
MLmetrics::MSE(efpc_eval_list_BASIS[[2]],efpc_eval_list_GRID[[2]])

MLmetrics::MSE(efpc_eval_list_BASIS[[3]],efpc_eval_list_BSPLINE[[3]])
MLmetrics::MSE(efpc_eval_list_BASIS[[3]],efpc_eval_list_GRID[[3]])

# PLOT -------------------------------------------------------------------------

minn1 = min(min(sapply(efpc_eval_list_BASIS[[1]],min)), min(sapply(efpc_eval_list_GRID[[1]],min)), min(sapply(efpc_eval_list_BSPLINE[[1]],min)))-0.5
maxx1 = max(max(sapply(efpc_eval_list_BASIS[[1]],max)), max(sapply(efpc_eval_list_GRID[[1]],max)), max(sapply(efpc_eval_list_BSPLINE[[1]],max)))+0.5
minn2 = min(min(sapply(efpc_eval_list_BASIS[[2]],min)), min(sapply(efpc_eval_list_GRID[[2]],min)), min(sapply(efpc_eval_list_BSPLINE[[2]],min)))-0.5
maxx2 = max(max(sapply(efpc_eval_list_BASIS[[2]],max)), max(sapply(efpc_eval_list_GRID[[2]],max)), max(sapply(efpc_eval_list_BSPLINE[[2]],max)))+0.5
minn3 = min(min(sapply(efpc_eval_list_BASIS[[3]],min)), min(sapply(efpc_eval_list_GRID[[3]],min)), min(sapply(efpc_eval_list_BSPLINE[[3]],min)))-0.5
maxx3 = max(max(sapply(efpc_eval_list_BASIS[[3]],max)), max(sapply(efpc_eval_list_GRID[[3]],max)), max(sapply(efpc_eval_list_BSPLINE[[3]],max)))+0.5

line_neg = 0
res = 400

col1="firebrick1"
col2="dodgerblue3"
col3="darkslategray4"

#x11()
png(file = paste0("D:/Poli/TESI/Pics/FPCA/FPCA_basis_1.png"), width = 6000, height = 2000, units = "px", res = res)
par(mfrow=c(1,3))
persp(x=x1.grid, y=x2.grid, z=efpc_eval_list_BASIS[[1]], col=col1,
      xlab="",ylab="",zlab="",
      zlim=c(minn1,maxx1),
      ticktype='detailed')
title("1st FPC", line = line_neg)
persp(x=x1.grid, y=x2.grid, z=efpc_eval_list_BASIS[[2]], col=col1,
      xlab="",ylab="",zlab="",
      zlim=c(minn2,maxx2),
      ticktype='detailed')
title("2nd FPC", line = line_neg)
persp(x=x1.grid, y=x2.grid, z=efpc_eval_list_BASIS[[3]], col=col1,
      xlab="",ylab="",zlab="",
      zlim=c(minn3,maxx3),
      ticktype='detailed')
title("3rd FPC", line = line_neg)
dev.off()

#x11()
png(file = paste0("D:/Poli/TESI/Pics/FPCA/FPCA_grid_1.png"), width = 6000, height = 2000, units = "px", res = res)
par(mfrow=c(1,3))
persp(x=x1.grid, y=x2.grid, z=efpc_eval_list_GRID[[1]], col=col2,
      xlab="",ylab="",zlab="",
      zlim=c(minn1,maxx1),
      ticktype='detailed')
title("1st FPC", line = line_neg)
persp(x=x1.grid, y=x2.grid, z=efpc_eval_list_GRID[[2]], col=col2,
      xlab="",ylab="",zlab="",
      zlim=c(minn2,maxx2),
      ticktype='detailed')
title("2nd FPC", line = line_neg)
persp(x=x1.grid, y=x2.grid, z=efpc_eval_list_GRID[[3]], col=col2,
      xlab="",ylab="",zlab="",
      zlim=c(minn3,maxx3),
      ticktype='detailed')
title("3rd FPC", line = line_neg)
dev.off()

#x11()
png(file = paste0("D:/Poli/TESI/Pics/FPCA/FPCA_bspline_1.png"), width = 6000, height = 2000, units = "px", res = res)
par(mfrow=c(1,3))
persp(x=x1.grid, y=x2.grid, z=efpc_eval_list_BSPLINE[[1]], col=col3,
      xlab="",ylab="",zlab="",
      zlim=c(minn1,maxx1),
      ticktype='detailed')
title("1st FPC", line = line_neg)
persp(x=x1.grid, y=x2.grid, z=efpc_eval_list_BSPLINE[[2]], col=col3,
      xlab="",ylab="",zlab="",
      zlim=c(minn2,maxx2),
      ticktype='detailed')
title("2nd FPC", line = line_neg)
persp(x=x1.grid, y=x2.grid, z=efpc_eval_list_BSPLINE[[3]], col=col3,
      xlab="",ylab="",zlab="",
      zlim=c(minn3,maxx3),
      ticktype='detailed')
title("3rd FPC", line = line_neg)
dev.off()

