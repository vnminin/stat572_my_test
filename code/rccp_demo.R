## An example function that is hard to vectorize in R

aux.mat = function(my.vec, my.scalar){

  space.size = length(my.vec)
  
  int.matrix = matrix(0,nrow = space.size, ncol = space.size)
  
  for (i in c(1:space.size)){
    for (j in c(1:space.size)){
      if (abs(my.vec[i] - my.vec[j]) < 10^-7){
        int.matrix[i,j] = exp(my.vec[j] * my.scalar)*my.scalar
      }else{
        int.matrix[i,j] = (exp(my.vec[i] * my.scalar) -
                           exp(my.vec[j] * my.scalar))/
                             (my.vec[i] - my.vec[j])
      }
    }
  }

  return(int.matrix)
}



## A classical way of calling C/C++ using .Call()
## Learn more by (re)reading "Writing R Extentions"

## 1. Show the source of "aux_mat.c"
## 2. Create a shared library: > R CMD SHLIB aux_mat.c

## 3. Load the shared library into R
dyn.load("aux_mat.so")

## 4. Call C function "aut_mat"
test.vec = runif(5)
test.scalar = runif(1)

test.output1 = aux.mat(test.vec, test.scalar)
test.output2 = .Call("aux_mat", test.vec, test.scalar)

## Now let's test the running times
system.time(for (i in 1:1000) test.outpu1=aux.mat(test.vec, test.scalar))
system.time(for (i in 1:1000) test.outpu2=.Call("aux_mat", test.vec, test.scalar))

## Now let's accomplish the same task with Rcpp and inline
## Note: inline requires Rcpp so the latter will be loaded automatically

library(Rcpp)
library(inline)

# and define our version in C++
rcpp.src <- '
  #include <cmath>
  
  NumericVector xx(x);
  double yy = as<double>(y);
  int n_xx = xx.size();
  NumericMatrix ans(n_xx, n_xx);

  for(int i = 0; i < n_xx; i++) {
    for(int j = 0; j < n_xx; j++){
      if ((xx[i]-xx[j]<10e-7) && (xx[i]-xx[j]>-10e-7)){
        ans(i,j) = exp(xx[j]*yy)*yy;
      }else{
	      ans(i,j) = (exp(xx[i]*yy) - exp(xx[j]*yy))/(xx[i] - xx[j]);
      }
    }
  }
  
  return wrap(ans);'
 
aux.mat.rcpp = cxxfunction(signature(x="numeric",
y="numeric"), body=rcpp.src, plugin="Rcpp")

aux.mat.rcpp(test.vec, test.scalar)

## Another running time comparison

system.time(for (i in 1:10000) test.output1=aux.mat(test.vec, test.scalar))
system.time(for (i in 1:10000) test.output2=.Call("aux_mat", test.vec, test.scalar))
system.time(for (i in 1:10000) test.output3=aux.mat.rcpp(test.vec, test.scalar))

## OK, now suppose we want to use some linear algebra within Rcpp code
## There are three options: RcppArmadillo, RcppEigen, and RcppGSL. 
## We will use the former.
## Learn more about armadillo at http://arma.sourceforge.net/

## Supppose we want to exponentiate a matrix A given its eigen decomposition A = UDV

library(RcppArmadillo)

## Let's create an example matrix and compute its eigen decomposition
A = matrix(runif(100),100,100)
## Make it symmetric for stability of eigen decomposition
A = A + t(A)

A.eigen = eigen(A)

## Here is an R function to compute the matrix exponential

matexp.eigen = function(eigen.list){
   return(eigen.list$vectors%*%
    diag(exp(eigen.list$values))%*%solve(eigen.list$vectors))
}

test.matexp1 = matexp.eigen(A.eigen)

## Now in Rccp using RccpArmadillo
## Note: You will need to install GNU Fortran (available from CRAN; look under tools) 
matexp.code  = '
Rcpp::List reigen(eigenlist);  

arma::colvec eigen_values	= Rcpp::as<arma::colvec>(reigen["values"]); 
arma::mat eigen_vectors	= Rcpp::as<arma::mat>(reigen["vectors"]);

return Rcpp::wrap(eigen_vectors*(arma::diagmat(arma::exp(eigen_values)))*(arma::inv(eigen_vectors)));'

matexp.eigen.rcpp = cxxfunction(signature(eigenlist="list"), body=matexp.code,plugin="RcppArmadillo")

test.matexp2 = matexp.eigen.rcpp(A.eigen)

system.time(for (i in 1:10000) test.matexp1=matexp.eigen(A.eigen))
system.time(for (i in 1:10000) test.matexp2=matexp.eigen.rcpp(A.eigen))

## 3-D arrays in R = cubes in Armadillo

src <- '
  using namespace Rcpp;

NumericVector vecArray(myArray);
IntegerVector arrayDims(myDims);

arma::cube cubeArray(vecArray.begin(), arrayDims[0], arrayDims[1], arrayDims[2], false);

//change one element in the array/cube
cubeArray(0,0,0) = 518;  

return(wrap(cubeArray));  
'

readCube = cxxfunction(signature(myArray="numeric", myDims="integer"),body=src, plugin="RcppArmadillo")

testArray = array(rnorm(18), dim=c(3,3,2))
print(testArray)
readCube(testArray,dim(testArray))
