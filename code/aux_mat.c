#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

SEXP aux_mat(SEXP x, SEXP y){
  int i, j, nx;
  double *vec, scalar;
  SEXP ans;
  
  nx = length(x);
  x = coerceVector(x, REALSXP); 
  y = coerceVector(y, REALSXP); 
  
  vec = REAL(x);
  scalar = REAL(y)[0];
    

  PROTECT(ans = allocMatrix(REALSXP, nx, nx));

  for(i = 0; i < nx; i++) {
    for(j = 0; j < nx; j++){
      if ((vec[i]-vec[j]<10e-7) && (vec[i]-vec[j]>-10e-7)){
	REAL(ans)[i + nx*j] = exp(vec[j]*scalar)*scalar;
      }else{
	REAL(ans)[i + nx*j] = (exp(vec[i]*scalar) - exp(vec[j]*scalar))/(vec[i] - vec[j]);
      }
    }
  }
  
  UNPROTECT(1);
  return(ans);
}
