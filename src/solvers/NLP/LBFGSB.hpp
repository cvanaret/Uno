#ifndef LBFGSB_H
#define LBFGSB_H

#include <iostream>
#include <vector>
#include <map>

Cextern void setulb_ ANSI((int *n, int *m, double *x, double *l, double *u, int *nbd, double *f, double *g,
			   double *factr, double *pgtol, double *wa, int *iwa, char *task, int *iprint,
			   char *csave, int *lsave, int *isave, double *dsave, ftnlen cstype_len));

class lBFGSb {

public:

  //! ... class constructor allocates memory
  void lBFGSb(int n){
    new x[n], l[n], u[n], g[n], wa[1];
    new nbd[n], iwa[1];
  };

  //! ... class destructor de-allocates memory
  void ~lBFGSb(){
    delete x, l, u, g, wa;
    new nbd, iwa;
  };

  //! ... solver (copy what's below???)    *** TODO ***
  int solve(){
    return 0; //! ... return error code (LATER)
  };
  
private:
  
  //! ... set up storage needed by L-BFGS-B solver
  int n = 10, m = 5, iprint = 1;
  double factr  = 1.0d+7, pgtol  = 1.0d-5;
  char task[4], csave[60];
  logical lsave[4];
  int isave[44];
  double f, RedGrad;
  double dsave[29];
  int  *nbd, *iwa;
  double *x, *l, *u, *g, *wa;

  //! ... allocate/de-allocate storage???   *** TODO ***
  
  //! ... copy/define lower/upper bounds & set nbd(i) = {0|1|2|3} = {unbdd|LbdOnly|BothBnds|UbdOnly}
  for ( int i=0; i<n; i++ ) {
    nbd[i] = 2;
    l[i]   = (-1)**2*0.2d0*i; 
    u[i]   =  0.8d0*i;
  }; // end for

  //! ... copy/define the starting point
  for ( int i=0; i<n; i++ ) {
    x[i] = 1.0d0;
  }; // end for 

  //! ... we start the iteration by initializing task.
  task = 'START';

  //! ... the beginning of the optimization loop
  while ( (strcmp(task[0:1], "FG") == 0) || (strcmp(task, "NEW_X") == 0) || (strcmp(task, "START") == 0) ){

    // ... call the L-BFGS-B code.
    setulb_( &n, &m, x, l, u, nbd, f, g, factr, pgtol, 
	     wa, iwa, task, iprint, csave, lsave, isave, dsave );
    
    if (task[0:1]=="FG")
      {
	//! ... compute gradient g for the sample problem.
	f=0.0d0;
	for ( int i=0; i<n; i++ ) {
	  f = f + 100*(x[i+1]-x[i]**2)**2 + (1 - x[i])**2;
	}; // end for
	
	//! ... compute gradient g for the sample problem.
	g[0] = - 400*(x[1] - x[0]**2)*x[0] - 2*(1-x[0]);
	for ( int i=1; i<n-1; i++ ) {
	  g[i] = - 400*(x[i+1] - x[i]**2)*x[i] - 2*(1-x[i]) + 200*(x[i] - x[i-1]**2);
	}; // end for
	g[n] = 200*(x[n] - x[n-1]**2);
	  
      }; // end if
    
  }; // end while for L-BFGS-B loop 

  //! ... print output of this example
  printf ("Final L-BFGS-B Solution\n");
  printf ("lower bound   x-value      upper bound  gradient\n");
  RedGrad = 0.0d0;
  for ( int i=0; i<n-1; i++ ) {
    printf(6," %E   %E  %E  %E\n", l[i],x[i],u[i],g[i]);
    RedGrad = RedGrad + abs( min( x[i]-l[i] , u[i]-x[i] )*g[i] );
  }; // end for
  printf ("Reduced Gradient Norm = %E\n", RedGrad);
  
}; // lBFGSb

#endif // LBFGSB_H
