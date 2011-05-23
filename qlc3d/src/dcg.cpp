  /************************************************************************
 ***************************************************************************
*****************************************************************************
*****					 				*****
*****		    Preconditioned Conjugate Gradient			*****
*****		    Double Precision (double) version.			*****
*****					 				*****
*****************************************************************************
 ***************************************************************************
  ************************************************************************/

#include <math.h>
#include <vector>
#include "meshtest.h"
//static char DCGsid[] = "@(#)dcg.c	1.2  5/19/86";

//#define NULL 0
double ddot(int n, double *sx, int incx, double *sy, int incy )
/*
    PURPOSE
        Forms the dot product of a vector.

    INPUT
        n       Number of elements to sum.
        sx      Address of first element of x vector.
        incx    Incrament for the x vector.
        sy      Address of first element of y vector.
        incy    incrament for the y vector.

    OUPUT
        sdot    Dot product x and y.  Double returned
		due to `C' language features.
*/
{
  register int i;
  double stemp = 0.0e0;

  if( n<1 ) return( stemp );
  if( incx == incy ) {
    if( incx == 1 ) {
      /* Both increments = 1 */
      for( i=0; i<n; i++, sx++, sy++ )
	stemp += (*sx)*(*sy);
      return( stemp );
    }
    if( incx>0 ) {
      /* Equal, positive, non-unit increments. */
      for( i=0; i<n; i++, sx+=incx, sy+=incx)
	stemp += (*sx)*(*sy);
      return( stemp );
    }
  }
  /* Unequal or negative increments. */
  if( incx < 0 ) sx += ((-n+1)*incx + 1);
  if( incy < 0 ) sy += ((-n+1)*incy + 1);
  for( i=0; i<n; i++,sx+=incx,sy+=incy ) 
    stemp += (*sx)*(*sy);
  return( stemp );
}				/* End of ---DDOT--- */
void dexopy(int n,double *v, double *x,double *y, int itype )
/*
  Purpose: v = x op y
    To operate on the vectors x and y.

  Input:
    n		Number of elements to scale.
    x		First operand vector.
    y		Second operand vector.
    itype	Type of operation to perform:
		itype = 1 => '+'
		itype = 2 => '-'

  Output:
    v		Result vector of x op y.
*/
{
  register int i;

  if( n<1 ) return;

  if( itype == 1 )	/* ADDITION. */
    for( i=0; i<n; i++, v++, x++, y++ )
      *v = *x + *y;
  else			/* SUBTRACTION. */
    for( i=0; i<n; i++, v++, x++, y++ )
      *v = *x - *y;
}

void dcopy( int n, double *sx, int incx, double *sy, int incy )
/*
    PURPOSE
        Copies vector sx into vector sy.
 
    INPUT
        n    Number of components to copy.
	sx   Source vector
	incx Index increment for sx.
        incy Index increment for sy.
 
    OUTPUT
        sy   Destination vector.
*/
{
  register int i;

  if( n<1  ) return;
  if( incx == incy ) {
    if( incx == 1 ) {
      /* Both increments = 1 */
      for( i=0; i<n; i++ )
	*(sy++) = *(sx++);
      return;
    }
    if( incx > 0 ) {
      /* Equal, positive, non-unit increments. */
      for( i=0; i<n; i++, sx+=incx, sy+=incx)
	*sy = *sx;
      return;
    }
  }
  /* Non-equal or negative increments. */
  if( incx < 0 ) sx += ((-n+1)*incx + 1);
  if( incy < 0 ) sy += ((-n+1)*incy + 1);
  for( i=0; i<n; i++,sx+=incx,sy+=incy ) 
    (*sx) = (*sy);
  return;
}// end void dcopy

void daxpyx(int n, double sa, double *sx,int incx, double *sy, int incy )
/*
  PURPOSE
    Vector times a scalar plus a vector.  sx = sy + sa*sx.

  INPUT
    n		Number of elements to multiply.
    sa		Scalar to multiply by.
    sx		Pointer to double vector to scale.
    incx	Storage incrament for sx.
    sy		Pointer to double vector to add.
    incy	Storage incrament for sy.

  OUTPUT
    sx		sx = sy + sa*sx
*/
{
  register int i;

  if( n<=0 || sa==0.0 ) return;
  if( incx == incy ) {
    if( incx == 1 ) {
      /* Both increments = 1 */
      for( i=0; i<n; i++, sx++, sy++ )
	*sx = *sy + sa*(*sx);
      return;
    }
    if( incx>0 ) {
      /* Equal, positive, non-unit increments. */
      for( i=0; i<n; i++, sx+=incx, sy+=incx)
	*sx = *sy + sa*(*sx);
      return;
    }
  }
  /* Unequal or negative increments. */
  if( incx < 0 ) sx += ((-n+1)*incx + 1);
  if( incy < 0 ) sy += ((-n+1)*incy + 1);
  for( i=0; i<n; i++,sx+=incx,sy+=incy ) 
    *sx = *sy + sa*(*sx);
}// end void daxpyx
void daxpy(int n, double sa, double *sx,int incx,double *sy,int incy )
/*
  PURPOSE
    Vector times a scalar plus a vector.  sy = sy + sa*sx.

  INPUT
    n		Number of elements to multiply.
    sa		Scalar to multiply by.
    sx		Pointer to double vector to scale.
    incx	Storage incrament for sx.
    sy		Pointer to double vector to add.
    incy	Storage incrament for sy.

  OUTPUT
    sy		sy = sy + sa*sx
*/
{
  register int i;

  if( n<=0 || sa==0.0 ) return;
  if( incx == incy ) {
    if( incx == 1 ) {
      /* Both increments = 1 */
      for( i=0; i<n; i++,sy++,sx++ )
	*sy += sa*(*sx);
      return;
    }
    if( incx>0 ) {
      /* Equal, positive, non-unit increments. */
      for( i=0; i<n; i++,sx+=incx,sy+=incx )
	*sy += sa*(*sx);
      return;
    }
  }
  /* Unequal or negative increments. */
  if( incx < 0 ) sx += ((-n+1)*incx + 1);
  if( incy < 0 ) sy += ((-n+1)*incy + 1);
  for( i=0; i<n; i++,sx+=incx,sy+=incy ) 
    *sy += sa*(*sx);
}// end void daxpy
 

//msolve( n, mptr, aptr, z, r, atimes )
int msolve(SparseMatrix *A, double *z, double *r)
{ // this function performs the preconditioninf using the Jacobi - diagonal method
 //i             nt msolve( n, mptr, z, r, atimes )
 //	        The user should write a routine that solves Mz=r for z, given
 //		a vector r.  Return non-zero if an error occures.
 // Here M is the diagonal of A
 
	if (A->rows != A->cols) {printf("dcg.cpp, msolve - error, A->rows != A->cols\n"); return(1);}
	int i;
	
	// z[i] = r[i] / A[i,i]
	for (i=0;i<A->rows;i++)
		z[i] = r[i] / A->sparse_get(i,i); // 
return (0);
}//end int msolve 
 
int dcg(SparseMatrix *aptr, double *x, double *b, int itmax, double eps, double *err )
/*
  Purpose:
    This routine does preconditioned conjugate gradient iteration
    on the symmetric positive definte system Ax = b.

  Input:
    n		Size of system to be iterated (ie A is nxn).
    aptr	Pointer to the structure defining the matrix A.
		This structure need not be known by this routine.
		The structure is passed on to the atimes routine which does
		the matrix multiply.  See below for the defintion of atimes.
    mptr	Pointer to the structure defining the preconditioning matrix 
	        M.  M should be chosen so that:
		  Mz = r is easy to solve for z, given r (compared to Ax = b).
		  Inv(M) is close to Inv(A).  Inv = inverse.
		Obviously, the better choice of M the faster the convergence.
    x		Initial guess of the solution.  If no information on the 
	        solution is available then use a zero vector or b.
    b		Right hand side.
    eps		Requested error tollerence.  System is iterated until
		||b-Ax||/||b|| < eps.  Normal choice is 1.0e-6.
    itmax	Maximum number of iterations the user is willing to allow.

  Output:
    dcg		Returns number of iterations taken.  Error Returns follow:
    		dcg = itmax 	May have not converged.  Check err.
		dcg = -1	System size too small (n<2).
		dcg = -2	Could not get memory for tempory vectors.
		dcg = -3	Error tollerence requested too tight.
		dcg = -4	Divergence of iteration detected.
		dcg = -10	User declared error in atimes.
		dcg = -11	User declared error in msolve.
    err		L2 norm of the relative residual error estimate 
	        ||b-Ax||/||b||.  This estimate of the true error is not always
		very accurate.  Note that err is a return value so the user
		must call dcg with a pointer to err.

  User Defined Functions and Structures.
    aptr	The user defines an external structure DCGAMAT and hands this
		routine a pointer to that structure.  The structure should
		contain all the information the user needs to calculate
		A times a vector.  See DCG.H for an example.

    mptr	The user defines an external structure DCGMMAT and hands this
		routine a pointer to that structure.  The structure should
		contain all the information the user needs to solve the system
		Mz=r for z.  See DCG.H for an example.

    int atimes( n, aptr, x, y )  
	        The user should write a routine that calculates y = Ax, given
		an x vector.  Return non-zero if an error occures.
    int msolve( n, mptr, z, r, atimes )
	        The user should write a routine that solves Mz=r for z, given
		a vector r.  Return non-zero if an error occures.
*/
{
  //char   *malloc();
  
  int n = aptr->rows;
  register int i;
  double  *r, *z, *p;
  //extern double ddot();
  double bnorm, bknum, bkden, bk;
  double akden, ak;
  double toobig;	/* When err est gets this big we signal divergence! */

  /* Check input data and allocate temporary storage. */
  if( n<2 ) return( -1 );
  if( (r = (double *)malloc( n*sizeof(double) )) == NULL ) return( -2 );
  if( (z = (double *)malloc( n*sizeof(double) )) == NULL ) return( -2 );
  if( (p = (double *)malloc( n*sizeof(double) )) == NULL ) return( -2 );
  if( eps < 1.0e-12 ) return( -3 );

  /* Calculate norm of right hand size. */
  bnorm = sqrt( ddot( n, b, 1, b, 1 ) );

  /* Calculate r = b - Ax and initinal error. */
  //if( atimes( n, aptr, x, r ) ) return( -10 );
  //int atimes( n, aptr, x, y )  
  // The user should write a routine that calculates y = Ax, given
  //an x vector.  Return non-zero if an error occures.
  aptr->MultiplyArray(x,r); //   r = Ax 
  dexopy( n, r, b, r, 2 );  //  r = b- r  , 2 is minus , 1 is plus 
	// now r = b-Ax
	
  *err = sqrt( ddot(n, r, 1, r, 1 ) )/bnorm;
  if( *err < eps ) return( 0 );
  toobig = *err * 1.0e2;

  /* Iterate!!! */
  for( i=0; i<itmax; i++ ) {

    /* Solve Mz = r. */
    //if( msolve( n, mptr, aptr, z, r, atimes ) ) return( -11 );
	if( msolve( aptr, z, r) ) return( -11 );  /// PRECONDITIONING 
    /* Calculate bknum = (z, Mz) and p = z (first iteration). */
    if( i == 0 ) {
      dcopy( n, z, 1, p, 1 );
      bknum = bkden = ddot( n, z, 1, r, 1 );
    } 
    else {
      /* Calculate bknum = (z, r), bkden and bk. */
      bknum = ddot( n, z, 1, r, 1 );
      bk    = bknum/bkden;
      bkden = bknum;

      /* Calculate p = z + bk*p, z = Ap, akden = (p, Ap) and ak. */
      daxpyx( n, bk, p, 1, z, 1 );
    }
   // if( atimes( n, aptr, p, z ) ) return( -10 );
    aptr->MultiplyArray(p,z);
	akden = ddot( n, p, 1, z, 1 );
    ak    = bknum/akden;

    /* Update x and r. Calculate error. */
    daxpy( n,  ak, p, 1, x, 1 );
    daxpy( n, -ak, z, 1, r, 1 );
    *err = sqrt( ddot( n, r, 1, r, 1 ) )/bnorm;
    if( *err < eps ) return( i+1 );
    if( *err > toobig ) return( -4 );
  }				/* end for loop. */
  return( itmax );
}				/* end dcg. */
