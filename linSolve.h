/* Library to solve linear problems valued in cln types. */

#include <math.h>
#include <algorithm>
#include <cln/cln.h>
#include <cln/number.h>
#include <cln/integer.h>
#include <cln/real.h>
#include <cln/float.h>
#include <cln/io.h>
#include <cln/float_io.h>
#include <cln/malloc.h>

using namespace cln;

void linsolvepE(cl_N * a,cl_N * b,cl_N x[],long int n);
/*
 * Function to solve a linear problem A x=b for complex valued matrices
 * This modifies the input matrices A and b
 * The memory for the result is preallocated and passed to the function as x
 * */

void linsolvepER(cl_R * a,cl_R * b,cl_R x[],long int n);
/*
 * Function to solve a linear problem A x=b for real valued matrices.

 * Uses Gauss elimination with pivoting to compute an LU-decomposition

 * This modifies the input matrices A and b
 * The memory for the result is preallocated and passed to the function as x
 * */

void linsolveHouseholder(cl_R * a,cl_R * b,cl_R * x,int n);
/*
 * Function to solve a linear problem A x=b for real valued matrices

 * Uses Householder reflections to compute a QR-decomposition

 * This modifies the input matrices A and b
 * The memory for the result is preallocated and passed to the function as x
 * */

void linsolveHouseholder_C(cl_N * a,cl_N * b,cl_N * x,int n);
/*
 * Function to solve a linear problem A x=b for complex valued matrices

 * Uses Householder reflections to compute a QR-decomposition
 
 * This modifies the input matrices A and b
 * The memory for the result is preallocated and passed to the function as x
 * */

void linOptimizeHouseholder(cl_R * a,cl_R * b,cl_R * x,int m,int n);
/*
 The matrix a has size m x n, b has m components, x has n components with 
 */
