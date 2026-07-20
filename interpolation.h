#include <cln/cln.h>
#include <cln/number.h>
#include <cln/integer.h>
#include <cln/real.h>
#include <cln/float.h>
#include <cln/io.h>
#include <cln/float_io.h>
#include <cln/malloc.h>

using namespace cln;

cl_R interpol_n (const cl_R * y,const cl_R * x, cl_R  target,float_format_t precision, int order, int n);
/*
 * x: Points at which we are given values
 * y: An array of pointers to the values at the points
 * target: we return the value of the interpolation polynomial at this point.
 * order: The order of the interpolation polynomial
 * n: The number of points at which we are given values
 */

cl_R interpol_n (const cl_R * y,const cl_R * x, cl_R  target,float_format_t precision, int n);
/*
 * x: Points at which we are given values
 * y: An array of pointers to the values at the points
 * target: we return the value of the interpolation polynomial at this point.
 * n: The number of points at which we are given values. We interpolate at order n-1
 */

cl_R interpol_2 (const cl_R * y[],const cl_R * x, cl_R  target,float_format_t precision);
/*
 * Second order polynomial interpolation
 *
 * x: Points at which we are given values
 * y: An array of pointers to the values at the points
 * target: we return the value of the interpolation polynomial at this point.
 *
 */

cl_R interpol_1 (const cl_R * y[],const cl_R * x, cl_R  target);
/*
 * Linear interpolation
 * 
 * x: Points at which we are given values
 * y: An array of pointers to the values at the points
 *
 */