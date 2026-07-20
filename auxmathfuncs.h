#include <cln/cln.h>
#include <cln/number.h>
#include <cln/integer.h>
#include <cln/real.h>
#include <cln/float.h>
#include <cln/io.h>
#include <cln/float_io.h>
#include <cln/malloc.h>

using namespace cln;


cl_R shifted_factorial (int n,int s); // n(n-1)...(n-s+1)
cl_R shifted_factorial (cl_RA n,int s); // n(n-1)...(n-s+1)

cl_R binomial ( int n, int k); // n(n-1)...(n-s+1)/(k!)
cl_R binomial ( cl_RA n, int k); // n(n-1)...(n-s+1)/(k!)

cl_R power(cl_I  i, int k);//To compute powers of 

cl_R d(int n, int l, const cl_R * g); //An auxiliary function for expanding the product of two \bP-functions at large u.

int max (int a, int b);
int min (int a, int b);
