#include <cln/cln.h>
#include <cln/number.h>
#include <cln/integer.h>
#include <cln/real.h>
#include <cln/float.h>
#include <cln/io.h>
#include <cln/float_io.h>
#include <cln/malloc.h>

#include "linSolve.h"

using namespace cln;

cl_R interpol_n (const cl_R * y,const cl_R * x, cl_R  target,float_format_t precision, int order, int n){
    /*
    order: The order of the interpolation polynomial,
    n: The number of points we use for the interpolation
    */
    if(n<order+1){
        std::cerr<<"Not enough data provided for interpolation! "<<order<<" "<<n;
        exit(2);
    }
    cl_R * m = new cl_R[n*(order+1)];
    cl_R * v = new cl_R[n];
    cl_R * poly=new cl_R[order+1];
    cl_R tmp;

    for( int i=0;i<n;i++){
        m[(order+1)*i]=cl_float(1,precision);
        for(int j=1;j<order+1;j++){
            m[(order+1)*i+j]=expt(x[i],j);
        }
        v[i]=y[i];
    }
    linOptimizeHouseholder(m,v,poly,n,order+1);
    tmp=poly[0];
    for(int i=1;i<order+1;i++){
        tmp+=expt(target,i)*poly[i];
    }
    //tmp=cl_float(tmp,precision);

    delete[] m;
    delete[] v;
    delete[] poly;

    return tmp;
}

cl_R interpol_n (const cl_R * y,const cl_R * x, cl_R  target,float_format_t precision, int n){
    cl_R * m = new cl_R[(n+1)*(n+1)];
    cl_R * v = new cl_R[n+1];
    cl_R * poly=new cl_R[n+1];
    cl_R tmp;

    for( int i=0;i<n+1;i++){
        m[(n+1)*i]=cl_float(1,precision);
        for(int j=1;j<n+1;j++){
            m[(n+1)*i+j]=expt(x[i],j);
        }
        v[i]=y[i];
    }
    linsolveHouseholder(m,v,poly,n+1);
    tmp=poly[0];
    for(int i=1;i<n+1;i++){
        tmp+=expt(target,i)*poly[i];
    }

    delete[] m;
    delete[] v;
    delete[] poly;

    return tmp;
}

cl_R interpol_2 (const cl_R * y[], const cl_R * x, cl_R target, float_format_t precision){
    cl_R * m = new cl_R[3*3];
    cl_R * v = new cl_R[3];
    cl_R * poly=new cl_R[3];
    cl_R tmp;

    for( int i=0;i<3;i++){
        m[3*i]=cl_float(1,precision);
        for(int j=1;j<3;j++){
            m[3*i+j]=expt(x[i],j);
        }
        v[i]=*y[i];
    }
    linsolveHouseholder(m,v,poly,3);
    tmp=poly[0];
    for(int i=1;i<3;i++){
        tmp+=expt(target,i)*poly[i];
    }

    delete[] m;
    delete[] v;
    delete[] poly;

    return tmp;
}

cl_R interpol_1 (const cl_R * y[],const cl_R * x, cl_R target){

    return ((*y[0])-(*y[1]))/(x[0]-x[1])*(target)+((x[0])*(*y[1])-(x[1])*(*y[0]))/((x[0])-(x[1]));

}