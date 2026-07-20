#include <cln/cln.h>
#include <cln/number.h>
#include <cln/integer.h>
#include <cln/real.h>
#include <cln/float.h>
#include <cln/io.h>
#include <cln/float_io.h>
#include <cln/malloc.h>

#include <cmath>
#include <iostream>

using namespace cln;

cl_R shifted_factorial(int n,int s){
    cl_I tmp=1;
    for(cl_I k=0;k<s;k++){
        tmp*=(n-k);
    };
    return tmp;
};

cl_R shifted_factorial(cl_RA n,int s){
    cl_RA tmp=1;
    for(cl_RA k=0;k<s;k++){
        tmp*=(n-k);
    };
    return tmp;
};

cl_R binomial(int n, int k){
    return shifted_factorial(n,k)/shifted_factorial(k,k);
};

cl_R binomial(cl_RA n, int k){
    return shifted_factorial(n,k)/shifted_factorial(k,k);
};

cl_R power(cl_I  i, int k){

    if(i==-1){
        if(k%2==0){
            return 1;
        }
        else{
            return -1;
        }
    }
    else{
        return expt(i,k);
    }

};


cl_R d(int n, int l,const cl_R * g){
    cl_R tmp=0;
    cl_R tmp2;

    for(int m=0;m<=l;m++){
        tmp2=0;
        for(int  k=0;k<=m;k++){
            tmp2+=power(-1,m-k)*binomial(m,k)*binomial(k*recip(cl_I(2)),l);
        }
        tmp+=power(2,-m)*binomial(n,m)*tmp2;
    }
    tmp*=expt(2*(*g),2*l);

    return tmp;
};

int max(int x,int y){
    if(x<y){
        return y;
    }
    else{
        return x;
    }
}

int min(int x,int y){
    if(x>y){
        return y;
    }
    else{
        return x;
    }
}