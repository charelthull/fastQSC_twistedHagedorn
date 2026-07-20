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

#include <iostream>

#include "linSolve.h"

using namespace cln;

void linsolvepE(cl_N * a,cl_N * b,cl_N x[],long int n)//decomposition with pivoting, then solving the equation
{
	long int d=1;
  	for(long int k=0;k<n-1;k++)
  	{
   		cl_N am = a[k*n+k]; //pivoting
	 	long int xm=k;
	 	for(long int i=k;i<n;i++)
	 	{
	   		if(abs(am)<abs(a[i*n+k]))
	  		{
	     		am=a[i*n+k];
	     		xm=i;
	   			cl_N c;
				for (long int i=0;i<n;i++)//switching lines
				{
					c=a[xm*n+i];
					a[xm*n+i]=a[k*n+i];
					a[k*n+i]=c;
	 			}
	 			c=b[xm];
	 			b[xm]=b[k];
	 			b[k]=c;
	   			d=-d;
	   		}
	 	}
     	for(long int i=k+1;i<n;i++)
     	{
	 		a[i*n+k]/=a[k*n+k];
	 		for(long int j=k+1;j<n;j++)
	 		{
	  			a[i*n+j]-=a[i*n+k]*a[k*n+j];
	 		}
     	}
  	}
	//inverting matrix
	cl_N y[n]; 
  	y[0]=b[0];
 	for(long int k=1;k<n;k++) //solving Ly=b
 	{
   		y[k]=b[k];
   		for(long int j=0;j<=k-1;j++)
   		{
     		y[k]-=a[k*n+j]*y[j];
   		}
   	}
   	x[n-1]=y[n-1]/a[(n-1)*n+n-1]; //solving Ux=y
   	for(long int k=n-2;k>=0;k--)
 	{
     	x[k]=complex(0,0);
     	for(long int j=k+1;j<n;j++)
     	{
			x[k]-=a[k*n+j]*x[j];
     	}
     	x[k]=(y[k]+x[k])/a[k*n+k];
   }
}

void linsolvepER(cl_R * a,cl_R * b,cl_R x[],long int n)//decomposition with pivoting, then solving the equation
{
	//long int d=1;
  	for(long int k=0;k<n-1;k++)
  	{
   		cl_R am = a[k*n+k]; //pivoting
	 	long int xm=k;
	 	for(long int i=k;i<n;i++)
	 	{
			if(k>n-10){std::cerr<<k<<" "<<i<<" "<<a[i*n+k]<<"\n";}
	   		if(abs(am)<abs(a[i*n+k]))
	  		{
	     		am=a[i*n+k];
	     		xm=i;
	   			cl_R c;
				for (long int i=0;i<n;i++)//switching lines
				{
					c=a[xm*n+i];
					a[xm*n+i]=a[k*n+i];
					a[k*n+i]=c;
	 			}
	 			c=b[xm];
	 			b[xm]=b[k];
	 			b[k]=c;
	   			//d=-d;
	   		}
	 	}
     	for(long int i=k+1;i<n;i++)
     	{
			//if(k>130){std::cerr<<k<<" "<<i<<"\n";};
	 		a[i*n+k]/=a[k*n+k];
	 		for(long int j=k+1;j<n;j++)
	 		{
	  			a[i*n+j]-=a[i*n+k]*a[k*n+j];
	 		}
     	}
  	}
	//inverting matrix
	cl_R y[n]; 
  	y[0]=b[0];
 	for(long int k=1;k<n;k++) //solving Ly=b
 	{
   		y[k]=b[k];
   		for(long int j=0;j<=k-1;j++)
   		{
     		y[k]-=a[k*n+j]*y[j];
   		}
   	}
   	x[n-1]=y[n-1]/a[(n-1)*n+n-1]; //solving Ux=y
   	for(long int k=n-2;k>=0;k--)
 	{
     	x[k]=0;
     	for(long int j=k+1;j<n;j++)
     	{
			x[k]-=a[k*n+j]*x[j];
     	}
     	x[k]=(y[k]+x[k])/a[k*n+k];
   }
}


void linsolveHouseholder(cl_R * a,cl_R * b,cl_R * x,int n){
	cl_R * v = new cl_R[n];
	cl_R nv;
	cl_R tmp;
	cl_R test;

	for(int i=0;i<n-1;i++){

		//Figure out the correct Householder transformation
		for(int j=i;j<n;j++){
			v[j]=a[j*n+i]/2;
		}

		nv=v[n-1]*v[n-1];
		for(int j=n-2;j>=i;j--){
			nv+=v[j]*v[j];
		}
		nv=sqrt(nv);

		if(v[i]>0 ){
			v[i]+=nv;
			a[i*n+i]=-2*nv;
			nv*=2*v[i];
		}else{
			v[i]-=nv;
			a[i*n+i]=2*nv;
			nv*=-2*v[i];
		}
		
		//Apply the Householder transformation to the lower right (n-i)x(n-i-1) block of a
		for(int k=i+1;k<n;k++){
			tmp=v[n-1]*a[(n-1)*n+k];
			for(int j=n-2;j>=i;j--){
				tmp+=v[j]*a[j*n+k];
			}

			tmp*=2/nv;

			for(int j=i;j<n;j++){
				a[j*n+k]-=v[j]*tmp;
			}
		}

		//Apply the Householder transformation to the lower (n-i) components of b
		tmp=v[n-1]*b[n-1];
		for(int j=n-2;j>=i;j--){
			tmp+=v[j]*b[j];
		}
		
		tmp*=2/nv;

		for(int j=i;j<n;j++){
			b[j]-=v[j]*tmp;
		}
	}

	x[n-1]=b[n-1]/a[(n-1)*n+n-1]; //solving Ux=y
   	for(int k=n-2;k>=0;k--)
 	{
     	x[k]=0;
     	for(int j=n-1;j>k;j--)
     	{
			x[k]-=a[k*n+j]*x[j];
     	}
     	x[k]=(b[k]+x[k])/a[k*n+k];
   }

   delete[] v;


}


void linsolveHouseholder_C(cl_N * a,cl_N * b,cl_N * x,int n){
	cl_N * v = new cl_N[n];
	cl_R nv;
	cl_N tmp;
	cl_N expphi;

	for(int i=0;i<n-1;i++){

		//Figure out the correct Householder transformation
		for(int j=i;j<n;j++){
			v[j]=a[j*n+i]/2;
		}
		nv=realpart(conjugate(v[i])*v[i]);
		for(int j=i+1;j<n;j++){
			nv+=realpart(conjugate(v[j])*v[j]);
		}
		nv=sqrt(nv);
		
		expphi=exp(complex(0,phase(a[i*n+i])));

		v[i]+=nv*expphi;
		a[i*n+i]=-2*nv*expphi;
		nv*=realpart(2*v[i]/expphi);
		
		//Apply the Householder transformation to the lower right (n-i)x(n-i-1) block of a
		for(int k=i+1;k<n;k++){
			tmp=conjugate(v[i])*a[i*n+k];
			for(int j=i+1;j<n;j++){
				tmp+=conjugate(v[j])*a[j*n+k];
			}
			
			tmp*=2/nv;

			for(int j=i;j<n;j++){
				a[j*n+k]-=v[j]*tmp;
			}
		}

		//Apply the Householder transformation to the lower (n-i) components of b
		tmp=conjugate(v[i])*b[i];
		for(int j=i+1;j<n;j++){
			tmp+=conjugate(v[j])*b[j];
		}
		tmp*=2/nv;
		for(int j=i;j<n;j++){
			b[j]-=v[j]*tmp;
		}
	}


	x[n-1]=b[n-1]/a[(n-1)*n+n-1]; //solving Rx=y
   	for(int k=n-2;k>=0;k--)
 	{
     	x[k]=a[k*n+n-1]*x[n-1];
     	for(int j=n-2;j>k;j--)
     	{
			x[k]+=a[k*n+j]*x[j];
     	}
     	x[k]=(b[k]-x[k])/a[k*n+k];
   }

   delete[] v;

}

void linOptimizeHouseholder(cl_R * a,cl_R * b,cl_R * x,int m,int n){
	/*
	The matrix a has size m x n, b has m components, x has n components*/
	
	cl_R * v = new cl_R[m];
	cl_R nv;
	cl_R tmp;
	cl_R test;

	for(int i=0;i<n-1;i++){

		//Figure out the correct Householder transformation
		for(int j=i;j<m;j++){
			v[j]=a[j*n+i]/2;
		}

		nv=v[m-1]*v[m-1];
		for(int j=m-2;j>=i;j--){
			nv+=v[j]*v[j];
		}
		nv=sqrt(nv);

		if(v[i]>0 ){
			v[i]+=nv;
			a[i*n+i]=-2*nv;
			nv*=2*v[i];
		}else{
			v[i]-=nv;
			a[i*n+i]=2*nv;
			nv*=-2*v[i];
		}
		
		//Apply the Householder transformation to the lower right (m-i)x(n-i-1) block of a
		for(int k=i+1;k<n;k++){
			tmp=v[m-1]*a[(m-1)*n+k];
			for(int j=m-2;j>=i;j--){
				tmp+=v[j]*a[j*n+k];
			}

			tmp*=2/nv;

			for(int j=i;j<m;j++){
				a[j*n+k]-=v[j]*tmp;
			}
		}

		//Apply the Householder transformation to the lower (m-i) components of b
		tmp=v[m-1]*b[m-1];
		for(int j=m-2;j>=i;j--){
			tmp+=v[j]*b[j];
		}
		
		tmp*=2/nv;

		for(int j=i;j<m;j++){
			b[j]-=v[j]*tmp;
		}
	}

	x[n-1]=b[n-1]/a[(n-1)*n+n-1]; //solving Rx=y: We forget about the additional components of b, they give the residue
   	for(int k=n-2;k>=0;k--)
 	{
     	x[k]=0;
     	for(int j=n-1;j>k;j--)
     	{
			x[k]-=a[k*n+j]*x[j];
     	}
     	x[k]=(b[k]+x[k])/a[k*n+k];
   }

   delete[] v;


}
