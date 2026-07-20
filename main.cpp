#include <cln/cln.h>
#include <cln/number.h>
#include <cln/integer.h>
#include <cln/real.h>
#include <cln/float.h>
#include <cln/io.h>
#include <cln/float_io.h>
#include <cln/malloc.h>

#include <cmath>
#include <ctime>
#include <filesystem>
#include <cstdlib> 
#include <iostream>
#include <string>

#include <fstream>
using std::ofstream;
using std::ifstream;

#include "QSC.h"
#include "interpolation.h"
#include "auxmathfuncs.h"

using namespace cln;

void writer(QSC_R * qq,ofstream & c,  ofstream & d){

        c<<qq->get_parameters()[0]<<" ";

        for(int i=0;i<qq->get_n_params();i++){
            c<<*(qq->get_fitCoeff(i))<<" ";
        }
        c<<std::endl;

        d<<double_approx(qq->get_parameters()[0])<<" "<<double_approx(*(qq->get_fitParameter()))<<std::endl;
        
}


void read_previous(int i, int n_omega, int n_terms, int n, int n_old, QSC_R * q_arr[],int precision, std::string loc ,ofstream * o,ofstream * outdata, ofstream * outdouble){
    ifstream input;
    ifstream input2;
    cl_R dummy;
    input.open(loc+"rangeSqrtG_omega"+std::to_string(double_approx(sin(pi(float_format(20))*(i+1)/(2*(n_omega+1)))))+".txt");
    input2.open(loc+"rangeSqrtG_omega"+std::to_string(double_approx(sin(pi(float_format(20))*(i+1)/(2*(n_omega+1)))))+"_onlyT_double.txt");
    cl_R * params_init=new cl_R[2*n_terms+1];
    cl_R params[2];

    //Take care of the data that we will not be using
    for(int k=0;k<n_old-n;k++){
        for(int l=0;l<2*(n_terms+1);l++){
            input>>dummy;
            *outdata<<dummy<<" ";
        }
        *outdata<<std::endl;
        input2>>dummy;
        *outdouble<<dummy<<" ";
        input2>>dummy;
        *outdouble<<dummy<<std::endl;
    }

    //Read the data that we will use to continue computations
    for(int k=0;k<min(n_old,n);k++){

        input>>params[0];
        *outdata<<params[0]<<" ";
        params[0]=cl_float(params[0],float_format(precision));
        params[1]=sin(pi(float_format(precision))*(i+1)/(2*(n_omega+1)));

        for(int l=0;l<2*n_terms+1;l++){
            input>>params_init[l];
            *outdata<<params_init[l]<<" ";
        }
        *outdata<<std::endl;

        //Only partially initialize the QSC instance: As it has already been optimized in a previous run we only need the QSC_R object to hold the data conveniently for later use
        q_arr[k]=new QSC_R(params_init,params,n_terms,precision,o,false);   

        input2>>dummy;
        *outdouble<<dummy<<" ";
        input2>>dummy;
        *outdouble<<dummy<<std::endl;

    }

    input.close();
    input2.close();
    delete[] params_init;

}




int main(int argc, char * argv[]){



    //Some parameters for the program
    int n_terms=16; //Number of terms we keep in Qai, P's have 1 more term
    int precision=160;
    int max_iter=25;
    int ytol=-18;
    int err_tol=-10;
    int max_interpol_order=7;

    int n_g = 20;   //We will compute for n_g values of the coupling
    cl_R * params = new cl_R[2*(n_g)];

    for(int i=0;i<n_g;i++){
        params[2*i]=sqrt(cl_float((i+1)*recip(cl_I(500)),float_format(precision))); //The values of coupling for which the Hagedorn temperature will get computed
        params[2*i+1]=sin(pi(float_format(precision))*(9+1)/(2*(19+1)));     //Set the value of the chemical potential
    }

    cl_R * params_init = new cl_R[2*n_terms+1];     //To store initial data we feed to the optimization algorithm

    cl_R * g2= new cl_R[max_interpol_order+1];      //Some temporary data for extrapolation to get initial data for the optimization
    cl_R tmp_param[max_interpol_order+1];


    QSC_R * q_arr[max_interpol_order+2];    //

    //Opening the files for data output and diagnostics.

    ofstream outstream;

    ofstream outcoeffs;

    ofstream outdouble;

    if(std::filesystem::exists("rangeSqrtG_omega"+std::to_string(double_approx(sin(pi(float_format(20))*(9+1)/(2*(19+1)))))+".txt")){
        time_t timestamp = time(NULL);
        std::filesystem::copy("rangeSqrtG_omega"+std::to_string(double_approx(sin(pi(float_format(20))*(9+1)/(2*(19+1)))))+".txt","rangeSqrtG_omega"+std::to_string(double_approx(sin(pi(float_format(20))*(9+1)/(2*(19+1)))))+"_backup"+ctime(&timestamp) +".txt");
    }

    outcoeffs.open("rangeSqrtG_omega"+std::to_string(double_approx(sin(pi(float_format(20))*(9+1)/(2*(19+1)))))+".txt");

    if(! outcoeffs){
        std::cerr<<"Error loading file.";
        outcoeffs.close();
        exit(1);
    }

    if(std::filesystem::exists("rangeSqrtG_omega"+std::to_string(double_approx(sin(pi(float_format(20))*(9+1)/(2*(19+1)))))+"_onlyT_double.txt")){
        time_t timestamp = time(NULL);
        std::filesystem::copy("rangeSqrtG_omega"+std::to_string(double_approx(sin(pi(float_format(20))*(9+1)/(2*(19+1)))))+"_onlyT_double.txt","rangeSqrtG_omega"+std::to_string(double_approx(sin(pi(float_format(20))*(9+1)/(2*(19+1)))))+"_onlyT_double_backup"+ctime(&timestamp) +".txt");
    }
            
    outdouble.open("rangeSqrtG_omega"+std::to_string(double_approx(sin(pi(float_format(20))*(9+1)/(2*(19+1)))))+"_onlyT_double.txt");

    if(! outdouble){
        std::cerr<<"Error loading file.";
        outdouble.close();
        exit(1);
    }

    outstream.open("messages_omega"+std::to_string(double_approx(sin(pi(float_format(20))*(9+1)/(2*(19+1)))))+".txt", std::ofstream::out | std::ofstream::app);
    
    if(! outstream){
        std::cerr<<"Error loading file.";
        outcoeffs.close();
        exit(1);
    }

    outstream<<"\nOmega="<<double_approx(sin(pi(float_format(20))*(9+1)/(2*(19+1))))<<std::endl;

    //If an optional argument is passed to the program, old data is loaded.
    int j=0;
    if(argc>1){
        j=std::stoi(argv[1]);
        read_previous(9,19,n_terms,max_interpol_order+1,j,q_arr,precision,"./Old/",&outstream,&outcoeffs,&outdouble);
    }

    //We treat the first value of coupling independently as we have to treat the case where we do not yet have data to interpolate from.
    outstream<<"\n\nPicking things back up at g="<<double_approx(params[2*j])<<"\n\n"<<std::endl;

    int index=min(j,max_interpol_order+1);

    if(index!=0){

        for(int l=0;l<index;l++){
            g2[l]=expt(q_arr[index-l-1]->get_parameters()[0],2);    //Grabbing the coupling values for the data we extrapolate
        }
        for(int k=0;k<2*n_terms+1;k++){
            for(int l=0;l<index;l++){
                tmp_param[l]= (*q_arr[index-l-1]->get_fitCoeff(k)); //Grabbing the parameter values we extrapolate
            }
            params_init[k]=interpol_n(tmp_param,g2,expt(params[2*j],2),float_format(precision),min(index-1,max_interpol_order),index);  //extrapolating
        }
    }else{
        params_init[0]=cl_float(2.002,float_format(precision));//(2+sqrt(cl_float(3,float_format(precision))))/2;
        for(int k=1;k<2*n_terms+1;k++){
            params_init[k]=cl_float(expt(cl_I(10),-1),float_format(precision)); //Initializing with some values as we have no data yet to extrapolate
        }
    }
  
    q_arr[index]=new QSC_R(params_init,params+2*j,n_terms,precision,&outstream); //Create a new instance
    q_arr[index]->optimize(max_iter,ytol,err_tol);  //Starting the optimizer
    writer(q_arr[index],outcoeffs,outdouble);    //Committing the results to the output files

    if(index==max_interpol_order+1){
        delete q_arr[0];    //Delete an old instance
        for(int k=1;k<(max_interpol_order+2);k++){
            q_arr[k-1]=q_arr[k];    
        }
    }

    for(j++;j<n_g;j++){

        outstream<<"\nCoupling g="<<double_approx(params[2*j])<<std::endl;

        index=min(j,max_interpol_order+1);

        for(int l=0;l<index;l++){
            g2[l]=expt(q_arr[index-l-1]->get_parameters()[0],2);
        }
        for(int k=0;k<2*n_terms+1;k++){
            for(int l=0;l<index;l++){
                tmp_param[l]=*q_arr[index-l-1]->get_fitCoeff(k);
            }
            params_init[k]=interpol_n(tmp_param,g2,expt(params[2*j],2),float_format(precision),min(index-1,max_interpol_order),index);
        }

        q_arr[index]=new QSC_R(params_init,params+2*j,n_terms,precision,&outstream);
        q_arr[index]->optimize(max_iter,ytol,err_tol);
        writer(q_arr[index],outcoeffs,outdouble);

        if(index==max_interpol_order+1){
            delete q_arr[0];
            for(int k=1;k<(max_interpol_order+2);k++){
                q_arr[k-1]=q_arr[k];
            }
        }

    }

    for(int k=0;k<(max_interpol_order+2);k++){
        delete q_arr[k];
    }

    outcoeffs.flush();
    outcoeffs.close();

    outdouble.flush();
    outdouble.close();

    outstream.flush();
    outstream.close();
    
    
    delete[] g2;
    delete[] params;
    delete[] params_init;

    return 0;
}
