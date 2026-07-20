
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
#include <chrono>
#include <fstream>
using std::ofstream;

#include "QSC.h"
#include "linSolve.h"
#include "auxmathfuncs.h"

using namespace cln;


const cl_R* Q_Function_R::get_bb(){
    update_leadCoeffs();
    return bb;
};

void Q_Function_R::update_leadCoeffs(){

    const cl_R y2w = realpart(expt((*fit_parameter),2*parameters[1]));

    bb[0]=1;
    bb[1] =expt(expt((*fit_parameter),2)+y2w,2) *expt(1+expt((*fit_parameter),2)*y2w,2)/(expt(y2w,2)*expt(expt(*fit_parameter,4)-1,2));
    bb[2]=bb[0];
    bb[3]=bb[1];

};

void Q_Function_R::update_Coeffs(){
    update_leadCoeffs();
};

Q_Function_R::Q_Function_R(const cl_R* fit_param_init , const cl_R * param_init):bb(new cl_R[4]), fit_parameter(fit_param_init),parameters(param_init) {
    update_Coeffs();
};

Q_Function_R::~Q_Function_R(){
    delete[] bb;
};

cl_R Q_Function_R::glue(int i, cl_R & u , float_format_t p){
    return power(-1,i%2)*exp(-2*pi(p)*u*(1-2*(i/2)));
};

void Q_Function_R::cut_mismatch(Qai_Function_R * qai, p_Function_R * p, cl_N * df , int n_points,float_format_t prec){

    cl_N * Qev = new cl_N[2*4];
    cl_N * pp = new cl_N[4];
    cl_N * pt = new cl_N[4];
    cl_N * dftmp = new cl_N[4*n_points];

    for(int n=0;n<n_points;n++){
        cl_R u=-2*parameters[0]*cos(pi(prec)*(n+recip(cl_I(2)))/n_points);
        qai->evaluate_Qai_nearCut(Qev,u);
        p->evaluate_pLower(pp,u );
        p->evaluate_ptLower(pt,u );
        for(int i=0;i<2;i++){
            df[4*n+i]=0;
            df[4*n+i+2]=0;
            for(int a=0;a<4;a++){
                df[4*n+i]-=power(-1,3-a)*pp[3-a]*Qev[2*a+i];
                df[4*n+i+2]-=power(-1,3-a)*pt[3-a]*Qev[2*a+i];
            };
            dftmp[4*(n_points-n-1)+i]=glue(i,u,prec)*df[4*n+i];
            dftmp[4*(n_points-n-1)+i+2]=glue(i,u,prec)*df[4*n+i+2];
        };

    };
    
    for(int n=0;n<4*n_points;n++){
        df[n]/=dftmp[n];
        df[n]-=1;
    }

    delete[] Qev;
    delete[] pp;
    delete[] pt;
    delete[] dftmp;


};








Qai_Function_R::Qai_Function_R ( p_Function_R * p, const cl_R * params_init, int n, float_format_t prec_init) : pFunc(p), parameters(params_init), n_terms(n), precision (prec_init) {
    
    b = new cl_R[2*4*(n_terms)];
    b_lead = new cl_R[2*4];
    QQmatrix = new cl_R[64*n_terms*n_terms];
    QQsource = new cl_R[2*4*n_terms];
    qFunc = new Q_Function_R(pFunc->get_fitParam(),parameters);

    dd= new cl_R[(n_terms+8)*((n_terms+8)/2)];

    for(int m=2;m>-n_terms-6;m--){
        for(int l=0;l<(n_terms+8)/2;l++){
            dd[(2-m)*((n_terms+8)/2) +l]=d(m,l,parameters);
        }
    }

}; 

Qai_Function_R::~Qai_Function_R() {

    delete[] b;
    delete[] b_lead;
    delete[] QQmatrix;
    delete[] QQsource;
    delete qFunc;
    delete[] dd;

};

cl_R Qai_Function_R::get_d(int m, int l){
    return dd[(2-m)*((n_terms+8)/2) +l];
}

void Qai_Function_R::generate_QQequations(){

    const cl_R * c = pFunc->get_fitCoeffs();
    const cl_R * ys = pFunc->get_fitParam();
    const int k_terms=pFunc->get_k();
    cl_R tmp;
    cl_R tmp2;
    cl_R tmp_bis;

    for(int i=0;i<64*n_terms*n_terms;i++){
        QQmatrix[i]=0;
    }

    for(int i=0;i<8*n_terms;i++){
        QQsource[i]=0;
    }

    for(int j=0;j<4;j++){
        QQmatrix[4+j]+=power(-1,j)* b_lead[4+(j ^ 1)];
    };
    

    //All the (2,i) contributions to the equations.
    for( int bb=0;bb<4;bb++){

        for(int m=2-(bb%2);m>-n_terms-(bb%2)+1;m--){

            tmp=0;
            for(int k=max(-1,-m-k_terms-(bb%2)+1);k<-m+2-(bb%2)&&k<k_terms;k++){
                tmp+=power(-1,(1-(bb/2))*(-m-k+1-(bb%2))+3-bb)*c[k_terms+2+k]*c[(1-bb%2)*(k_terms+2)-m-k];
            }

            tmp/=expt(parameters[0],2*m);
            tmp*=realpart(expt(*ys,parameters[1]*(bb/2)*2));

            for(int i=0;i<4;i++){

                tmp_bis=power(-1,(bb/2))*b_lead[(bb%2)*4+(i ^((bb/2)<<1))]*tmp;

                for(int n=(bb%2)+(i%2); n>-1 && n>-n_terms-m+1+(i%2) ; n--){

                    tmp2=binomial((bb%2)+(i%2),-n+(bb%2)+(i%2))*power(-2,n-(bb%2)-(i%2))*tmp_bis;
                    if(1+(i%2)>=n+m){
                        QQsource[(1+i*n_terms+1+(i%2)-n-m)]-=get_d(m,0)*tmp2;
                    }

                    for(int l=1;2*l-n-m<n_terms-(i%2)-1;l++){
                        if(1+(i%2)>=n+m-2*l){
                            QQsource[(1+i*n_terms+1+(i%2)-n-m+2*l)]-=get_d(m,l)*tmp2;
                        }
                    }
                    

                    for(int p=1;p<(bb%2)+(i%2)+1 &&p<-n+(bb%2)+(i%2)+1 && p<n_terms+1;p++){

                        tmp2=power(-1,(bb/2)*(p+1))*binomial(-p+(bb%2)+(i%2),-p-n+(bb%2)+(i%2))*power(-2,n+p-(bb%2)-(i%2))*tmp;
                        if(1+(i%2)>=n+m){
                            QQmatrix[(1+i*n_terms+1+(i%2)-n-m)*8*n_terms+(bb%2)*4+(i ^((bb/2)<<1))+(p-1)*8]+=get_d(m,0)*tmp2; 
                        }

                        for(int l=1;2*l-n-m<n_terms-(i%2)-1;l++){
                            if(1+(i%2)>=n+m-2*l){
                                QQmatrix[(1+i*n_terms+1+(i%2)-n-m+2*l)*8*n_terms+(bb%2)*4+(i ^((bb/2)<<1))+(p-1)*8]+=get_d(m,l)*tmp2;
                            }
                        }

                    }

                }

                for(int n=(bb%2)+(i%2);n>-n_terms-m+1+(i%2);n--){

                    for(int p=(bb%2)+(i%2)+1;p<-n+(bb%2)+(i%2)+1 && p<n_terms+1;p++){

                        tmp2=power(-1,(bb/2)*(p+1))*binomial(-p+(bb%2)+(i%2),-p-n+(bb%2)+(i%2))*power(-2,n+p-(bb%2)-(i%2))*tmp;
                        if(n+m<=1+(i%2)){
                            QQmatrix[(1+i*n_terms+1+(i%2)-n-m)*8*n_terms+(bb%2)*4+(i ^((bb/2)<<1))+(p-1)*8]+=get_d(m,0)*tmp2; 
                        }

                        for(int l=1;2*l-n-m<n_terms-(i%2)-1;l++){
                            if(n+m-2*l<=1+(i%2)){
                                QQmatrix[(1+i*n_terms+1+(i%2)-n-m+2*l)*8*n_terms+(bb%2)*4+(i ^((bb/2)<<1))+(p-1)*8]+=get_d(m,l)*tmp2;
                            }
                        }

                    }

                }
            }
        }
    }

    for(int i=0;i<4;i++){
        tmp=realpart(expt(*ys,2*(parameters[1]+(1-2*(i/2)))));

        for(int m=(i%2);m>-n_terms+(i%2)+1;m--){
            for(int l=1;l<-m+2+(i%2)&&l<=n_terms;l++){
                QQmatrix[(1+i*n_terms+1+(i%2)-m)*8*n_terms+8*(l-1)+4+i]+=binomial(-l+1+(i%2),-l-m+1+(i%2))*(power(-1,-l-m+1+(i%2))+tmp)*power(2,l+m-1-(i%2));
            }
        }

        for(int m=1+(i%2);m>-1;m--){

            QQsource[(1+i*n_terms+1+(i%2)-m)]-=binomial(1+(i%2),-m+1+(i%2))*(power(-1,-m+1+(i%2))+tmp)*power(2,m-1-(i%2))*b_lead[4+i];
            
        }

    }

    
    //All the (4,i) contributions to the equations.
    for( int bb=0;bb<4;bb++){
        for(int m=2-(bb%2);m>-n_terms-(bb%2);m--){

            tmp=0;
            for(int k=max(-1,-m-k_terms-(bb%2)+1);k<-m+2-(bb%2)&&k<k_terms;k++){
                tmp+=power(-1,(1-bb/2)*(-m-k+1-(bb%2))+k+1+3-bb)*c[k_terms+2+k]*c[(1-bb%2)*(k_terms+2)-m-k];
            }
            tmp/=expt(parameters[0],2*m);
            tmp*=realpart(expt(*ys,parameters[1]*(-1+bb/2)*2));


            for(int i=0;i<2;i++){

                tmp_bis=power(-1,(bb/2))*b_lead[(bb%2)*4+(i ^((bb/2)<<1))]*tmp;

                for(int n=(bb%2)+(i%2); n>-1 && n>-n_terms-m+(i%2) ;n--){

                    tmp2=binomial((bb%2)+(i%2),-n+(bb%2)+(i%2))*power(-2,n-(bb%2)-(i%2))*tmp_bis;
                    if((i%2)-n-m>=0){
                        QQsource[(1+(4+i)*n_terms+(i%2)-n-m)]-=get_d(m,0)*tmp2; 
                    }

                    for(int l=1;2*l-n-m<n_terms-(i%2);l++){
                        if((i%2)-n-m+2*l>=0){
                            QQsource[(1+(4+i)*n_terms+(i%2)-n-m+2*l)]-=get_d(m,l)*tmp2;
                        }
                    }

                    for(int p=1;p<(bb%2)+(i%2)+1&&p<-n+(bb%2)+(i%2)+1 &&p<n_terms+1;p++){

                        tmp2=power(-1,(bb/2)*(p+1))*binomial(-p+(bb%2)+(i%2),-p-n+(bb%2)+(i%2))*power(-2,n+p-(bb%2)-(i%2))*tmp;
                        if((i%2)-n-m>-1){
                            QQmatrix[(1+(4+i)*n_terms+(i%2)-n-m)*8*n_terms+(bb%2)*4+(i ^((bb/2)<<1))+(p-1)*8]+=get_d(m,0)*tmp2; 
                        }

                        for(int l=1;2*l-n-m<n_terms-(i%2);l++){
                            if((i%2)-n-m+2*l>-1){
                                QQmatrix[(1+(4+i)*n_terms+(i%2)-n-m+2*l)*8*n_terms+(bb%2)*4+(i ^((bb/2)<<1))+(p-1)*8]+=get_d(m,l)*tmp2;
                            }
                        }
                    }
                }

                for(int n=(bb%2)+(i%2);n>-n_terms-m+(i%2);n--){

                    for(int p=(bb%2)+(i%2)+1;p<-n+(bb%2)+(i%2)+1&&p<n_terms+1;p++){

                        tmp2=power(-1,(bb/2)*(p+1))*binomial(-p+(bb%2)+(i%2),-p-n+(bb%2)+(i%2))*power(-2,n+p-(bb%2)-(i%2))*tmp;
                        if((i%2)-n-m>-1){
                            QQmatrix[(1+(4+i)*n_terms+(i%2)-n-m)*8*n_terms+(bb%2)*4+(i ^((bb/2)<<1))+(p-1)*8]+=get_d(m,0)*tmp2; 
                        }

                        for(int l=1;2*l-n-m<n_terms-(i%2);l++){
                            if((i%2)-n-m+2*l>-1){
                                QQmatrix[(1+(4+i)*n_terms+(i%2)-n-m+2*l)*8*n_terms+(bb%2)*4+(i ^((bb/2)<<1))+(p-1)*8]+=get_d(m,l)*tmp2;
                            }
                        }
                    }
                }
            }
        }
    }

    for(int i=0;i<2;i++){

        tmp=realpart(expt(*ys,2*(-parameters[1]+(1-2*(i/2)))));

        for(int m=(i%2);m>-n_terms+(i%2);m--){
            for(int l=1;l<-m+2+(i%2)&&l<=n_terms;l++){
                QQmatrix[(1+(4+i)*n_terms+(i%2)-m)*8*n_terms+8*(l-1)+4+(i ^ 2)]+=binomial(-l+1+(i%2),-l-m+1+(i%2))*(power(-1,-l-m+1+(i%2))+tmp)*power(2,l+m-1-(i%2))*power(-1,l+1);
            }
        }

        for(int m=(i%2);m>-1;m--){
            
            QQsource[(1+(4+i)*n_terms+(i%2)-m)]+=binomial(1+(i%2),-m+1+(i%2))*(power(-1,-m+1+(i%2))+tmp)*power(2,m-1-(i%2))*b_lead[4+(i ^ 2)];
            
        }

    }

    
    //All contributions for (1,i)
    for( int bb=0;bb<4;bb++){
        for(int m=1-(bb%2);m>-n_terms-(bb%2)-4;m--){

            tmp=0;
            for(int k=max(0,-m-k_terms-(bb%2)+1);k<-m+2-(bb%2)&&k<k_terms+1;k++){
                tmp+=power(-1,(1-bb/2)*(-m-k+1-(bb%2))+3-bb)*c[k]*c[(1-bb%2)*(k_terms+2)-m-k];
            }
            tmp/=expt(parameters[0],2*m);
            tmp*=realpart(expt(*ys,parameters[1]*(bb/2)*2));

            for(int i=0;i<2;i++){

                tmp_bis=power(-1,(bb/2))*b_lead[(bb%2)*4+(i ^((bb/2)<<1))]*tmp;

                for(int n=(bb%2)+(i%2);n>-1&&n>-n_terms-m-4+(i%2);n--){

                    tmp2=binomial((bb%2)+(i%2),-n+(bb%2)+(i%2))*power(-2,n-(bb%2)-(i%2))*tmp_bis;
                    if(-n-m>3){
                        QQsource[(1+(6+i)*n_terms-4-n-m)]-=get_d(m,0)*tmp2; 
                    }

                    for(int l=1;2*l-n-m<n_terms+4-(i%2);l++){
                        if(-n-m+2*l>3){
                            QQsource[(1+(6+i)*n_terms-4-n-m+2*l)]-=get_d(m,l)*tmp2;
                        }
                    }

                    for(int p=1;p<+(bb%2)+(i%2)+1&&p<-n+(bb%2)+(i%2)+1 &&p<n_terms+1;p++){

                        tmp2=power(-1,(bb/2)*(p+1))*binomial(-p+(bb%2)+(i%2),-p-n+(bb%2)+(i%2))*power(-2,n+p-(bb%2)-(i%2))*tmp;

                        if(-n-m>3){
                            QQmatrix[(1+(6+i)*n_terms-4-n-m)*8*n_terms+(bb%2)*4+(i ^((bb/2)<<1))+(p-1)*8]+=get_d(m,0)*tmp2;
                        }

                        for(int l=1;2*l-n-m<n_terms+4-(i%2);l++){
                            if(-n-m+2*l>3){
                                QQmatrix[(1+(6+i)*n_terms-4-n-m+2*l)*8*n_terms+(bb%2)*4+(i ^((bb/2)<<1))+(p-1)*8]+=get_d(m,l)*tmp2;
                            }
                        }
                    }
                }

                for(int n=(bb%2)+(i%2);n>-n_terms-m-4+(i%2);n--){

                    for(int p=(bb%2)+(i%2)+1;p<-n+(bb%2)+(i%2)+1&&p<n_terms+1;p++){

                        tmp2=power(-1,(bb/2)*(p+1))*binomial(-p+(bb%2)+(i%2),-p-n+(bb%2)+(i%2))*power(-2,n+p-(bb%2)-(i%2))*tmp;

                        if(-n-m>3){
                            QQmatrix[(1+(6+i)*n_terms-4-n-m)*8*n_terms+(bb%2)*4+(i ^((bb/2)<<1))+(p-1)*8]+=get_d(m,0)*tmp2;
                        }

                        for(int l=1;2*l-n-m<n_terms+4-(i%2);l++){
                            if(-n-m+2*l>3){
                                QQmatrix[(1+(6+i)*n_terms-4-n-m+2*l)*8*n_terms+(bb%2)*4+(i ^((bb/2)<<1))+(p-1)*8]+=get_d(m,l)*tmp2;
                            }
                        }
                    }
                }
            }
        }
    }

    for(int i=0;i<2;i++){

        tmp=realpart(expt(*ys,2*(parameters[1]+(1-2*(i/2)))));

        for(int m=-4;m>-n_terms-4+(i%2);m--){
            for(int l=1;l<-m+1+(i%2)&&l<=n_terms;l++){
                QQmatrix[(1+(6+i)*n_terms-4-m)*8*n_terms+8*(l-1)+i]+=binomial(-l+(i%2),-l-m+(i%2))*(power(-1,-l-m+(i%2))+tmp)*power(2,l+m-(i%2));
            }
        }

    }

}

void Qai_Function_R::test_QQ(){

    cl_N * qaip=new cl_N[16];
    cl_N * qaim=new cl_N[16];
    cl_N * p=new cl_N[4];

    cl_R u=0;

    evaluate_Qai_farFromCut(qaip,u,100);
    evaluate_Qai_farFromCut(qaim,u,99);
    cl_N u2=complex(0,cl_float(100,precision));
    std::cerr<<"evaluating p\n";
    pFunc->evaluate_pLower(p,u2);
    std::cerr<<"evaluated p\n";

    cl_N tmp=0;
    for(int a=0;a<4;a++){
        for(int aa=0;aa<4;aa++){
            tmp=0;
            for(int i=0;i<4;i++){
                tmp+=qaip[4*a+i]*qaip[4*(3-aa)+3-i]*power(-1,3-aa+3-i);
            }
            std::cerr<<cl_float(abs(tmp),float_format(15))<<"  ";
        }
        std::cerr<<"\n\n";
    }
    for(int i=0;i<4;i++){
        for(int j=0;j<4;j++){
            tmp=0;
            for(int a=0;a<4;a++){
                tmp+=qaip[4*a+i]*qaip[4*(3-a)+3-j]*power(-1,3-a+3-j);
            }
            std::cerr<<cl_float(abs(tmp),float_format(15))<<"  ";
        }
        std::cerr<<"\n\n";
    }

    cl_N qi;

    for(int i=0;i<4;i++){
        qi=0;
        for(int bb=0;bb<4;bb++){
            qi-=power(-1,3-bb)*p[3-bb]*qaip[bb*4+i];
        }
        for(int a=0;a<4;a++){
            std::cerr<<cl_float(abs(1+(-p[a]*qi-qaim[a*4+i])/qaip[4*a+i]),float_format(10))<<" ";//<<cl_float(realpart(1+(-p[a]*qi-qaim[a*4+i])/qaip[4*a+i]),float_format(10))<<"+"<<cl_float(imagpart(1+(-p[a]*qi-qaim[a*4+i])/qaip[4*a+i]),float_format(15))<<"i ";
        }
        std::cerr<<"\n";
    }

    delete[] qaip;
    delete[] qaim;
    delete[] p;

}

void Qai_Function_R::update_Coeffs(){

    update_leadCoeffs();
    //std::cerr<<"generating equations"<<std::endl;
    generate_QQequations();


    //std::cerr<<"solving equations"<<std::endl;
    linsolveHouseholder(QQmatrix,QQsource,b,8*n_terms);

    //test_QQ();

};

void Qai_Function_R::update_leadCoeffs(){
    
    const cl_R* ys=pFunc->get_fitParam();
    const cl_R* bb=qFunc->get_bb();
    const cl_R* aa=pFunc->get_aa();

    cl_R yw = realpart(expt(*ys,parameters[1]+1));
    
    cl_R tmp = recip(yw+recip(yw));
    
    for(int i=0;i<2;i++){
        for(int a=0;a<2;a++){
            b_lead[(a*4+i)]=-tmp*aa[a]*bb[i];
        };
    };
    
    yw = realpart(expt(*ys,parameters[1]-1));
    tmp = recip(yw+recip(yw));
    
    for(int i=2;i<4;i++){
        for(int a=0;a<2;a++){
            b_lead[(a*4+i)]=tmp*aa[a]*bb[i];
        };
    };
    
};

cl_N Qai_Function_R::twisty (cl_N  & iu, int i, float_format_t p){
    return exp(-complex(0,pi(p))*iu)*expt(*pFunc->get_fitParam(),2*iu*(1-2*(i/2)));
}

cl_N Qai_Function_R::twistx (cl_N & iu, int a){
    return expt(*pFunc->get_fitParam(),2*iu*(1-2*(a/2))*parameters[1]);
}

void Qai_Function_R::evaluate_Qai_farFromCut(cl_N * Qaiu , cl_R & u){
    
    cl_N ipt=complex(0,u)-cl_RA(2*uu+1)/2;

    for(int i=0;i<2;i++){
        for(int a=0;a<4;a++){

            Qaiu[2*a+i]=0;

            for(int n=n_terms;n>0;n--){

                Qaiu[2*a+i]+=power(-1,(a/2)*(n+1))*expt(ipt,(a%2)+(i%2)-n)*b[8*(n-1)+4*(a%2)+(i^((a/2)<<1))];

            };

            Qaiu[2*a+i]+=power(-1,a/2)*expt(ipt,(a%2)+(i%2))*b_lead[4*(a%2)+(i^((a/2)<<1))];

            Qaiu[2*a+i]*=complex(0,1)*(twisty(ipt,i,precision))*(twistx(ipt,a));

        };
    };

};

void Qai_Function_R::evaluate_Qai_farFromCut(cl_N * Qaiu , cl_R & u, int x){
    
    cl_N ipt=complex(0,u)-cl_RA(2*x+1)/2;

    for(int i=0;i<4;i++){
        for(int a=0;a<4;a++){

            Qaiu[4*a+i]=0;

            for(int n=n_terms;n>0;n--){

                Qaiu[4*a+i]+=power(-1,(a/2)*(n+1))*expt(ipt,(a%2)+(i%2)-n)*b[8*(n-1)+4*(a%2)+(i^((a/2)<<1))];

            };
            Qaiu[4*a+i]+=power(-1,a/2)*expt(ipt,(a%2)+(i%2))*b_lead[4*(a%2)+(i^((a/2)<<1))];

            Qaiu[4*a+i]*=complex(0,1)*(twisty(ipt,i,precision))*(twistx(ipt,a));

        };
    };

};

void Qai_Function_R::evaluate_Qai_nearCut(cl_N * Qaiu,cl_R & u){

    cl_N * p= new cl_N[4];
    cl_N qi;
    cl_N pt=u+complex(0,uu+1);

    evaluate_Qai_farFromCut(Qaiu,u);

    for(int v=uu;v>0;v--){
        pt-=complex(0,1);
        pFunc->evaluate_pLower(p,pt);
        for (int i=0;i<2;i++){
            qi=0;
            for(int bb=0;bb<4;bb++){
                qi-=power(-1,3-bb)*p[3-bb]*Qaiu[2*bb+i]; 
            };
            for(int a=0;a<4;a++){
                Qaiu[2*a+i]-=p[a]*qi;
            };
        };
    };

    delete[] p;

};

Q_Function_R * Qai_Function_R::get_qFunc(){
    return qFunc;
};







QSC_R::QSC_R(){

};

QSC_R::QSC_R (const cl_R * params_fit_init , const cl_R * params_init , int n_init , int prec_init,ofstream * o):  n_terms(n_init) , precision(float_format(prec_init)) ,outstream(o) {
    

    parameters = new cl_R[2];
    for(int i=0;i<2;i++){
        parameters[i]=params_init[i];
    };
    
    pFunc = new p_Function_R(params_fit_init,parameters,n_terms+1,precision);

    QaiFunc = new Qai_Function_R(pFunc, parameters ,n_terms,precision);
    pFunc->set_QaiFunc(QaiFunc);
    qFunc=QaiFunc->get_qFunc();

    n_params=2*pFunc->get_k()-1;
    n_points=(n_params/2)*2+2;
    dx=expt(cl_float(10,precision),-100);
    err = new cl_N[4*n_points];

    QaiFunc->update_Coeffs();
    *outstream<<"Computing initial mismatch"<<std::endl;
    qFunc->cut_mismatch(QaiFunc,pFunc,err,n_points,precision);
    
};

QSC_R::QSC_R (const cl_R * params_fit_init , const cl_R * params_init , int n_init , int prec_init,ofstream * o,bool initialize):  n_terms(n_init) , precision(float_format(prec_init)) ,outstream(o) {
    

    parameters = new cl_R[2];
    for(int i=0;i<2;i++){
        parameters[i]=params_init[i];
    };
    
    pFunc = new p_Function_R(params_fit_init,parameters,n_terms+1,precision);

    QaiFunc = new Qai_Function_R(pFunc, parameters ,n_terms,precision);
    pFunc->set_QaiFunc(QaiFunc);
    qFunc=QaiFunc->get_qFunc();

    n_params=2*pFunc->get_k()-1;
    n_points=(n_params/2)*2+2;
    dx=expt(cl_float(10,precision),-100);
    err = new cl_N[4*n_points];

    if(initialize){
        QaiFunc->update_Coeffs();
        *outstream<<"Computing initial mismatch"<<std::endl;
        qFunc->cut_mismatch(QaiFunc,pFunc,err,n_points,precision);
        cl_R tmp;
        err_norm(tmp);
    }
    
};

QSC_R::~QSC_R (){
    delete[] parameters;
    delete pFunc;
    delete QaiFunc;

    delete[] err;

};

void QSC_R::initialize(const cl_R * params_fit_init ,const cl_R * params_init , int n_init , int prec_init,ofstream * o){
    n_terms=n_init;
    precision=float_format(prec_init);
    outstream=o;

    parameters = new cl_R[2];
    for(int i=0;i<2;i++){
        parameters[i]=params_init[i];
    };
    
    pFunc = new p_Function_R(params_fit_init,parameters,n_terms+1,precision);

    QaiFunc = new Qai_Function_R(pFunc, parameters ,n_terms,precision);
    pFunc->set_QaiFunc(QaiFunc);
    qFunc=QaiFunc->get_qFunc();

    n_params=2*pFunc->get_k()-1;
    n_points=(n_params/2)*2+2;
    dx=expt(cl_float(10,precision),-100);
    err = new cl_N[4*n_points];

    QaiFunc->update_Coeffs();
    *outstream<<"Computing initial mismatch"<<std::endl;
    qFunc->cut_mismatch(QaiFunc,pFunc,err,n_points,precision);

};

const cl_R * QSC_R::get_parameters(){
    return parameters;
};

const int QSC_R::get_n_params(){
    return n_params;
}

const cl_R * QSC_R::get_fitParameter(){
    return pFunc->get_fitParam();
};

const cl_R * QSC_R::get_fitCoeffs(){
    return pFunc->get_fitCoeffs();
};

const cl_R * QSC_R::get_fitCoeff(int i){
    return pFunc->get_fitCoeff(i);
}

void QSC_R::gradient(cl_N * df){


    for(int i = 0; i<n_params ; i++){

        //*outstream<<"\rComputing gradient: "<<i<<"/"<<n_params<<"  ";
        pFunc->shift_coeff(dx,i);

        QaiFunc->update_Coeffs();
        
        qFunc->cut_mismatch(QaiFunc,pFunc,df+i*4*n_points,n_points,precision);

        for(int j=0;j<4*n_points;j++){
            df[i*4*n_points+j]-=err[j];
            df[i*4*n_points+j]/=(dx);
        };

        pFunc->shift_coeff_minus(dx,i);

    }

};

void QSC_R::hessian(cl_N * h, cl_N * df, cl_R damper){
    
    for(int i=0 ; i< n_params;i++){
        for(int j=0 ; j<n_params;j++){
            h[i*n_params+j]=0;
            for(int k=0 ; k<4*n_points;k++){
                h[i*n_params+j]+=conjugate(df[i*4*n_points+k])*df[j*4*n_points+k];
            };
        };
        h[i*n_params+i]*=cl_float(1+damper,precision);
        
    };

};

void QSC_R::grad_dot_err(cl_N * df, cl_N * gde){
    for(int i=0;i<n_params;i++){
        gde[i]=0;
        for(int j=0;j< 4*n_points;j++){
            gde[i]-=conjugate(df[i*4*n_points+j])*err[j];
        };
    };
};

void QSC_R::generate_step(cl_N * h, cl_N * gde , cl_N * s){

    linsolveHouseholder_C(h,gde,s,n_params);

};

void QSC_R::take_step(cl_N * df, cl_N * h, cl_N * gde , cl_N * s){

    auto start = std::chrono::high_resolution_clock::now();

    gradient(df);
    *outstream<<"\rGradient computation duration "<<(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start)).count()<<"ms"<<std::endl;
    
    
    cl_R tmp=0;
    for(int i=0;i<4*n_points;i++){
        tmp+=realpart(conjugate(err[i])*err[i]);
    }
    if(tmp<0.0001*(1-parameters[1])*1.5){//0.1*(1-parameters[1])*1.5){
        tmp=0;
    }else{
        if(tmp<0.01*(1-parameters[1])*1.5){//1*(1-parameters[1])*1.5){
            tmp*=2*expt(cl_I(10),-10); //0.0000001;
        }else{
            tmp*=2*expt(cl_I(10),-6); //0.00001;
        }
    }
    
    hessian(h,df,tmp);

    grad_dot_err(df,gde);

    generate_step(h,gde,s);//Compute the step to take now


    for(int i =0; i<n_params;i++){
        //Take the step
        tmp=realpart(s[i]);
        pFunc->shift_coeff(tmp,i);
    };

    *outstream<<cl_float(1/(4*realpart(log(*(pFunc->get_fitParam())))),float_format(25))<<std::endl;

    QaiFunc->update_Coeffs();
    qFunc->cut_mismatch(QaiFunc,pFunc,err,n_points,precision); //Compute the new error vector.

};

void QSC_R::err_norm(cl_R & tmp){
    tmp=0;
    for(int i=0;i<4*n_points;i++){
        tmp+=realpart(conjugate(err[i])*err[i]);
    }
    *outstream<<"Error norm: "<<cl_float(tmp,float_format(20))<<std::endl;
}

void QSC_R::optimize(int max_steps, int toly , int tolnorm){

    cl_N * df = new cl_N[n_params*4*n_points];
    cl_N * h = new cl_N[n_params*n_params];
    cl_N * gde = new cl_N[n_params];
    cl_N * s = new cl_N[n_params];

    cl_R err_norm_new;
    err_norm(err_norm_new);

    take_step(df,h,gde,s);
    err_norm(err_norm_new);

    int step_count=1;

    while(!((abs(s[0]/(*pFunc->get_fitCoeff(0)))<expt(cl_RA(10),toly)) && (err_norm_new<expt(cl_RA(10),tolnorm))) && (step_count<max_steps)){
        take_step(df,h,gde,s);
        err_norm(err_norm_new);
        ++step_count;
    }


    delete[] df;
    delete[] h;
    delete[] gde;
    delete[] s;
};








p_Function_R::p_Function_R(const cl_R * params_init , const cl_R* params_fix  , int K , float_format_t prec_init ): parameters(params_fix), k_terms(K), precision(prec_init), aa(new cl_R[4]) {
    
    fit_parameter = new cl_R;
    *fit_parameter = params_init[0];
    
    c = new cl_R[2*(k_terms+1)];

    for(int i=0;i<k_terms;i++){

        c[i+1]=params_init[i+1];

    };

    for(int i=0;i<k_terms-2;i++){

        c[k_terms+i+4]=params_init[i+k_terms+1];

    };
    
    update_leadCoeffs();
    
};

p_Function_R::~p_Function_R() {
    delete fit_parameter;
    delete[] c;
    delete[] aa;
};

void p_Function_R::set_QaiFunc(Qai_Function_R * q){
    QaiFunc=q;
};

void p_Function_R::update_leadCoeffs(){

    const cl_R yw2 = realpart(expt(*fit_parameter,2*parameters[1]));
    const cl_R ys = *fit_parameter;
    const cl_R ys2 = expt(*fit_parameter,2);
    const cl_R yw4 = expt(yw2,2);
    const cl_R ys4 = expt(*fit_parameter,4);
    const cl_R yw6 = expt(yw2,3);
    const cl_R ys6 = expt(*fit_parameter,6);
    const cl_R yw8 = expt(yw2,4);
    const cl_R ys8 = expt(*fit_parameter,8);
    const cl_R yw10 = expt(yw2,5);
    const cl_R yw12 = expt(yw2,6);

    aa[0] = (ys2+yw2)*(1+ys2*yw2)/(ys2*(expt(yw2,2)-1));

    for(int a=1;a<4;a++){
        aa[a]=aa[0];
    };

    for(int a=0;a<2;a++){
        c[a*(k_terms+1)]=expt(parameters[0],2*(a%2))*aa[a];
    };

    c[k_terms+2]= -(-ys2+yw2+yw6+ys4*yw2+6*ys2*yw4+ys4*yw6-ys2*yw8)/(2*ys2*expt(yw4-1,2))+expt(parameters[0],2)*c[1];
    c[k_terms+3]= 2*yw2*(1+ys4+yw4+4*ys2*yw2+ys4*yw4)*(expt(ys2+yw2+ys4*yw2+ys2*yw4,2)+expt(parameters[0],2)*ys4*expt(yw4-1,2)*expt(c[1],2));
    c[k_terms+3]-=2*expt(parameters[0],2)*ys2*(yw4-1)*(ys4+(3*(yw2+yw10)+10*yw6)*(ys2+ys6)+(yw4+yw8)*(2+11*ys4+2*ys8)+ys4*yw12)*c[2];
    c[k_terms+3]+=c[1]*(-ys2*yw2*((yw2+yw10)*(-1+10*ys4-ys8)+(1+7*(yw4+yw8)+yw12)*(ys2+ys6)+2*yw6*(3+2*ys4+3*ys8))-4*expt(parameters[0],4)*ys4*expt(yw4-1,3)*(ys2+yw2+ys4*yw2+ys2*yw4)*c[2]);
    c[k_terms+3]+=4*expt(parameters[0],2)*ys2*expt(yw4-1,2)*expt(ys2+yw2+ys4*yw2+ys2*yw4,2)*(expt(parameters[0],2)*c[3]-c[k_terms+1+3]);
    c[k_terms+3]/=(ys2+yw2)*(1+ys2*yw2)*(-2*ys2*yw2*(yw4-1)*(4*ys2*yw2+(1+yw4)*(1+ys4))-4*expt(parameters[0],2)*ys4*expt(-1+yw4,3)*c[1]);

    //std::cerr<< ys<<std::endl;
    //std::cerr<<c[k_terms+2]<<std::endl;
    //std::cerr<<c[k_terms+3]<<std::endl;
    
};

cl_N p_Function_R::twistx(cl_N pt,int a){
    return expt(*fit_parameter,complex(0,2)*pt*parameters[1]*(1-2*(a/2)));
};

cl_N p_Function_R::zhukovsky(cl_N & u){
    return (u+complex(0,1)*sqrt(4*parameters[0]*parameters[0]-(u*u)))/(2*parameters[0]);
};

cl_N p_Function_R::invZhukovsky(cl_N & u){
    return (u+complex(0,-1)*sqrt(4*parameters[0]*parameters[0]-(u*u)))/(2*parameters[0]);
};

void p_Function_R::evaluate_pLower(cl_N * Pl, cl_N & u){

    cl_N xx = complex(0,1)*zhukovsky(u)/parameters[0];
    cl_N tmp;
    cl_N tmp2;

    Pl[0]=c[0];
    Pl[1]=c[k_terms+1]*xx;
    Pl[2]=c[0];
    Pl[3]=Pl[1];

    tmp=1;

    for(int n=1;n<k_terms+1;n++){
        tmp2=c[k_terms+1+n]*tmp;
        Pl[1]+=tmp2;
        Pl[3]+=power(-1,n)*tmp2;
        tmp/=xx;
        tmp2=c[n]*tmp;
        Pl[0]+=tmp2;
        Pl[2]+=power(-1,n)*tmp2;
    }
    for(int a=0;a<4;a++){
        Pl[a]*=twistx(u,a);
    }

};

void p_Function_R::evaluate_ptLower(cl_N * Ptl, cl_N & u){

    cl_N xx = complex(0,1)*invZhukovsky(u)/parameters[0];
    cl_N tmp;
    cl_N tmp2;

    Ptl[0]=c[0];
    Ptl[1]=c[k_terms+1]*xx;
    Ptl[2]=c[0];
    Ptl[3]=Ptl[1];

    tmp=1;

    for(int n=1;n<k_terms+1;n++){
        tmp2=c[k_terms+1+n]*tmp;
        Ptl[1]+=tmp2;
        Ptl[3]+=power(-1,n)*tmp2;
        tmp/=xx;
        tmp2=c[n]*tmp;
        Ptl[0]+=tmp2;
        Ptl[2]+=power(-1,n)*tmp2;
    }
    for(int a=0;a<4;a++){
        Ptl[a]*=twistx(u,a);
    }

};

void p_Function_R::shift_coeff(cl_R & s, int n){
    if(n!=0){
        if(n<=k_terms){
            c[n]+= (s);
            if(n<4){
                update_leadCoeffs();
            }
        }
        else{
            c[n+3]+= (s);
            if(n==k_terms+1){
                update_leadCoeffs();
            }
        }
    }
    else{
        *fit_parameter+=(s);
        update_leadCoeffs();
    }
}

void p_Function_R::shift_coeff_minus(cl_R & s, int n){
    if(n!=0){
        if(n<=k_terms){
            c[n]-= (s);
            if(n<4){
                update_leadCoeffs();
            }
        }
        else{
            c[n+3]-= (s);
            if(n==k_terms+1){
                update_leadCoeffs();
            }
        }
    }
    else{
        *fit_parameter-=(s);
        update_leadCoeffs();
    }
}

const cl_R* p_Function_R::get_fitCoeffs(){
    return c;
};

const cl_R * p_Function_R::get_fitCoeff(int i){
    if(i!=0){
        if(i<=k_terms){
            return c+i;
        }else{
            return c+i+3;
        }
    }else{
        return fit_parameter;
    }
}

const cl_R* p_Function_R::get_fitParam(){

    return fit_parameter;

};

const cl_R* p_Function_R::get_aa(){
    
    return aa;

};

const int p_Function_R::get_k(){
    return k_terms;
}
