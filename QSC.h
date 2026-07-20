
#include <cln/cln.h>
#include <cln/number.h>
#include <cln/integer.h>
#include <cln/real.h>
#include <cln/float.h>
#include <cln/io.h>
#include <cln/float_io.h>
#include <cln/malloc.h>

#include <fstream>
using std::ofstream;

using namespace cln;

class p_Function_R; // Just define the P-function class, so the Qai can have a pointer to the corresponding P.
class Qai_Function_R; //Define the Qai class so that it is known for the whole file.


class Q_Function_R {
    /*
     * So far the only reason to have these is to compute the bb coefficients.
     * We scale coefficients for u^n with \io^n to get the convention we use in Mathematica.
     */
private:
    cl_R* bb; //The leading coefficient of the \bQ-functions at large u
    const cl_R* fit_parameter; //A pointer to the exponential of the Hagedorn temperature.
    const cl_R* parameters; //The coupling and the twist parameter
    
    void update_leadCoeffs(); //This updates the coefficients bb based on the value of fit_parameter.

    cl_R glue(int i, cl_R & u,float_format_t p); //Evaluate the diagonal gluing matrix at -u.
    
public:
    Q_Function_R (const cl_R* fit_param_init , const cl_R * param_init);
    ~Q_Function_R ();
    
    void update_Coeffs();
    
    const cl_R* get_bb(); //This also calls the updating.

    void cut_mismatch (Qai_Function_R * qai , p_Function_R * p , cl_N * df , int n_points, float_format_t prec); /*
        This computes the mismatch between Q and \wt{Q} on n_points points along the cut
    */
    
};



class Qai_Function_R {
    /*
     * An object of this class represents the Qai-functions for a QSC.
     * Upon initialization the memory space for the coefficients b in the large u expansion is allocated.
     * The coefficients are computed from the coefficients of the P-functions only if update_Coeffs is called.
     * We scale coefficients for u^n with \io^(n+1) to get the convention we use in Mathematica.
     */
    
private:
    
    cl_R * b; //Upon creation we allocate memory for all these coefficients. Ordered as b_{a|i,n} -> b[4*a+i+8*n].
    cl_R * b_lead; //Upon creation we allocate memory for all these coefficients
    p_Function_R * pFunc; //Pointing to the corresponding \bP-function
    Q_Function_R * qFunc; //Pointing to the corresponding \bQ-function
    const cl_R * parameters; //contains (g,\omega)
    int n_terms; //There is n_terms subleading terms
    float_format_t precision;
    int uu=60;

    cl_R * QQmatrix; 
    /*
     * The matrix for the linear system for the b coefficients.
     */

    cl_R * QQsource; // The source for the linear system.

    cl_R * dd; //An auxiliary to save results of get_d so we avoid recomputing these factors.
    
    void update_leadCoeffs(); //The leading order coefficients have explicit expressions in terms of the parameters of the problem.

    void generate_QQequations(); //This sets the values of the QQmatrix and QQsource.

    cl_R get_d(int m, int l); //An auxiliary function

    void test_QQ();
    /*
     This tests how well the QQ-relations are solved.
    */
    
    //Functions to help evaluate close to the cut

    cl_N twisty (cl_N & iu, int i, float_format_t p); //The y^{iu t_i} factor
    cl_N twistx (cl_N & iu, int a); //The y^{iu \omega t_a} factor
    
    void evaluate_Qai_farFromCut (cl_N * Qaiu , cl_R & u); //evaluates at u+i(2*uu+1)/2 and saves it into Qaiu
    void evaluate_Qai_farFromCut (cl_N * Qaiu , cl_R & u, int x); //evaluates at u+i(2*x+1)/2 and saves it into Qaiu
    
public:
    
    Qai_Function_R ( p_Function_R * p, const cl_R * params_init, int n, float_format_t prec_init);
    
    ~Qai_Function_R ();
    
    void update_Coeffs(); //First calls update_lead, then solves the QQ-relation for the b coefficients.
    
    void evaluate_Qai_nearCut (cl_N * Qaiu, cl_R & u); //evaluates at u+i/2 and saves it into Qaiu

    Q_Function_R * get_qFunc ();

};




class p_Function_R {
    /*
     * This represents the P-functions for a QSC.
     * We scale coefficients for u^n with \io^n to get the convention we use in Mathematica.
     * Upon initialization it allocates memory for the coefficients c and the fit_parameter.
     * The coefficients ct are computed using Fourier after we compute the Qai near the cut and then compute the corresponding P's. Updating ct first updates the coefficients b of Qai if not otherwise specified.
     */
    
private:
    
    cl_R * fit_parameter; //Memory is allocated for this
    const cl_R * parameters; //The coupling and the twist parameter
    cl_R * c; //Memory is allocated for these
    cl_R* aa; //A container for the leading coefficients to skip allocating memory later on.
    int k_terms; //There is k_terms subleading terms to the \bP-functions.
    float_format_t precision;
    Qai_Function_R * QaiFunc;
    
    //This function takes the full array of QSC fit parameters and splits it into the c coefficients and the other fit parameter
    
    void update_leadCoeffs (); // Computes the values of the leading coefficients which have explicit expressions, saved both as part of c and in aa.

    cl_N zhukovsky(cl_N & u); //Gives the Zhukovsky variable. At the cut it gives the value above the cut.
    cl_N invZhukovsky(cl_N & u); //Gives the inverse Zhukovsky variable.
    
public:
    
    p_Function_R (const cl_R* params_init, const cl_R* params_fix , int K , float_format_t prec_init );
    
    ~p_Function_R ();
    
    void set_QaiFunc (Qai_Function_R * q); //To tell the \bP-function where the corresponding Qai-function is.

    //cl_N twisty (cl_N  pt, int i, float_format_t p);
    cl_N twistx (cl_N pt, int a); //y^{i pt \omega t_a}
    
    void evaluate_pLower (cl_N * Pl, cl_N & u); //Evaluates the P_a(u)
    void evaluate_ptLower (cl_N * Ptl, cl_N & u); //Evaluates the tilde(P)_a(u)
    
    void shift_coeff(cl_R & s , int n); //Adds the value of s to the n-th parameter (indexing of c starts at 0). If n equals coeff_count it shifts fit_parameter.
    void shift_coeff_minus(cl_R & s , int n); //Subtracts the value of s to the n-th parameter (indexing of c starts at 0). If n equals coeff_count it shifts fit_parameter.
    
    const cl_R* get_fitCoeffs();//This function returns a pointer to the current value of the fitted coefficients c
    const cl_R* get_fitCoeff(int i);
    const cl_R* get_fitParam();//This function returns a pointer to the parameter we are trying to fit

    const cl_R* get_aa(); //This function returns a pointer to the leading coefficient aa
    const int get_k(); //Returns the value of k_terms
    
};






class QSC_R {
    /*
     * An object of this class represents the quantum spectral curve at one fixed value of coupling, twists, charges, etc..
     * Upon initializ it creates its own instances of the P- and Qai-functions.
     * It can be asked to perform a step with the Newton method and how well the current values of the parameters work as a solution. The caller decides whether to take more steps.
     * The QSC object handles varying all the parameters, the P-function only knows their current value and cannot edit them.
     */
    
private:
    ofstream * outstream;

    cl_R * parameters; // This contains the coupling and other parameters that are held fixed, like charges and twist parameters. Memory is allocated for these
    int n_terms; // The number of subleading terms we want to keep in the Qai-functions
    float_format_t precision; //The number of digits we keep for floating point numbers
    p_Function_R * pFunc; //We set k_terms to n_terms+1
    Qai_Function_R * QaiFunc;
    Q_Function_R * qFunc;

    int n_params; //The number of parameters to fit: n_params= 2*n_terms+1
    int n_points; //The number of points on the cut that we want to consider 
    cl_R dx; //The small shift we use to compute the gradient.
    cl_N * err; // The current value of the mismatches on the cut.
    
    void gradient(cl_N * df); //Computes the gradient of cut_mismatch by shifting each parameter by dx
    void hessian(cl_N * h, cl_N * df, cl_R damper);
    /*
        A Hessian matrix computed by contracting the gradient with itself.
        We multiply the diagonal elements by (1+damper) to damp the optimization problem.
    */
    void grad_dot_err(cl_N * df, cl_N * gde); //The right-hand side of the optimization problem
    void generate_step(cl_N * h, cl_N * gde , cl_N * s); //Calls the linear system solver to compute the step to take.

    void take_step(cl_N * df, cl_N * h, cl_N * gde , cl_N * s); //Computes the gradient, figures out the step and takes the step.

    void err_norm(cl_R & tmp); //Computes the norm squared of the cut-mismatch.

    
public:


    QSC_R ();
    QSC_R (const cl_R * params_fit_init ,const cl_R * params_init , int n_init , int prec_init,ofstream * o);
    QSC_R (const cl_R * params_fit_init ,const cl_R * params_init , int n_init , int prec_init,ofstream * o,bool initialize);
    
    ~QSC_R ();

    void initialize(const cl_R * params_fit_init ,const cl_R * params_init , int n_init , int prec_init,ofstream * o);

    void optimize(int max_steps, int ytol, int errtol); //This automatizes taking the step. Iterates until both tolerances are met or the uip to max_steps.

    const cl_R * get_parameters(); //Returns a pointer to the coupling and the twist parameter

    const int get_n_params(); //Returns the number of fit parameters.
    
    const cl_R * get_fitParameter(); //Returns the fitted value of y
    const cl_R * get_fitCoeffs(); //Returns the array of fitted coefficients of \bP.
    const cl_R * get_fitCoeff(int i); //Returns the i-th fitted parameter.

};


