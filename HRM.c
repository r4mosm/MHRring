/*

    MHR Model

    solver: gsl Runge-Kutta-Fehlberg (4, 5) method
    intermediate results are printed after each successful step

    input:
        [k] [sigma] [total time]
 
    output:
        file 1: x(t)

    author: Marlon Ramos, 02 Jul 2023
    email: marlon.ramos@fat.uerj.br
 
 */


#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <math.h>

int Nosc=100;

const gsl_rng_type * T;
gsl_rng * r;

unsigned int seed_dev(void);
double seed;

long double *x, *y, *z, *phi, *C, *J, *theta, *SM;

int dfunc(double t, const double y[], double f[], void *params_ptr);

void order(FILE *fp, double t);
void SI(FILE *fp, double t);

int main(int argc, char *argv[]){
    
    
    int i;
    double a, b, d, rm, sm, x0, I, k0, k1, k2, sgm, alpha, beta, M;
    double t, tf;
    //int sample=1;
    

    if(!(argc==3)){
        printf("ERROR: invalid argument usage:\n");
        printf("%s [k] [sigma] [total time]\n", argv[0]);
        return 0;
    }
    
    
    a=1;
    b=3;
    d=5;
    
    rm=0.001;
    sm=4;
    
    x0=-1.618;
    I=3.25;
    
    k0=atof(argv[1]);
    k1=1.0;
    k2=0.5;
    sgm=atof(argv[2]);
    
    alpha=0.1;
    beta=0.06;
    
    M=25;
    
    tf=atof(argv[3]);
    
    FILE *outf1;
    char str[1024];
    
    sprintf(str, "x_potential_k_%lf_sigma_%lf-tf-%.0lf-r.dat", k0, sgm, tf);
    outf1=fopen(str, "w");

    
    /* create a generator chosen by the
     environment variable GSL_RNG_TYPE */
    
    gsl_rng_env_setup();
    
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    
    seed = seed_dev();
    //seed = sample;
    gsl_rng_set(r, seed);
    
    
    /* gsl ode */
    const gsl_odeiv2_step_type * odeT;
    gsl_odeiv2_step * s;
    gsl_odeiv2_control * c;
    gsl_odeiv2_evolve * e;
    gsl_odeiv2_system sys;
    
    const int dimension = 4*Nosc; /* number of differential equations */
    int status;                /* status of driver function */
    const double eps_abs = 1.e-8;  /* absolute error requested  */
    const double eps_rel = 0; /* relative error requested  */
    
    /* initialise the ode solver */
    odeT = gsl_odeiv2_step_rkf45; /* solver: Runge-Kutta-Fehlberg (4, 5) method */
    s = gsl_odeiv2_step_alloc(odeT, dimension);
    c = gsl_odeiv2_control_y_new(eps_abs, eps_rel);
    e = gsl_odeiv2_evolve_alloc(dimension);
    
    double h;
    h = 1e-2; /* initial step-size */
    t = 0.0;
    
    double odeparams[14];
    
    odeparams[0] = a;       /* problem parameters */
    odeparams[1] = b;
    odeparams[2] = d;
    odeparams[3] = rm;
    odeparams[4] = sm;
    odeparams[5] = x0;
    odeparams[6] = I;
    odeparams[7] = k0;
    odeparams[8] = k1;
    odeparams[9] = k2;
    odeparams[10] = sgm;
    odeparams[11] = alpha;
    odeparams[12] = beta;
    odeparams[13] = M;

    
    x = (long double *)calloc(Nosc, sizeof(long double)); /* initial positions */
    y = (long double *)calloc(Nosc, sizeof(long double)); /* initial positions */
    z = (long double *)calloc(Nosc, sizeof(long double)); /* initial positions */
    phi = (long double *)calloc(Nosc, sizeof(long double)); /* initial positions */
    
    J = (long double *)calloc(Nosc, sizeof(long double)); /* initial positions */
    C = (long double *)calloc(Nosc, sizeof(long double)); /* initial positions */
    
    theta = (long double *)calloc(Nosc, sizeof(long double));

    
    double *v;
    v = (double *)calloc(4*Nosc, sizeof(double)); /* initial positions */
    
    for (i=0; i<Nosc; i++) {
    
        
        x[i]=gsl_rng_uniform(r);
        y[i]=gsl_rng_uniform(r);
        z[i]=gsl_rng_uniform(r);
        phi[i]=gsl_rng_uniform(r);

    }

    for (i=0; i<Nosc; i++) {
        v[4*i+0]=x[i];
        v[4*i+1]=y[i];
        v[4*i+2]=z[i];
        v[4*i+3]=phi[i];
    }
    
    sys = (gsl_odeiv2_system) {dfunc, NULL, dimension, &odeparams};

    
    while (t < tf ){
        
    
        status = gsl_odeiv2_evolve_apply(e, c, s, &sys, &t, tf, &h, v);
        if (status != GSL_SUCCESS)
        break;

        for (i=0; i<Nosc; i++){
            fprintf(outf1, "%lf %Lf\n", t, x[i]);
        }
        
    }
    
    
    fclose(outf1);

    gsl_rng_free (r);

    gsl_odeiv2_evolve_free(e);
    gsl_odeiv2_control_free(c);
    gsl_odeiv2_step_free(s);
    
    free(x);
    free(y);
    free(z);
    
    free(v);
    
        
    return 0;
    
}


int dfunc (double t, const double v[], double f[], void *params_ptr) {
    
    /* get parameter(s) from params_ptr */
    double *lparams = (double *) params_ptr;

    int i, j, k;
    double a, b, d, rm, sm, x0, I, k0, k1, k2, sgm, alpha, beta, M;
    double S1=0;

    
    a = lparams[0];
    b = lparams[1];
    d = lparams[2];
    
    rm = lparams[3];
    sm = lparams[4];
    
    x0 = lparams[5];
    
    I = lparams[6];
    
    k0 = lparams[7];
    k1 = lparams[8];
    k2 = lparams[9];
    sgm = lparams[10];
    
    alpha = lparams[11];
    beta = lparams[12];
    
    M = lparams[13];
    

    for (i=0; i<Nosc; i++) {
        x[i]=v[0+4*i];
        y[i]=v[1+4*i];
        z[i]=v[2+4*i];
        phi[i]=v[3+4*i];
    }
        

    for (i=0; i<Nosc; i++) {

        J[i]=-k0*x[i]*(alpha+3*beta*phi[i]*phi[i]);
        
        S1=0;
        j=fmod((i-M+Nosc),Nosc);

        for (k=0; k<2*M+1; k++) {
            
            S1+=x[j]-x[i];
            j=fmod((j+1),Nosc);
        }
        C[i]=(sgm/(2*M))*S1;
    
    }
    
    for (i=0; i<Nosc; i++) {
        f[0+4*i]=y[i]-a*pow(x[i],3)+b*x[i]*x[i]-z[i]+I+J[i]+C[i];
        f[1+4*i]=1-d*x[i]*x[i]-y[i];
        f[2+4*i]=rm*(sm*(x[i]-x0)-z[i]);
        f[3+4*i]=k1*x[i]-k2*phi[i];
    }
    
    
    return GSL_SUCCESS; /* GSL_SUCCESS in gsl/errno.h as 0 */
}




unsigned int seed_dev(void){
    
    unsigned int randval;
    FILE *arq;
    
    arq = fopen("/dev/random", "r");
    fread(&randval, sizeof(randval), 1, arq);
    fclose(arq);
    
    
    return randval;
}
