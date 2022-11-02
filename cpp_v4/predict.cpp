#include "predict.h"

void predict( double t0, double t1, double* y0, prms* ppc, 
              gsl_odeiv_step* s, gsl_odeiv_control* c, gsl_odeiv_evolve* e)
{
    int i,ac;
    gsl_odeiv_system sys = {func, jac, DIM, ppc};
    double t = t0;    
    double h = 1e-6;
    double ic[DIM];
    int counter=0;
    
    // copy initial conditions
    for(i=0;i<DIM;i++) ic[i]=y0[i];
    
    while (t < t1)
    {
        int status = gsl_odeiv_evolve_apply(e, c, s, &sys, &t, t1, &h, y0);

	//if(counter%40==0) printf("%1.1f \t %1.4f \t %1.4f\n", t, y0[0], y0[1] );
	
	counter++;
	
        if (status != GSL_SUCCESS)
            break;
    }
    
    return;
}



