#include "essentials.h"
#include "prms.h"

// this function integrates from time t0 to time t1
void predict( double t0, double t1, double* y0, prms* ppc, 
              gsl_odeiv_step* s, gsl_odeiv_control* c, gsl_odeiv_evolve* e);
 




