#include <math.h>
#include <stdio.h>
#include <ctype.h>
#include <assert.h>
#include <stdlib.h>
#include <ctime>
#include <sys/time.h>

#include <gsl/gsl_rng.h> // random number generators from Gnu Scientific Library
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

#include <gsl/gsl_sf_gamma.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>

#include <map>
#include <vector>
#include <algorithm>


// number of locations 
#define NUMLOC 1

// number OF types/subtypes of influenza (this will always be three - H1, H3, and B) 
// for generality (and to avoid constantly having to specify type/subtype) we call this serotypes
#define NUMSEROTYPES 3

// number of R stages in the recovery class
#define NUMR 6

// the start index of the infected classes; the I-classes start after all of the R-classes
// have been listed
#define STARTI NUMLOC*NUMSEROTYPES*NUMR // 12

// the start index of the cumulative infected classes, i.e. the J-classes
#define STARTJ STARTI + NUMLOC*NUMSEROTYPES // 15

// the start index of the suceptible (S) classes; there will be one of these for every location
#define STARTS STARTJ + NUMLOC*NUMSEROTYPES // 18

// this is the dimensionality of the ODE system
#define DIM STARTS+NUMLOC // 19

// this is the number of days of simulation that will be sent to standard output (and used for model fitting)
#define NUMDAYSOUTPUT 3650*2 // use this to define "cycle" lengths

// Two population sizes: main population and all outer populations have the same size
#define POPSIZE_MAIN 1000000.00
#define POPSIZE_OUT 100000.00

//#define HAVE_INLINE 	// this will make GSL vector operations faster




#ifndef ESSENTIALS
#define ESSENTIALS

using namespace std;

//
//
// this function contains the ode system
//
int func(double t, const double y[], double f[], void *params);



// void* jac;	// do this for C-compilation
//
// for C++ compilation we are replacing the void* declaration above with
// the inline dummy declaration below
inline int jac(double a1, const double* a2, double* a3, double* a4, void* a5)
{
    return 0;	
};

inline double popsum( double yy[] )
{
    double sum=0.0;
    for(int i=0; i<DIM; i++) sum += yy[i];
    
    for(int i=STARTJ; i<STARTJ+NUMLOC*NUMSEROTYPES; i++) sum -= yy[i];
    
    return sum;
}



inline void write_to_file( const char* szFilename, vector< vector<double> >& vvDATA )
{
    FILE* fp = fopen( szFilename, "w" );
    int nr = vvDATA.size();	// number of rows
    int nc = vvDATA[0].size();
    
   for(int rr=0;rr<nr;rr++)
    {
        for(int cc=0;cc<nc;cc++)
        {
            fprintf(fp, "%1.3f \t", vvDATA[rr][cc] );
        }
        fprintf(fp, "\n");
    }

    fclose(fp);
    return;
}




#endif // ESSENTIALS
