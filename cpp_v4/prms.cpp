//#include <iostream>
//#include <string>
//#include <cstdlib>

#include "assert.h"
#include "prms.h"

// constructor
prms::prms( )
{
    v.insert( v.begin(), num_params, 0.0 );
    assert( v.size()==num_params );
    
    //dim = num_params;
    
    for(int i=0; i<NUMSEROTYPES; i++)
    {
        for(int j=0; j<NUMSEROTYPES; j++)
        {
            sigma[i][j] = 0.0; // this initializes the system with full cross-immunity among serotypes
        }
    }

    for(int i=0; i<NUMLOC; i++)
    {
        for(int j=0; j<NUMLOC; j++)
        {
            if(i==j) eta[i][j] = 1.0;  // the diagonal elements have to be 1.0, as a population mixes fully with itself
            if(i!=j) eta[i][j] = 0.0;  // these are initialized to 0.0 indicating that the different sub-populations do not mix at all
        }
    }
    
    // set initial values for the population sizes, 1 million for first location and 100K for others
    N[0] = POPSIZE_MAIN;
    for(int i=1; i<NUMLOC; i++) N[i] = POPSIZE_OUT;
    
    os=NULL;
    oc=NULL;
    oe=NULL;
}

// destructor
prms::~prms()
{
}


double prms::seasonal_transmission_factor( double t )
{
    
    /*
        modulus (denoted with %) only works with integers, so need the code below to
            create ten-year model cycles
     
     */
    
    // This is needed to create the 10-year "cycles" in transmission.
     
    int x = (int)t; // This is now to turn a double into an integer
    double remainder = t - (double)x;
    int xx = x % 3650; // int xx = x % NUMDAYSOUTPUT;
    double yy = (double)xx + remainder;
    // put yy into the sine function, let it return the beta value
    t = yy;
    
    
    
    
    
    
    double sine_function_value = 0.0;
    
    for(int i=0; i<phis.size(); i++)
    {
        // if( fabs( t - phis[i] ) < 91.25 )
        if( fabs( t - phis[i] ) < (v[i_epidur] / 2))
        {
            // sine_function_value = sin( 2.0 * 3.141592653589793238 * (phis[i]-t+91.25) / 365.0 );
            sine_function_value = sin( 2.0 * 3.141592653589793238 * (phis[i]-t+(v[i_epidur] / 2)) / (v[i_epidur] * 2));
            // printf("\n\t\t\t %1.3f %1.3f %1.3f \n\n", t, phis[i], sine_function_value );
        }
    }

    return 1.0 + v[i_amp] * sine_function_value;    
}





