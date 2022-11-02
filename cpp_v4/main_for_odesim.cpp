#include "essentials.h"
#include "predict.h"
#include "readdata.h"
#include "prms.h"

// global variables - CLO means it is a command-line option
double G_CLO_BETA1 = 1.20;
double G_CLO_BETA2 = 1.40;
double G_CLO_BETA3 = 1.60;

double G_CLO_SIGMA12 = 0.70;
double G_CLO_SIGMA13 = 0.30;
double G_CLO_SIGMA23 = 0.70;

double G_CLO_AMPL = 0.10;

double G_CLO_NU_DENOM = 5;
double G_CLO_RHO_DENOM = 900;

double G_CLO_EPIDUR = 365;

// this is the absolute population size threshold, for each subtype separately, at which the I-variable
// for that subtype will be set to 0.0
double G_CLO_ZEROPREVLEVEL = -1.0;
bool G_EXTINCTION_BEHAVIOR_ON[NUMSEROTYPES]; // this tells us whether the extinction behavior is currently on or off for each serotype

int G_CLO_NUMDAYS_BTW_IMMIGRATIONS = 1000000; // every G_CLO_NUMDAYS_BTW_IMMIGRATIONS days, we will immigrate a new serotype into the population

bool G_CLO_CHECKPOP_MODE = false;

//double G_CLO_OUTPUT_ALL_TRAJ = false;

double G_CLO_ADJUST_BURNIN = -100.0; // add this number of days to the t0 variable (currently set at -3000.0)
double G_CLO_ISS = 10.0; // this is the initial step size in the iterator

string G_CLO_STR_PARAMSFILE;
int G_CLO_INT_PARAMSFILE_INDEX = 0;



bool G_PHIS_INITIALIZED_ON_COMMAND_LINE = false;


// random number generator
gsl_rng *G_RNG;		



// FUNCTION DECLARATIONS


// this function parses all the command-line arguments
void ParseArgs(int argc, char **argv, prms* ppc);






int main(int argc, char* argv[])
{

    
    //NOTE this is commented out since there is no stochastic component right now
    //
    // get random number seed from current time
    //struct timeval pTV[1];
    //gettimeofday( pTV, NULL );
    //int seed = ((int)pTV->tv_sec) + ((int)pTV->tv_usec);  // this adds seconds to microseconds
    //
    // make random number generator (RNG) the Mersenne Twister which has period 2^19937 - 1
    //const gsl_rng_type *TT_RAND = gsl_rng_mt19937;
    //G_RNG = gsl_rng_alloc(TT_RAND);
    //gsl_rng_set( G_RNG, seed ); // seed the RNG        

    
    // must initialize these to false
    for(int vir=0; vir<NUMSEROTYPES; vir++) G_EXTINCTION_BEHAVIOR_ON[vir] = false;
    
    
    //
    // ###  1.  ALLOCATE SPACE FOR A PARAMETERS CLASS -- ppc=pointer to a parameters class
    //
    prms* ppc = new prms();  assert( ppc );
        
    // initialize the ODE system structures
    const gsl_odeiv_step_type* T = gsl_odeiv_step_rkf45;    // this tells you it's a 4th/5th ordre Runge-Kutta integration scheme
    ppc->os = gsl_odeiv_step_alloc(T, DIM);                 // this allocates memory for the various data structures
    ppc->oc = gsl_odeiv_control_y_new(1e-6, 0.0);           // this sets some parameters like error tolerance
    ppc->oe = gsl_odeiv_evolve_alloc(DIM);                  // this allocates some more memory (for some reason)


    
    
    //
    // ###  2.  PARSE THE COMMAND-LINE ARGUMENTS
    //

    ParseArgs( argc, argv, ppc );

    if (ppc->sigma[0][1] > 1) {
        fprintf(stderr,"\n\n\tWARNING : Sigma can't be over 1. You input %1.3f\n\n", ppc->sigma[0][1]); // %1.3f is a placeholder for what is being printed
        exit(-1);
    }
    

    //
    // ###  3.  INITIALIZE PARAMETERS - these are the default/starting values
    //
    
    // if the phi-parameters are not initialized on the command line
    if( !G_PHIS_INITIALIZED_ON_COMMAND_LINE )
    {
        for(int i=0;  i<10; i++) ppc->v[i] = ((double)i)*365.0 + 240.0; // sets the peak epidemic time in late August for the first 10 years
        for(int i=10; i<20; i++) ppc->v[i] = -99.0; 
    }
    
    ppc->v[ i_amp ]   = G_CLO_AMPL;
    ppc->beta[0] = G_CLO_BETA1 / POPSIZE_MAIN;    // NOTE this is in a density-dependent transmission scheme
    ppc->beta[1] = G_CLO_BETA2 / POPSIZE_MAIN;
    ppc->beta[2] = G_CLO_BETA3 / POPSIZE_MAIN;
    
    ppc->sigma[0][1] = G_CLO_SIGMA12; // 0.7; // the level of susceptibility to H1 if you've had B
    ppc->sigma[1][0] = G_CLO_SIGMA12; // 0.7; // and vice versa

    ppc->sigma[1][2] = G_CLO_SIGMA23; // 0.7; // the level of susceptibility to H3 if you've had B
    ppc->sigma[2][1] = G_CLO_SIGMA23; // 0.7; // and vice versa

    ppc->sigma[0][2] = G_CLO_SIGMA13; // 0.3; // the level of susceptibility to H3 if you've had H1
    ppc->sigma[2][0] = G_CLO_SIGMA13; // 0.3; // and vice versa
    
    ppc->sigma[0][0] = 0;
    ppc->sigma[1][1] = 0;
    ppc->sigma[2][2] = 0;
    
    ppc->v[ i_nu ]    = 1 / G_CLO_NU_DENOM;                // recovery rate
    ppc->v[ i_immune_duration ] = G_CLO_RHO_DENOM;    // 2.5 years of immunity to recent infection'
    
    ppc->v[ i_epidur ] = G_CLO_EPIDUR;
    
    //
    // ###  4.  SET INITIAL CONDITIONS FOR ODE SYSTEM
    //
    
    // declare the vector y that holds the values of all the state variables (at the current time)
    // below you are declaring a vector of size DIM
    double y[DIM];   

    
    //srand(time(NULL));
    
    for(int loc=0; loc<NUMLOC; loc++)
    {
        // put half of the individuals in the susceptible class
        y[ STARTS + loc ] = 0.5 * ppc->N[loc];
        
        // put small number (but slightly different amounts each time) of individuals into the infected classes
        // double r = rand() % 50 + 10;
        //double x = r / 1000.0; // double x = 0.010;
        double x = 0.010;
        double sumx = 0.0;
        for(int vir=0; vir<NUMSEROTYPES; vir++)
        {
            // double r = rand() % 50 + 1;
            // x = r / 1000.0;
            
            if (vir == 0) { x = 35 / 1000.0; }
            if (vir == 1) { x = 25 / 1000.0; }
            if (vir == 2) { x = 21 / 1000.0; }
            
            // fprintf(stderr, "r = %1.4f, x = %1.6f", r, x);
            
            sumx += x;
            y[ STARTI + NUMSEROTYPES*loc + vir ] = x * ppc->N[loc];
            y[ STARTJ + NUMSEROTYPES*loc + vir ] = 0.0;     // initialize all of the J-variables to zero
            
            x += 0.001;
        }
        x=0.010; // reset x
        
        // distribute the remainder of individuals into the different recovered stages equally
        for(int vir=0; vir<NUMSEROTYPES; vir++)
        {
            double z = (0.5 - sumx)/((double)NUMR*NUMSEROTYPES);  // this is the remaining fraction of individuals to be distributed 
            for(int stg=0; stg<NUMR; stg++)
            {
                y[ NUMSEROTYPES*NUMR*loc + NUMR*vir + stg ] = z * ppc->N[loc];
            }
        }
    }

    //
    // ###  5.  INTEGRATE ODEs AND OUTPUT STATE VARIABLES BY DAY
    //

    
    
    //
    //
    //BEGIN -- OUTPUT TRAJECTORIES TO STDOUT
    //
    //
    
    // start time and end time
    double t0=0.0; 
    double tf=NUMDAYSOUTPUT;
    
    int immigration_counter = 0;
    
    while( t0 < tf )
    {

        // print trajectory to stdout
        printf("%1.1f\t", t0); //fflush(stdout); // print time to stdout
        // printf("%1.5f\t", ppc->v[i_epidur]);
        printf("   %1.5f   \t", ppc->seasonal_transmission_factor(t0) );
        // printf("%1.2f \t", G_CLO_BETA1);
        for(int i=0; i<DIM; i++) printf("%1.1f\t", y[i]); //fflush(stdout);
	    printf("  %1.5f  \n", popsum(y) );
	    
        // integrate ODEs one day forward
        predict( t0, t0+1.0, y, ppc, ppc->os, ppc->oc, ppc->oe );
        
        
        // ###  5.1  ###  check if the I-variables are below their thresholds, and set them to zero if they are
        
        for(int loc=0; loc<NUMLOC; loc++)
        {
            for(int vir=0; vir<NUMSEROTYPES; vir++)
            {
                if( y[ STARTI + NUMSEROTYPES*loc + vir ] < G_CLO_ZEROPREVLEVEL  &&  G_EXTINCTION_BEHAVIOR_ON[vir] )
                {
                    y[ STARTS + loc ] += y[ STARTI + NUMSEROTYPES*loc + vir ];  // add a small number to the susceptibles
                    y[ STARTI + NUMSEROTYPES*loc + vir ] = 0.0;                 // substract this same small number from this serotype's I class
                    G_EXTINCTION_BEHAVIOR_ON[vir] = false;
                }
                
                if( y[ STARTI + NUMSEROTYPES*loc + vir ] > G_CLO_ZEROPREVLEVEL  &&  !G_EXTINCTION_BEHAVIOR_ON[vir] )
                {
                    G_EXTINCTION_BEHAVIOR_ON[vir] = true;
                }
            }
        }
        
        
        // ###  5.2  ###  immigrate an individual into the I-class
        
        int t_int = (int)t0;
        
        if( t_int%G_CLO_NUMDAYS_BTW_IMMIGRATIONS == 0 )
        {
            // this tells us which serotype will be immigrated
            int serotype_index = immigration_counter % NUMSEROTYPES;
            
            // loop through all locations
            for(int loc=0; loc<NUMLOC; loc++)
            {
                // add one person to the I class, and subtract one person from the S class
                if( y[ STARTS + loc ] > 1.1 ) // this should always be true
                {
                    y[ STARTI + NUMSEROTYPES*loc + serotype_index ] += 1.0;  
                    y[ STARTS + loc ] -= 1.0;
                }
                
            }  
            
            immigration_counter++;
        }
        
        
        
        
        

        // increment time by one day
        t0 += 1.0; 
    }
    // printf("\n\n\t %d \n\n", DIM);
    //
    //
    //END -- OUTPUT TRAJECTORIES TO STDOUT
    //
    //

    
    
    
    
    
    // free memory
    gsl_odeiv_evolve_free( ppc->oe );
    gsl_odeiv_control_free( ppc->oc );
    gsl_odeiv_step_free( ppc->os );
    

    //delete[] y0;
    //delete[] incidence;
    delete ppc;
    //delete[] pContactMatrix;
    return 0;
}


bool isFloat( string myString ) {
    std::istringstream iss(myString);
    float f;
    iss >> noskipws >> f; // noskipws considers leading whitespace invalid
    // Check the entire string was consumed and if either failbit or badbit is set
    return iss.eof() && !iss.fail(); 
}


// parses command line arguments
void ParseArgs(int argc, char **argv, prms* ppc)
{
    string str;
    int i, start;
    i=1;    // this is the position where you start reading command line options
            // you skip "i=0" which is the text string "odesim"

    /*if( argc<start )
    { 
        PrintUsageModes(); 
        exit(-1);
    }*/

    
    // read in options from left to right
    while(i<argc)
    {
        str = argv[i]; // read the ith text string into the variable "str"
        
        //BEGIN MAIN IF-BLOCK BELOW
        
        // ### 1 ### IF BLOCK FOR PHI
        if( str == "-phi" )
        {
            ppc->phis.clear();
            i++;
            
            //BEGIN LOOPING THROUGH THE LIST OF BETAS
            while(i<argc)
            {
                string s( argv[i] );    // convert argv[i] into a normal string object
                if( isFloat(s) )        // if the current string is a floating point number, write it into the phis array
                {
                    // if the command line argument is <0, just set it back to zero
                    double d = atof( argv[i] );
                    if( d < 0.0 ){
                        d = 0.0;
                        //TODO print warning here and probably should exit
                        fprintf(stderr, "\n\n \t Don't make phis less than zero! \n\n");
                    }
                    
                    ppc->phis.push_back( d );
                                        
                    // increment and move on in this sub-loop
                    i++;
                }
                else
                {
                    // if the current string is NOT a float, set the counter back down (so the string can be handled by the next if block
                    // and break out of the loop
                    i--;
                    break;
                }
                 
            } 
            //END OF LOOPING THROUGH THE LIST OF BETAS

            // make sure at least one phi-value was read in
            if( ppc->phis.size() == 0 )
            {
                fprintf(stderr,"\n\n\tWARNING : No phi-values were read in after the command-line option \"-phi\".\n\n");
            }
        }
        // ### 2 ### BLOCKS FOR FOR THE OTHER NON-PHI COMMAND-LINE OPTIONS
        else if( str == "-checkpop" )
        {
            G_CLO_CHECKPOP_MODE = true;
        }
        else if( str == "-beta1" )              {     G_CLO_BETA1 = atof( argv[++i] );        }
        else if( str == "-beta2" )              {     G_CLO_BETA2 = atof( argv[++i] );        } //atof is character to floating point
        else if( str == "-beta3" )              {     G_CLO_BETA3 = atof( argv[++i] );        } //atoi changes it to integer
        else if( str == "-sigma12" )            {   G_CLO_SIGMA12 = atof( argv[++i] );      } //atoi changes it to integer
        else if( str == "-sigma13" )            {   G_CLO_SIGMA13 = atof( argv[++i] );      }
        else if( str == "-sigma23" )            {   G_CLO_SIGMA23 = atof( argv[++i] );      }
        else if( str == "-amp" )                {   G_CLO_AMPL          = atof( argv[++i] );         }
        else if( str == "-nu_denom")            {   G_CLO_NU_DENOM      = atof( argv[++i] );     }
        else if( str == "-rho_denom")           {   G_CLO_RHO_DENOM     = atof( argv[++i]);      }
        else if( str == "-epidur")              {   G_CLO_EPIDUR        = atof( argv[++i]);        }
        else if( str == "-zeroprevlevel")       {   G_CLO_ZEROPREVLEVEL = atof( argv[++i]);   }
        else if( str == "-immig")               {   G_CLO_NUMDAYS_BTW_IMMIGRATIONS     = atoi( argv[++i]);   } // this tells us that we will immigrate a new 
                                                                                                // serotype into the population every G_CLO_NUMDAYS_BTW_IMMIGRATIONS days 
        

        else
        {
            fprintf(stderr, "\n\tUnknown option [%s] on command line.\n\n", argv[i]);
            exit(-1);
        }
        //END OF MAIN WHILE-LOOP BLOCK; INCREMENT AND MOVE ON TO THE NEXT COMMAND-LINE ARGUMENT
 
        // increment i so we can look at the next command-line option
        i++;
    }

    return;
}





