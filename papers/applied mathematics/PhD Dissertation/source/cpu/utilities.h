//
//  FILE:   utilities.h
//  MODULE: diskDynamics
//
//  DESCRIPTION:
//  File contains global variable and function declarations for running generalized collision simulations.
//
//  REVISION HISTORY:
//  Dinius, J.       Created                              10/08/11
//  Dinius, J.       Added batch averaging and tumbling   03/15/12
//  Dinius, J.       Added generalized collision rule     02/18/13
//  Dinius, J.       Added localization width and average 04/01/13 
//                   angle calculations
//  Dinius, J.       Added velocity-velocity              04/09/13
//                   autocorrelation calculation
//  Dinius, J.       Comments added for v1.0 release      05/27/13
//  Dinius, J.       Added capability for computing       11/02/13 
//                   covariant Lyapunov vectors.
//  Dinius, J.       Cleaned up for dissertation release  11/25/13
//

// BEGIN INCLUDES/DEFINITIONS
#ifndef diskDynamics_utilities_h
#define diskDynamics_utilities_h
#include "utilitiesTypes.h"
#include <fstream>
#include <vector>
#include <random>

using namespace std;

// END INCLUDES/DEFINITIONS

/* BEGIN GLOBAL VARIABLE DECLARATIONS
   NOTE ON GLOBAL DECLARATIONS: The statement (variable made global for memory allocation convenience) is made to indicate a global variable that is declared for convenient memory allocation ONLY.  These variables need not be global, as they are not used in multiple functions, however they are made global to prevent the re-allocation of memory each time functions are called that use these variables */

/*extern int DIM;                   // spatial dimension of the configuration space 

extern int nlya;                  // number of Lyapunov exponents to compute

extern int nDisks;                // number of disks

extern int phaseDim;              // dimension of total phase space (=2*DIM*nDisks)

extern int NOCOLL;                // collision partner index when no collision is possible (cannot be an integer 
                                  // in [0,nDisks) )

extern int maxFlight;             // collision partner index when the disk would travel a distance defined by a 
                                  // distance constraint before any possible collisions (cannot be an integer in 
                                  // [0,nDisks) )

extern int iSeed;                 // value of random seed to initialize random number generator

extern int i_coll;                // counter for number of collisions that have occurred

extern int nOrtho;                // number of simulation steps in between successive applications of the Gram-
                                  // Schmidt orthnormalization process to tangent vectors

extern int countCLV;              // counter for storing CLV vectors

extern int modCLV;                // counter for storing CLV vectors

extern double *y;                 // array containing all of the state values (indices 0-(4N-1) ) and values for 
                                  // the nlya tangent vectors (i = 1,nlya: indices (i*(4N)+0-(i+1)*(4N)-1) )

extern double boxSize;            // size of the simulation box (assumes square box)

extern double *maxFlightDblArray; // distance condition for declaring maxFlight as a disk's collision partner

extern double *cum;               // vector storing accumulated logarithms of tangent vector stretching factors

extern double density;            // volume density (packing fraction) of disks in the simulation box

extern double stepSize;           // simulation time step

extern COLL   *collArray;         // array of structures storing own index, collision partner index and  
                                  // associated collision time for each disk (relative to current simulation time
                                  // )

extern double bigTime;            // (large) time to reinitialize disk next collision time after a collision 
                                  // involving the disk has occurred

extern double Time;               // simulation time

extern double vq;                  // dot product of relative momentum vector with relative position vector of 
                                   // disks undergoing collision

extern double vv;                  // dot product of relative momentum vector with itself of disks undergoing 
                                   // collision

extern double *dq1;                // relative position vector connecting centers of disks undergoing collision

extern double *lSpec;              // spectrum of Lyapunov exponents

extern double *norm;               // tangent vector stretching factors

extern double *yi;                 // position vector of ith disk (variable made global for memory allocation 
                                   // convenience)

extern double *yj;                 // position vector of jth disk (variable made global for memory allocation 
                                   // convenience)

extern double *vi;                 // momentum vector of ith disk (variable made global for memory allocation 
                                   // convenience)

extern double *vj;                 // momentum vector of jth disk (variable made global for memory allocation 
                                   // convenience)

extern double *v;                  // relative momentum vector between disks undergoing colliskion (variable made          
                                   // global for memory allocation convenience)

extern double *dq;                 // relative tangent position vector between disks undergoing collision 
                                   // (variable made global for memory allocation convenience)

extern double *dv;                 // relative tangent momentum vector between disks undergoing collision   
                                   // (variable made global for memory allocation convenience)

extern double *dqc;                // difference in relative position vectors between collision points in 
                                   // reference (at q) and perturbed (at q+dqc) trajectories (variable made 
                                   // global for memory allocation convenience)

extern double *qi;                 // orthogonal projector column vector for modified Gram-Schmidt process 
                                   // (variable made global for memory allocation convenience)

extern double *mean;               // computed mean of disk momenta after initialization by random draw (variable 
                                   // made global for memory allocation convenience)
extern int lastPass;               //

extern int processCLVs;            // flag to determine whether or not to compute covariant vectors (CLVs)

extern ofstream clvData;                    // file handle for covariant vector output file

extern vector< vector<double> > Rinv;       // vector holding time history of Gram-Schmidt scaling factors for
                                            // CLV computation

extern vector< vector<double> > GS;         // vector holding time history of Gram-Schmidt orthonormal bases
                                            // for CLV computation

extern vector< vector<double> > traj;       // vector holding time history of phase-space trajectories
                                            // for CLV computation

extern double clvSettleTime;                // amount of time before before beginning Ginelli method computations
*/
/* END GLOBAL VARIABLE DECLARATIONS */

/* BEGIN FUNCTION DECLARATIONS */
void      hardStep(double t);
COLL      nextColl(void);
void      freeFlight(double dt);
void      boxSet(void);
void      updateTimes(double tLast, int d1, int d2);
double    binTime(int , int );
double    image(double y1,double y2, double boxSize, int dim, int DIM);
void      collision(int , int );
int       checkOverlap(void);
int       overlap(int , int );
void      initialize(void);
double    randn(void);
double    ranf(void);
void      mgsr(void);
void      update(void);
/*vector<double> invertRmat(vector<double> Rvec);
void      computeCLVs(void);
vector<double> matrixMultiply(vector<double> v1, vector<double> v2);
double columnNorm(vector<double> v, int col);
vector<double> rescaleCol(vector<double> v, int col, double nrm);
void initializeCLVcollTimes(void);
vector<double> rowToColumn(vector<double> V);*/
/* END FUNCTION DECLARATIONS */

#endif
