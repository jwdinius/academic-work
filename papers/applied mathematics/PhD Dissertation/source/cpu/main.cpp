//
//  FILE:   main.cpp
//  MODULE: diskDynamics
//
//  DESCRIPTION:
//  File contains the main executive for running generalized collision simulations.  Input is read from a file, time integration and collision updates of both phase space and tangent space vectors is performed, output is recorded to user-defined files and allocated memory at run-time is freed.
//
//  REFERENCE:
//  Dinius, J.  "Dynamical Properties of a Generalized Collision Rule for Multi-Particle Systems" Ph.D dissertation.
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

// INCLUDES/DECLARATIONS
#include <iostream>
#include <string>
#include "utilities.h"
#include <cmath>
#include "stdlib.h"

using namespace std;
// END INCLUDES/DECLARATIONS

// VARIABLES IN GLOBAL SCOPE
// see "utilities.h" for variable description
/*int                      DIM, nlya, nDisks, phaseDim, NOCOLL, maxFlight, iSeed, i_coll, nOrtho, countCLV, modCLV, lastPass;
double                   *y, *maxFlightDblArray, *cum;
double                   density, stepSize, bigTime, Time, TimeStartCollection, vq, vv, beta, alpha, a_alpha, b_alpha, c_alpha, clvSettleTime, boxSize;
COLL                     *collArray;
ofstream                 clvData;
vector< vector<double> > Rinv, GS, traj;
int                      processCLVs;

// variables declared to avoid memory allocation issues
double *dq1, *lSpec, *norm, *yi, *yj, *vi, *vj, *dqc, *dv, *dq, *v, *qi;
*/
// END VARIABLES IN GLOBAL SCOPE

//
//  FUNCTION: main
//  MODULE:   diskDynamics
//
//  DESCRIPTION:
//  Function is the main executive routine for generalized collision simulations.
//
//  INPUTS:
//  argc - standard form for C/C++ main routine for integer input (currently unused)
//  argv - standard form for C/C++ main routine for character array input (used for input filename to read initialization data from)
//
//  OUTPUTS:
//  0 - normal execution
//  1 - clean exit due to user-defined exception declared
//
//  GLOBAL INPUTS (UNCHANGED DURING CALL):
//  (none)
//  GLOBAL OUTPUTS (UPDATED DURING CALL):
//  stepSize         (value declared)
//  nOrtho           (value declared)
//  nDisks           (value declared)
//  density          (value declared)
//  nlya             (value declared)
//  DIM              (value declared)
//  phaseDim         (value declared)
//  iSeed            (initialized)
//  NOCOLL           (value declared)
//  maxFlight        (value declared)
//  bigTime          (value declared)
//  i_coll           (initialized)
//  Time             (initialized and incremented)
//  y                (memory allocated and specified states written to file)
//  lSpec            (memory allocated)
//  collArray        (memory allocated)
//  boxSize          (memory allocated and values declared)
//  maxFlightDblArray (memory allocated and values declared)
//  norm             (memory allocated)
//  cumBatch         (memory allocated)
//  avgCumBatch      (memory allocated)
//  qi               (memory allocated)
//  yi               (memory allocated)
//  yj               (memory allocated)
//  vi               (memory allocated)
//  vj               (memory allocated)
//  v                (memory allocated)
//  dq               (memory allocated)
//  dv               (memory allocated)
//  dqc              (memory allocated)
//  modCLV           (initialized)
//  lastPass         (initialized to false)
//  clvData          (value declared, file opened for writing)
//
//  REVISION HISTORY:
//  Dinius, J.       Created                              10/08/11
//  Dinius, J.       Comments added for v1.0 release      05/27/13   
//  Dinius, J.       Cleaned up for dissertation release  11/25/13
//  Dinius, J.       Streamlined for gpu comparison       03/01/14
//
int main (int argc, const char * argv[])
{
    // DECLARATIONS FOR INPUT
    ifstream inData;          // file handle to input file
    string   dummyStr;        // dummy string for holding temporary data read-in from input file
    string   str1 ("*/");     // string for end of line comparison during input read
    int      nSteps;          // number of time steps to run simulation
    int      nModes;          // number of tangent vectors to write out to file
    int      nOutput;         // number of steps between successive writes to main output file
    int      i,ii,j;          // loop counters (initialized once for convenience)
    // END DECLARATIONS FOR INPUT

	//OPEN INPUT FILE
	inData.open(argv[1]);
    
	if (!inData){
		// INPUT FILE COULD NOT BE OPENED, SO EXIT
		cerr << "Error: input file could not be opened.  Exiting..." << endl;
		exit(1);
	} // END if (!inData)
	
	// DECLARATIONS FOR OUTPUT
	ofstream outData;                // file handle to output file
	char *outFile = new char [256];  // character array to hold output filename
	// END DECLARATIONS FOR OUTPUT

	// READ IN DATA FROM INPUT FILE AND OUTPUT VALUES TO SCREEN
	cout << "mdcpu-serial: Read-in parameters:" << endl;
	inData >> stepSize;
	inData >> dummyStr;
	cout << "Stepsize = " << stepSize << endl;
	while (dummyStr.compare(str1) != 0){
		inData >> dummyStr;
	} // END while (dummyStr.compare(str1) != 0)

	inData >> nSteps;
	inData >> dummyStr;
	cout << "No. of Steps = " << nSteps << endl;
	while (dummyStr.compare(str1) != 0){
		inData >> dummyStr;
	} // END while (dummyStr.compare(str1) != 0)

	inData >> nOrtho;
	inData >> dummyStr;
	cout << "No. of Steps between Orthonormalization Operations = " << nOrtho << endl;
	while (dummyStr.compare(str1) != 0){
		inData >> dummyStr;
	} // END while (dummyStr.compare(str1) != 0)

	inData >> nOutput;
	inData >> dummyStr;
	cout << "No. of Steps between File Output = " << nOutput << endl;
	while (dummyStr.compare(str1) != 0){
		inData >> dummyStr;
	} // END while (dummyStr.compare(str1) != 0)

	inData >> nDisks;
	inData >> dummyStr;
	cout << "No. of Disks = " << nDisks << endl;
	while (dummyStr.compare(str1) != 0){
		inData >> dummyStr;
	} // END while (dummyStr.compare(str1) != 0)

	inData >> density;
	inData >> dummyStr;
	cout << "Density of disks in box = " << density << endl;
	while (dummyStr.compare(str1) != 0){
		inData >> dummyStr;
	} // END while (dummyStr.compare(str1) != 0)
	
	inData >> outFile;
	outData.open(outFile,ios::out | ios::binary);
	cout << "File to record output: " << outFile << endl;
	inData >> dummyStr;
	while (dummyStr.compare(str1) != 0){
		inData >> dummyStr;
	} // END while (dummyStr.compare(str1) != 0)

	inData >> processCLVs;
	cout << "Compute covariant vectors? " << processCLVs << endl;
	inData >> dummyStr;
	while (dummyStr.compare(str1) != 0){
		inData >> dummyStr;
	} // END while (dummyStr.compare(str1) != 0)

	if (processCLVs){
		inData >> clvSettleTime;
		cout << "Time to run to settle GS transients? " << clvSettleTime << endl;
		inData >> dummyStr;
		while (dummyStr.compare(str1) != 0){
			inData >> dummyStr;
		} // END while (dummyStr.compare(str1) != 0)
	}
	else {
		// IF NOT PROCESSING CLVs, SET THE SETTLE TIME TO SOMETHING VERY LARGE
		clvSettleTime = 1.0e9;
	}

	// CLOSE THE INPUT FILE
	inData.close();
	// END READ IN DATA FROM INPUT FILE AND OUTPUT VALUES TO SCREEN

	if (processCLVs){
		string   str2 (outFile);
		string   clvFile("clv_");
		clvFile.append(str2);
		clvData.open(clvFile,ios::out | ios::binary);
	} // END if (processCLVs)

	// SET PROGRAM CONSTANTS/DEFAULTS
	DIM       = 2;
	phaseDim  = 2*DIM*nDisks;
	nlya      = phaseDim;
	iSeed     = 112324;
	NOCOLL    = -1;
	maxFlight = -2;
	bigTime   = 1.0e11;
	double pi = 4.0*atan(1.0);
	modCLV  = 25;
	lastPass = 0;

	// INITIALIZE COUNTERS AND TIME
	i_coll      = 0;
	Time        = 0.0;
	
	// ALLOCATE MEMORY FOR GLOBAL VARIABLES
	y                 = new double [(nlya+1)*phaseDim];    
	collArray         = new COLL [nDisks];
	maxFlightDblArray = new double [DIM];
	cum               = new double [nlya];
	dq1               = new double [DIM];
	lSpec             = new double [nlya];
	norm              = new double [nlya];
	qi                = new double [phaseDim];
	yi                = new double [DIM];
	yj                = new double [DIM];
	vi                = new double [DIM];
	vj                = new double [DIM];
	v                 = new double [DIM];
	dq                = new double [DIM];
	dv                = new double [DIM];
	dqc               = new double [DIM];

	// SET SIZE OF SIMULATION BOX
	// SOLVES N*pi*R^2 = density * A == density*(aspRatio*L^2) FOR L, DISKS HAVE RADIUS 0.5 BY CONSTRUCTION
	boxSize = sqrt( ((double) nDisks) * pi / 4.0 / density );
	// SET SIZE OF MAXIMUM FLIGHT CONDITION
	// THE MAXIMUM FLIGHT CONDITION IS, ARBITRARILY, SET TO 1/4 THE BOX LENGTH IN EACH DIMENSION (MINUS THE DISK RADIUS)
	for (int i = 0; i < DIM; i++){
		maxFlightDblArray[i] = (boxSize / 2.0 - 1.0) / 2.0;
	} // END for (int i = 0; i < DIM; i++)

	// CALL INITIALIZE ROUTINE TO POPULATE PHASE SPACE AND TANGENT VECTORS, ALONG WITH COLLISION TIMES
	cout << "mdcpu-serial: initialize()" << endl;
	initialize();
	// END INITIALIZE
	
	/////////////////////////////////////////////////////
	//              MAIN LOOP                          //
	/////////////////////////////////////////////////////
	cout << "mdcpu-serial: update()" << endl;
	for (i = 0; i <= nSteps; i++){
	
		if ( (i % nOrtho) == 0 ){
			// PERFORM GRAM-SCHMIDT REORTHONORMALIZATION OF TANGENT VECTORS
#ifdef LAPACK
#elif  QR
#else
			mgsr();
#endif
			// UPDATE COMPUTATION OF STRETCHING FACTORS AND LYAPUNOV EXPONENTS
			update();

		} // END if ( (i % nOrtho) == 0 )

		if ( (i % nOutput) == 0 ){
			// OUTPUT CURRENT TIME AND LARGEST LYAPUNOV EXPONENT VALUES TO THE SCREEN AS A STATUS CHECK
			cout << "\t" << Time << "\t" <<lSpec[0] << endl;

			// OUTPUT CURRENT TIME TO FILE
			outData.write((char *) &Time, sizeof(double));

			// OUTPUT POSITION VECTOR AT CURRENT TIME TO FILE (SEQUENTIALLY BY VECTOR INDEX)
			for (j = 0; j < phaseDim; j++){
				outData.write((char *) &y[j], sizeof(double));
			} // END for (j = 0; j < phaseDim; j++)

			
			// OUTPUT LYAPUNOV SPECTRUM AT CURRENT TIME
			for (ii = 0; ii < nlya; ii++){
				outData.write((char *) &lSpec[ii], sizeof(double));
			} // END for (ii = 0; ii < nlya; ii++)

			// OUTPUT NUMBER OF COLLISIONS
			outData.write((char *) &i_coll, sizeof(int));

		} // END if ( (i % nOutput) == 0 )

		// PERFORM TIME STEP PROPAGATION OF SYSTEM
		hardStep(stepSize);

		// INCREMENT TIME
		Time += stepSize;
		
		if (checkOverlap()){
			// IF DISKS OVERLAP, EXIT
			cerr << "Disks overlap... exiting" << endl;
		} // END if (checkOverlap())

	} // END for (i = 0; i < nSteps; i++)

	// FINALIZE (OUTPUT END-OF-RUN STATS AND PERFORM CLEANUP)
	cout << "mdcpu-serial: finalize()" << endl;

	// OUTPUT MEAN FREE TIME (AVERAGE TIME BETWEEN SUCCESSIVE COLLISIONS FOR AN INDIVIDUAL DISK)
	cout << "Mean free time = " << (double)(nDisks)*Time / ((double)(i_coll)) << endl;

	// CLOSE MAIN OUTPUT FILE
	outData.close();

	// IF INDICATED FROM INPUT FILE, COMPUTE COVARIANT LYAPUNOV VECTORS
	// USING GINELLI'S METHOD
	if (processCLVs){
		computeCLVs();
	}

	// CLEAR ALLOCATED MEMORY
	delete [] y;
	delete [] collArray;
	delete [] maxFlightDblArray;
	delete [] cum;
	delete [] norm;
	delete [] dq1;
	delete [] lSpec;
	delete [] yi;
	delete [] yj;
	delete [] vi;
	delete [] vj;
	delete [] v;
	delete [] dqc;
	delete [] dq;
	delete [] dv;
	delete [] qi;
	delete [] outFile;
	nModes = 0;
	
    // PROGRAM RAN TO COMPLETION, EXIT NORMALLY
    return 0;
    
} // END int main (int argc, const char * argv[])

