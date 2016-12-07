//
//  FILE:   utilitiesTypes.h
//  MODULE: diskDynamics
//
//  DESCRIPTION:
//  File contains basic class definitions for running generalized collision simulations.
//
//  REVISION HISTORY:
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
//

// BEGIN INCLUDES/DEFINITIONS
#ifndef diskDynamics_utilitiesTypes_h
#define diskDynamics_utilitiesTypes_h
// END INCLUDES/DEFINITIONS

// BEGIN CLASS DEFINITIONS:

// CLASS COLL DEFINITION
class COLL {
public:
    double time;    // next collision time of current disk
    int    index;   // current disk index
    int    partner; // collision partner of current disk at next collision
};

// END CLASS DEFINITIONS
#endif
