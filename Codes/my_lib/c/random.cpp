#ifndef randomguard
#define randomguard

//=========================================================================================================
//Random library:
#include "./randomc/randomc.h"

#ifndef MULTIFILE_PROJECT
// If compiled as a single file then include these cpp files, 
// If compiled as a project then compile and link in these cpp files.
   // Include code for the chosen random number generator:
   #include "./randomc/mersenne.cpp"
   // define system specific user interface:
   #include "./randomc/userintf.cpp"
#endif

int seed = (int)time(0); // random seed
CRandomMersenne Mersenne (seed);

double randX(double const & min = 0., double const & max = 1.0)
{	
	return (min + Mersenne.Random()*(max-min));	
};

int randIX(int const & min= 0, int const & max= 1)
{
	return (Mersenne.IRandomX(min,max));
};

//===========================================================================================================
#endif