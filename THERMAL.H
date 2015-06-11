//  _________________________________________________________
// |
// |   Thermal.h    Thermal conductivity header
// |
// |
// |   (C) 1998-99  Columbia University, MSME
// |   (C) 2000-01	Harvard University, DEAS
// |__________________________________________________________
#ifndef _THERMALH
#define _THERMALH

#include <stdarg.h>		       //for va_list etc

#define SISOLID 0
#define SILIQUID 1
#define SILICA 2
#define SUBSTRATE 3

#define K_MIN ((double)0.0)    //minimum thermal conductivity W/cm3K
#define K_MAX ((double)100.0)  //maximum thermal conductivity W/cm3K
#define K_ISVACUUM ((double)1E-12)  //threshold below which is considered vacuum

#define CP_MIN ((double)0.0)	//minimum heat capacity J/cm3K
#define CP_MAX ((double)100.0)  //maximum heat capacity J/cm3K

//_________________________________________
// Public Functions
//_________________________________________
int HeatFlow (void);		//Redistributes heat
void HeatFlowDump();		//dumps calculations from current clock
void HeatFlowFinalize();	//finalizes changes to T matrix
//inline double Interpolate2D (double t1, double t2, double t3, double t4, double fracI, double fracJ);
double InterpolateT (int i, int j, int k, double iFrac, double jFrac, double kFrac);
double InterpolateT(const int nodeA, const int dirEdge);
void StabilityCheck();		//Check that temperatures are stable
void ThermalInit (void);	//Initializes arrays
void ThermalCleanup (void);	//Frees arrays

//_________________________________________________
// Public variables
//_________________________________________________


#endif

