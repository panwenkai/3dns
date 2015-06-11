//  _________________________________________________________
// |
// |   Error.h
// |
// |   (C) 1998-99  Columbia University, MSME
// |   (C) 2000-01	Harvard University, DEAS
// |__________________________________________________________
#ifndef _ERRORH
#define _ERRORH

//____ execution modes _____________
#define INFO_OFF		((int) 0)	
#define INFO_ON			((int) 1)	//InfoMsg to console & log

#define ERR_NONE		((int) 0)	//Supress errors
#define ERR_CTRLC		((int) 1)	//ErrorMsg wait for keypress
#define ERR_BATCH		((int) 2)	//ErrorMsg to console & log then return
#define ERR_MATHLINK	((int) 3)	//MathLink error reporting
	
//_________________________________________
// Public functions
//_________________________________________
void ErrorMsg (char *fmt,...);
void InfoMsg (char *fmt,...);
bool InRange (int param, int iMin, int iMax);
bool InRange (double param, double vMin, double vMax);
bool InRange (double *array, int iStart, int iStop, double vMin, double vMax);
bool InRange (int *array, int iStart, int iStop, int vMin, int vMax);
void SetErrMode(int mErr);
void SetInfoMode(int mInfo);
void WarningMsg (char *fmt,...);
void ZeroStructures ();

//_________________________________________
// Public variables
//_________________________________________

#endif