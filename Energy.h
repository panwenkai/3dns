//  _________________________________________________________
// |
// |   Energy.h    
// |
// |
// |   (C) 1998-99  Columbia University, MSME
// |   (C) 2000-01	Harvard University, DEAS
// |__________________________________________________________
#ifndef _LASERH
#define _LASERH

//_________________________________________________
// Public functions
//_________________________________________________
void LaserInit ();
void LaserInput ();
void LaserCleanup ();

//_________________________________________________
// Public variables
//_________________________________________________

#ifdef EXT_LEVEL
#define EXTERN
#else
#define EXTERN extern
#endif
EXTERN DMATRIX QLaser;		// laser energy deposited in each node
#undef EXTERN

#endif