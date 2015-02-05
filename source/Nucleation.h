//  _________________________________________________________
// |
// |   Nucleation.h
// |
// |
// |   (C) 1998-99  Columbia University, MSME
// |   (C) 2000-01	Harvard University, DEAS
// |__________________________________________________________
#ifndef _NUCLEATIONH
#define _NUCLEATIONH

#define HET_RATE_MAX  ((double)1E32) // #/(m2 sec) maximum heterogeneous rate
#define HOM_RATE_MAX  ((double)1E41) // #/(m3 sec) maximum homogeneous rate
#define NUCRATE_MIN	((double)1E-40)


//_________________________________________________
// Public functions
//_________________________________________________
void Nucleation (void);
void NucleationCleanup (void);
void NucleationInit (void);
void NucleateHeterogeneous(CELL &cellNew, CELL &cellOld, unsigned int nodeLiquid, unsigned int nodeOther, const int dirSoltoLiq, const double tInterface);
void NucleateHomogeneous(CELL &cellNew, CELL &cellOld);
void NucleationEventReport(int i, int j, int k, double nucNumber);


//_________________________________________________
// Public variables
//_________________________________________________

#ifdef EXT_LEVEL
#define EXTERN
#else
#define EXTERN extern
#endif

#undef EXT_LEVEL
EXTERN BMATRIX CanHetNucleate;
EXTERN BMATRIX CanHomNucleate;
EXTERN DMATRIX DensityHetNuc;		//heterogeneous nucleation density
EXTERN DMATRIX DensityHomNuc;       //homogeneous nucleation density
#undef EXTERN

#endif
