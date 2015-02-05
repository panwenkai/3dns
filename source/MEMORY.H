//  _________________________________________________________
// |
// |  memory.h
// |
// |
// |   (C) 1998-99  Columbia University, MSME
// |   (C) 2000-01	Harvard University, DEAS
// |__________________________________________________________
#ifndef _MEMORYH
#define _MEMORYH

#define MAX_MATRICES 	64		    //reusable matrices

//_________________________________________
// Typedefs
//_________________________________________
typedef double ***DMATRIX;
typedef int ***IMATRIX;
typedef bool ***BMATRIX;

//_________________________________________
// PublicFunctions
//_________________________________________
void **Create2DArray (long x, long y, int size);
void ***Create3DArray (long x, long y, long z, int size);
void Free2DArray (void **a);
void Free3DArray (void ***a);
int IJKToIndex(const int i, const int j, const int k);
void IndexToIJK(const int a, int &i, int &j, int &k);
void MatrixCopy (DMATRIX mDest, DMATRIX mSource, int dShift=0);
void MatrixCopy (IMATRIX mDest, IMATRIX mSource, int iShift=0);
void MatrixFree (DMATRIX curArg);
void MatrixFree (IMATRIX curArg);
void MatrixFree (BMATRIX curArg);
void MatrixRecycle (DMATRIX curArg);
void MatrixRecycle (IMATRIX curArg);
void MatrixRecycle (BMATRIX curArg);
DMATRIX MatrixNew (DMATRIX * matxPtr);
IMATRIX MatrixNew (IMATRIX * matxPtr);
BMATRIX MatrixNew (BMATRIX * matxPtr);
int MatrixUsedD();
int MatrixUsedI();
int MatrixUsedB();
DMATRIX MatrixZero (DMATRIX matx);
IMATRIX MatrixZero (IMATRIX matx);
BMATRIX MatrixZero (BMATRIX matx);
void MemoryInit ();

//______________________________________________________
//
//	VAL
//		General retrieval of value from array using
//		single index notation
//______________________________________________________
inline double &VAL(DMATRIX d, int nodeIndex) 
	{ return (*(**d + nodeIndex)); };

inline int &VAL(IMATRIX i, int nodeIndex) 
	{ return (*(**i + nodeIndex)); };

inline bool &VAL(BMATRIX b, int nodeIndex) 
	{ return (*(**b + nodeIndex)); };


//_________________________________________
// PublicVariables
//_________________________________________
#endif