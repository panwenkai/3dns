//  _________________________________________________________
// |
// |   memory.cpp    Memory handling routines
// |
// |
// |   (C) 1998-99  Columbia University, MSME
// |   (C) 2000-01	Harvard University, DEAS
// |__________________________________________________________

#include "3dns.h" //also includes memory.h

//_________________________________________
// PrivateVariables
//_________________________________________
bool isAllocatedDmatrix[MAX_MATRICES];
bool isAllocatedImatrix[MAX_MATRICES];
bool isAllocatedBmatrix[MAX_MATRICES];

DMATRIX dmatrix[MAX_MATRICES];
IMATRIX imatrix[MAX_MATRICES];
BMATRIX bmatrix[MAX_MATRICES];

// ______________________________________________________________
// Create2DArray
//              Allocate array in a contiguous segment of memory
//              --cool trick from comp.lang.c faq
// ______________________________________________________________
void **Create2DArray (long r, long c, int size)
{
	long k;
	void **a;

	a = (void **) malloc (sizeof (void *) * r);
	
	if (a == NULL)
		ErrorMsg ("OUT OF MEMORY");
  
	a[0] = (void *) malloc (size * r * c);
  
	if (a[0] == NULL)
		ErrorMsg ("OUT OF MEMORY");
  
	for (k = 1; k < r; k++)
		a[k] = (void *) ((long) a[0] + k * c * size);

	return a;
}; //endfunc


// ______________________________________________________________
// Create3DArray
//              Allocate array in a contiguous segment of memory
// ______________________________________________________________
void ***Create3DArray (long x, long y, long z, int size)
{
	static int test=0;
	long i, j;
	void ***a;

	test++;
	a = (void ***) Create2DArray (x, y, sizeof (void *));
	**a = (void *) malloc (size * x * y * z);
	
	if (**a == 0)
		ErrorMsg ("OUT OF MEMORY");
    
	for (i = 0; i < x; i++)
		for (j = 0; j < y; j++)
			a[i][j] = (void *)((long) a[0][0] + (i * (y * z) + (j * z)) * size);
	
	return a;
}; //endfunc


// ______________________________________________________________
// Free2DArray
//              De-allocate memory in 2D array
// ______________________________________________________________
void Free2DArray (void **a)
{
	free (*a);
	free (a);
};


// ______________________________________________________________
// Free3DArray
//              De-allocate memory in 3D array
// ______________________________________________________________
void Free3DArray (void ***a)
{
	free (**a);		//dumps 3D contiguous info
	Free2DArray((void **)a);
}; //endfunc


//___________________________________________________
// IJKToIndex
//        Returns integer index into an array
//___________________________________________________
int IJKToIndex(const int i, const int j, const int k)
{
	return (((i * J_LAST * K_LAST) + (j * K_LAST) + k) );
}; //endfunc


//___________________________________________________
// IndexToIJK
//		Rebuilds ijk coordinates from integer index
//		expensive-- use sparingly !
//___________________________________________________
void IndexToIJK(const int a, int &i, int &j, int &k)
{
	int aTemp;

	i =  a / (J_LAST * K_LAST);		
	aTemp = a - i * (J_LAST * K_LAST);
	j = aTemp / K_LAST;
	k = aTemp - j * K_LAST;

	return;
}; //endfunc


//___________________________________________________
// MatrixCopy
//        Copy one matrix to another, with a shift
//___________________________________________________
void MatrixCopy (DMATRIX mDest, DMATRIX mSource, int dShift)
{
	if (dShift < 0)		       //shift left
    {
		memcpy (**mDest, **mSource - dShift,
			sizeof (double) * (I_LAST * J_LAST * K_LAST + dShift));

		memcpy (**mDest + (I_LAST * J_LAST * K_LAST + dShift), **mSource,
			sizeof(double) * (-dShift));
    }
  
	else				       //shift right
    {
		memcpy (**mDest + dShift, **mSource,
			sizeof (double) * (I_LAST * J_LAST * K_LAST - dShift));

		memcpy (**mDest, **mSource + (I_LAST * J_LAST * K_LAST - dShift),
			sizeof (double) * dShift);
    };

	return;
}; //endfunc


void MatrixCopy(IMATRIX mDest, IMATRIX mSource, int iShift)
{
	if (iShift < 0)		       //shift left
	{
		memcpy (**mDest, **mSource - iShift,
			sizeof (int) * (I_LAST * J_LAST * K_LAST + iShift));

		memcpy (mDest[0] + (I_LAST * J_LAST * K_LAST + iShift), mSource[0],
			sizeof (int) * iShift * -1);
    }	
  
	else				       //shift right
    {
		memcpy (**mDest + iShift, **mSource,
			sizeof (int) * (I_LAST * J_LAST * K_LAST - iShift));

		memcpy (**mDest, **mSource + (I_LAST * J_LAST * K_LAST - iShift),
			sizeof (int) * iShift);
    };

	return;
}; //endfunc


//___________________________________________________
// MatrixFree
//        De-allocate matrix, and the matrix associated with it
//___________________________________________________
void MatrixFree (DMATRIX curArg)
{
  int i;

  if (curArg != NULL)
    {
      i = 0;
      while ((curArg != dmatrix[i]) && (i < MAX_MATRICES))
	{
	  i++;
	};
      if (i < MAX_MATRICES)
	{
	  isAllocatedDmatrix[i] = false;	//no longer allocated
	  Free2DArray ((void **) dmatrix[i]);
	  dmatrix[i] = (DMATRIX) NULL;
	}
    };				       //endif

  return;
};

void MatrixFree (IMATRIX curArg)
{
  int i;

  if (curArg != NULL)
    {
      i = 0;
      while ((curArg != imatrix[i]) && (i < MAX_MATRICES))
	{
	  i++;
	};
      if (i < MAX_MATRICES)
	{
	  isAllocatedImatrix[i] = false;	//no longer allocated

	  Free2DArray ((void **) imatrix[i]);
	  imatrix[i] = (IMATRIX) NULL;
	}
    };				       //endif
  return;
};

void MatrixFree (BMATRIX curArg)
{
  int i;

  if (curArg != NULL)
    {
      i = 0;
      while ((curArg != bmatrix[i]) && (i < MAX_MATRICES))
	{
	  i++;
	};
      if (i < MAX_MATRICES)
	{
	  isAllocatedImatrix[i] = false;	//no longer allocated

	  Free3DArray ((void ***) bmatrix[i]);
	  bmatrix[i] = (BMATRIX) NULL;
	}
    };				       //endif
  return;
};

//___________________________________________________
// MatrixNew
//        Get pointer to unused matrix, allocate memory
//              if necessary
//___________________________________________________
DMATRIX MatrixNew (DMATRIX * matxPtr)
{
	int i = 0;

	while (isAllocatedDmatrix[i] && (i < MAX_MATRICES))
    {
      i++;
    };
	
	if (i == MAX_MATRICES)
		ErrorMsg ("Out of Memory in MatrixNew");	//will exit
	
	if (dmatrix[i] == (DMATRIX) NULL)
    {				       //I x J+1 matrix for some xtra memory
		dmatrix[i] = (DMATRIX) Create3DArray (I_LAST, J_LAST, K_LAST, sizeof (double));
      
		if (dmatrix[i] == (DMATRIX) NULL)
			ErrorMsg ("Out of Memory in MatrixNew");	//will exit
    };				       //endif
  
	isAllocatedDmatrix[i] = true;

	*(matxPtr) = dmatrix[i];
	return (*(matxPtr));		       //pointer to new 2D array

};				       //endfunc


IMATRIX MatrixNew (IMATRIX * matxPtr)
{
	int i = 0;

	while (isAllocatedImatrix[i] && (i < MAX_MATRICES))
    { //nn: move the the current matrix, which has yet
		//to be allocated
		i++;
    };

	if (i == MAX_MATRICES)
		ErrorMsg ("Out of Memory in MatrixNew");	//will exit

	 if (imatrix[i] == (IMATRIX) NULL)
    {				       //I x J+1 matrix for some xtra memory
		imatrix[i] = (IMATRIX) Create3DArray (I_LAST, J_LAST, K_LAST, sizeof (int));
		
		if (imatrix[i] == (IMATRIX) NULL)
		ErrorMsg ("Out of Memory in MatrixNew");	//will exit
    };				       //endif

	isAllocatedImatrix[i] = true;

	*(matxPtr) = imatrix[i];
	
	return (*(matxPtr));		       //pointer to new 2D array
}; //endfunc


BMATRIX MatrixNew (BMATRIX * matxPtr)
{
	int b = 0;

	while (isAllocatedBmatrix[b] && (b < MAX_MATRICES))
    {
		b++;
    };

	if (b == MAX_MATRICES)
		ErrorMsg ("Out of Memory in MatrixNew");	//will exit

	 if (bmatrix[b] == (BMATRIX) NULL)
    {				       //I x J+1 matrix for some xtra memory
		bmatrix[b] = (BMATRIX) Create3DArray (I_LAST, J_LAST, K_LAST, sizeof (bool));
		
		if (bmatrix[b] == (BMATRIX) NULL)
		ErrorMsg ("Out of Memory in MatrixNew");	//will exit
    };				       //endif

	isAllocatedBmatrix[b] = true;

	*(matxPtr) = bmatrix[b];
	
	return (*(matxPtr));		       //pointer to new 2D array
}; //endfunc	


//___________________________________________________
// MatrixRecycle
//        De-allocate matrix, but not memory associated with it
//___________________________________________________
void MatrixRecycle (DMATRIX curArg)
{
	int i;

	if (curArg != NULL)
    {
		i = 0;
		
		while ((curArg != dmatrix[i]) && (i < MAX_MATRICES))
		{
			i++;
		};
      
		if (i < MAX_MATRICES)
			isAllocatedDmatrix[i] = false;
    }; //endif
}; //endfunc

void MatrixRecycle (IMATRIX curArg)
{
	int i;

	if (curArg != NULL)
    {
		i = 0;
		
		while ((curArg != imatrix[i]) && (i < MAX_MATRICES))
		{
			i++;
		};
      
		if (i < MAX_MATRICES)
			isAllocatedImatrix[i] = false;

    }; //endif
}; //endfunc


void MatrixRecycle (BMATRIX curArg)
{
	int b;

	if (curArg != NULL)
    {
		b = 0;
		
		while ((curArg != bmatrix[b]) && (b < MAX_MATRICES))
		{
			b++;
		};
      
		if (b < MAX_MATRICES)
			isAllocatedBmatrix[b] = false;

    }; //endif
}; //endfunc


//___________________________________________________
// MatrixUsed
//        Calculate how many matrices are in use
//        good for diagnosing memory leaks
//___________________________________________________
int MatrixUsedD ()
{
	int i = 0, numUsed = 0;

	while (i < MAX_MATRICES)
		numUsed += isAllocatedDmatrix[i++];

	return (numUsed);
}; //endfunc

int MatrixUsedI ()
{
	int i = 0, numUsed = 0;

	while (i < MAX_MATRICES)
		numUsed += isAllocatedImatrix[i++];
	return (numUsed);
}; //endfunc

int MatrixUsedB ()
{
	int b = 0, numUsed = 0;

	while (b < MAX_MATRICES)
		numUsed += isAllocatedImatrix[b++];
	return (numUsed);
}; //endfunc


//___________________________________________________
// MatrixZero
//        Set all elements of a matrix to a specific value
//___________________________________________________
DMATRIX MatrixZero (DMATRIX matx)
{
	memset (matx[0][0], (int) 0, I_LAST * J_LAST * K_LAST * sizeof (double));
	return(matx);
}; //endfunc

IMATRIX MatrixZero (IMATRIX matx)
{
	memset (matx[0][0], (int) 0, I_LAST * J_LAST * K_LAST * sizeof (int));
	return(matx);
}; //endfunc

BMATRIX MatrixZero (BMATRIX matx)
{
	memset (matx[0][0], (int) 0, I_LAST * J_LAST * K_LAST * sizeof (bool));
	return(matx);
}; //endfunc


//___________________________________________________
// MemoryInit
//        Initialize memory related variables
//___________________________________________________
void MemoryInit ()
{
	memset(dmatrix, 0, sizeof(DMATRIX) * MAX_MATRICES);
	memset(imatrix, 0, sizeof(IMATRIX) * MAX_MATRICES);
	memset(bmatrix, 0, sizeof(BMATRIX) * MAX_MATRICES);

	return;
}; //endfunc



