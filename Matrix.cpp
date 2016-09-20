//  _________________________________________________________
// |
// |   Matrix.cpp
// |
// |	Matrix C++ class structures and development
// |
// |   (C) 1998-99  Columbia University, MSME
// |   (C) 2000-01	Harvard University, DEAS
// |__________________________________________________________

#include "compat.h"
#include "Matrix.h"		//matrix class prototypes
#include "Error.h"

const int XMAP = 0;
const int YMAP = 1;

MATRIX2D NULL_MATX2D;
MATRIX3D NULL_MATX3D;

//_________________________________________
// Constructor
//	Explicit size specifications
//_________________________________________
MATRIX2D::MATRIX2D(unsigned int r, unsigned int c, double val) //constructor
{
	memset(this, 0, sizeof(MATRIX2D));
	data = (double **)BlockCreate(r, c, val);

	nameStr.Zero();
	numI = r;
	numJ = c;

}; //endmethod


//_________________________________________
// Destruction operator
//_________________________________________
MATRIX2D::~MATRIX2D()
{
	//nameStr.Zero();  ~string is automatically called on exit

	if (data) delete *data;
	delete data;
	data = 0;		//force NULL evaluation

};


//_________________________________________
// Assignment operator
//	Deep copy with smart memory conservation
//_________________________________________
MATRIX2D &MATRIX2D::operator=(const MATRIX2D &rhs)
{
	if (this == &rhs) return *this;		//for self-assignment
	
	if (this->data)					//clear out any old data
	{
		delete this->data[0];
		delete this->data;
	};

	this->nameStr.Zero();
	memset(this, 0, sizeof(MATRIX2D));

	this->dataType = rhs.dataType;
	this->numI = rhs.numI;
	this->numJ = rhs.numJ;
	this->nameStr = rhs.nameStr;

	if (rhs.data)					//if something to copy
	{
		this->data = BlockCreate(numI, numJ);	//allocate new block
		BlockCopy(this->data[0], rhs.data[0], 
			rhs.numI * rhs.numJ * sizeof(double));		//copy data into new block
	};
		
	return *this;
		
}; //endmethod


//_________________________________________
// Copy constructor
//	Uses assignment operator
//_________________________________________
MATRIX2D::MATRIX2D(const MATRIX2D &rhs)
{
	memset(this, 0, sizeof(*this));		//zero everything
	operator=(rhs);
}; //endmethod
	

//___________________________________________________
// Index operator
//	  Return direct pointer to row, doesnt check
//___________________________________________________
double *MATRIX2D::operator[] (unsigned int r)
{
	if (this->numI) return(this->data[r]);
	else return(0);

}; //endmethod


//_________________________________________
// BlockCreate 
//	Allocates continguous block of memory
//_________________________________________
double **MATRIX2D::BlockCreate(unsigned int r, unsigned int c, double val)
{
	unsigned int k;
	void **a;
		
	if((r==0) || (c==0))
		return(NULL);

	a = (void **) new (void *[r]);
	
	if (a == NULL)
		ErrorMsg ("OUT OF MEMORY");
  
	a[0] = (void *) new double[r * c];
  
	if (a[0] == NULL)
		ErrorMsg ("OUT OF MEMORY");
  
	for (k = 1; k < r; k++)
		a[k] = (void *) ((long) a[0] + k * c * sizeof(double));

	if (val == 0.0)
		memset (a[0], (int) 0, r * c * sizeof (double));

	else if (val < DOUBLE_INFINITY)
		for (k=0; k<r*c; k++) *((double *)a[0]+k) = val;

	return((double **)a);
}; //endmethod
	

//___________________________________________________
// BlockCopy
//	Copy one memory block to another, with a shift
//___________________________________________________
void MATRIX2D::BlockCopy(double *mDest, double *mSource, long int numBytes, long int numShift)
{
	if (mSource==NULL || mDest==NULL) 
		return;
	
	if (numShift < 0)		       //shift left
    {
		memcpy (mDest, (char *)mSource - numShift,
			numBytes + numShift );

		memcpy (mDest + numBytes + numShift, mSource,
			-numBytes);
    }
  
	else				       //shift right
    {
		memcpy ((char *)mDest + numShift, mSource,
			numBytes - numShift );

		memcpy (mDest, (char *)mSource + numBytes - numShift,
			numShift);
    };
}; //endmethod


//___________________________________________________
// MATRIX2D::Columns
//	Returns number of columns
//___________________________________________________
int MATRIX2D::Columns()
{
	return (this->Valid() ? this->numJ : 0);

}; //endmethod


//___________________________________________________
// MATRIX2D::CopyRow
//	Copies a row of data into a double array
//___________________________________________________
double *MATRIX2D::CopyRow(unsigned int rowNum)
{
	double *retArray;

	if (!Valid() || (rowNum >= numI) ) 
		return(0);

	retArray = new double[numJ];

	memcpy(retArray, this->data+rowNum, numJ*sizeof(double));

	return(retArray);
};

	

/*
//___________________________________________________
// MATRIX2D::LookUp
//	Look up a value in a special 2-series matrix
//	Assumes that x series is constant or 
//  monotonically increasing
//___________________________________________________
double MATRIX2D::Interpolate(double x)
{
	unsigned int iHi=0;
	unsigned int iLo=0;

	if (!(this->Valid()) )
		return(0.0);

	if ((this->numI != 2) || (this->numJ < 1))
		return(0.0);

	while ( (iHi < numJ) && (this->data[0][iHi]) < x )
		iHi++;

	if (iHi == numJ)		//beyond end
		return (this->data[1][numJ-1]);	//return last value

	else if (iHi == 0)
		return (this->data[1][0]);

	else // (iHi > 0)
	{
		iLo = iHi-1;
		
		if (data[0][iLo] == data[0][iHi])	//prevents divide-by-zero
			return (data[0][iLo]);

		else
			return (data[1][iLo] + (x - data[0][iLo]) * 
					(data[1][iHi] - data[1][iLo]) / 
					(data[0][iHi] - data[0][iLo]) );
	}; //endif

	return(0.0);
}; //endfunc
 
*/

//______________________________________________________
//
// Interpolate
//	 Returns an interpolated value from a 3D surface
//	 and checks the validity of matrices
//
//______________________________________________________
double MATRIX2D::Interpolate(double x, double y)
{
	unsigned int iHi, jHi;		//node coordinates for interpolation
	unsigned int iNew, iOld, jNew, jOld;
	unsigned int iFirst, jFirst;
		
	double iFrac, jFrac;		//fractional coordinates;

	#define yMap(i) data[i + iFirst][0]
	#define xMap(j) data[0][j + jFirst]		//transposition implied
	#define	dataMap(i, j) data[i + iFirst][j+jFirst]

	if (!(this->Valid()) )
		return(0.0);

	if ((this->numI < 2) || (this->numJ < 1))  //not a typical matrix
		return(0.0);

	jFirst = ((numI == 2) ? 0 : 1);		//first j coordinate is different for 1D and 2D	
	jNew = this->jLo;

	do
	{
		jOld = jNew;

		if (jNew > 0)
			if (xMap(jNew) > x)
				jNew--;

		if (jNew < (this->numJ) - jFirst - 1)
			if (xMap(jNew+1) <= x)
				jNew++;

	} while (jOld != jNew);

	jLo = jNew;		//save this position for next call
	jHi = jLo + ((jLo < (this->numJ) - jFirst - 1) ? 1 : 0);

	jFrac = ((xMap(jLo) >= x) || (xMap(jLo) == xMap(jHi))) ? 0.0 : 
			(xMap(jHi) < x) ? 1.0 : (x - xMap(jLo)) / (xMap(jHi) - xMap(jLo));

	if (this->numI > 2)
	{
		iFirst = 1;			//begin data in second i-series
		iNew = this->iLo;

		do
		{
			iOld = iNew;

			if (iNew > 0)
				if (yMap(iNew) > y)
					iNew--;

			if (iNew < (this->numI) - iFirst - 1)
				if (yMap(iNew+1) <= y)
					iNew++;

		} while (iOld != iNew);

		iLo = iNew;		//save this position for next call		
		iHi = iLo + ((iLo < (this->numI) - iFirst - 1) ? 1 : 0);

		iFrac = ((yMap(iLo) >= y) || (yMap(iLo) == yMap(iHi))) ? 0.0 : 
				(yMap(iHi) < y) ? 1.0 : (y - yMap(iLo)) / (yMap(iHi) - yMap(iLo));

		return( dataMap(iLo, jLo) * (1.0 - iFrac) * (1.0 - jFrac) +
			dataMap(iLo, jHi) * (1.0 - iFrac) * jFrac +
			dataMap(iHi, jLo) * iFrac * (1.0 - jFrac) +
			dataMap(iHi, jHi) * iFrac * jFrac );

	}
	
	else	//1D only

		return( data[1][jLo+jFirst] * (1.0 - jFrac) +
			data[1][jHi+jFirst] * jFrac );

}; //endfunc


//___________________________________________________
// MATRIX2D::Name
//	Public access to string identifier
//___________________________________________________
string& MATRIX2D::Name()
{
	return(this->nameStr);
}; //endmethod

//___________________________________________________
// MATRIX2D::Overlay
//	Copy one memory block to another, with a shift
//___________________________________________________
bool MATRIX2D::Overlay(MATRIX2D &matx, int rowNum, int colNum)
{
	int r, c;

	if ((!this->Valid()) || (!matx.Valid()) )
		return(false);

	if ((rowNum >= this->Rows()) || (colNum >= this->Columns()) )
		return(false);

	for (r=rowNum; ((r < this->Rows()) && (r-rowNum < matx.Rows()) ); r++)
		for (c=colNum; ((c < this->Columns()) && (c-colNum < matx.Columns()) ); c++)
			
			this->data[r][c] = matx.data[r-rowNum][c-colNum];

	return(true);
};

bool MATRIX2D::Overlay(double *array, int numElements, int rowNum, int colNum)
{
	int c;

	if ((!this->Valid()) || (array == NULL) )
		return(false);

	if ((rowNum >= this->Rows()) || (colNum >= this->Columns()) )
		return(false);

	for (c=colNum; ((c < this->Columns()) && (c-colNum < numElements) ); c++)
			
			this->data[rowNum][c] = array[c-colNum];

	return(true);
};


//___________________________________________________
// MATRIX2D::Pick
//	Non-contiguous submatrix of the parent in both directions
//___________________________________________________
MATRIX2D& MATRIX2D::Pick(int *iPick, unsigned int iNum, int *jPick, unsigned int jNum)
{
	if ((iNum > 0) && (jNum > 0))		//pick both directions
		return((this->PickI(iPick, iNum)).PickJ(jPick, jNum));
	
	else if (iNum > 0)		//pick rows only
		return(this->PickI(iPick, iNum));

	else if (jNum > 0)		//pick columns only
		return(this->PickJ(jPick, jNum));

	else					//copy whole matrix
		return(*this);
}; //endfunc


//___________________________________________________
// MATRIX2D::PickI
//	Non-contiguous submatrix of the parent
//___________________________________________________
MATRIX2D& MATRIX2D::PickI(int *iPick, unsigned int iNum)
{
	static MATRIX2D newMatx;
	
	unsigned int i, iNew;

	if (!this->Valid() )	//bad matrix
		return(NULL_MATX2D);

	if ((iNum <= 0) || (iPick == NULL))		//invalid iPick
		return(NULL_MATX2D);

	newMatx.Reset(iNum, this->Columns());	//clear space for new data

	for (i = 0, iNew = 0; i < iNum; i++)
	{
		if ((*(iPick+i) < 0) || (*(iPick+i) >= (int)this->numI) )
			continue;		//out of range

		newMatx.Overlay(this->data[*(iPick+i)], this->numJ, iNew, 0);
		iNew++;
	}; //endloop i

	return(newMatx.SubMatrix(0, 0, iNew, this->Columns() ));
}; //endmethod


//___________________________________________________
// MATRIX2D::PickJ
//	Non-contiguous submatrix of the parent
//___________________________________________________
MATRIX2D& MATRIX2D::PickJ(int *jPick, unsigned int jNum)
{
	static MATRIX2D newMatx;
	
	unsigned int i, j, jNew;

	if (!this->Valid() )	//bad matrix
		return(NULL_MATX2D);

	if ((jNum <= 0) || (jPick == NULL))		//invalid iPick
		return(NULL_MATX2D);

	newMatx.Reset(this->Rows(), jNum);	//clear space for new data

	for (j = 0, jNew = 0; j < jNum; j++)
	{
		if ((*(jPick+j) < 0) || (*(jPick+j) >= (int)this->numJ) )
			continue;		//out of range

		for (i = 0; i < (unsigned int)this->Rows(); i++)
			newMatx.data[i][jNew] = this->data[i][*(jPick+j)];

		jNew++;
	}; //endloop j

	return(newMatx.SubMatrix(0, 0, this->Rows(), jNew ));
}; //endmethod


//___________________________________________________
// MATRIX2D::QueryI
//	Selects i-rows based on query condition 
//  in column jQuery
//___________________________________________________
MATRIX2D& MATRIX2D::QueryI(unsigned int jQuery, double *dBound, unsigned int dNum)
{
	unsigned int i, iNum, d;
	static int *iPick = 0;

	if (iPick != NULL)			//delete data from previous
	{
		delete [] iPick;
		iPick = 0;
	}

	if (!this->Valid() )		//bad matrix
		return(NULL_MATX2D);

	if ((jQuery < 0) || (jQuery >= this->numJ))		//invalid query index
		return(NULL_MATX2D);

	iPick = new int[this->numI];				//space for pick array

	for (i=0, iNum=0; i < this->numI; i++)	//build pick array
	{
		for (d=0; d < dNum-1; d+=2)
		{
			if ((this->data[i][jQuery] >= dBound[d]) && 
				(this->data[i][jQuery] < dBound[d+1]) )
			{
				iPick[iNum++] = i;
				break;
			}; //endif
		}; //endloop d
	}; //endloop i
	
	return(this->PickI(iPick, iNum));	//call explicit PickI
}; //endmethod


//___________________________________________________
// MATRIX2D::QueryJ
//	Selects j-columns based on query condition 
//  in column iQuery
//___________________________________________________
MATRIX2D& MATRIX2D::QueryJ(unsigned int iQuery, double *dBound, unsigned int dNum)
{
	unsigned int j, jNum, d;
	static int *jPick = 0;

	if (jPick != NULL)			//delete data from previous
	{
		delete [] jPick;
		jPick = 0;
	}

	if (!this->Valid() )		//bad matrix
		return(NULL_MATX2D);

	if ((iQuery < 0) || (iQuery >= this->numI))		//invalid query index
		return(NULL_MATX2D);

	jPick = new int[this->numJ];				//space for pick array

	for (j=0, jNum=0; j < this->numJ; j++)	//build pick array
	{
		for (d=0; d < dNum-1; d+=2)
		{
			if ((this->data[iQuery][j] >= dBound[d]) && 
				(this->data[iQuery][j] < dBound[d+1]) )
			{
				jPick[jNum++] = j;
				break;
			}; //endif
		}; //endloop d
	}; //endloop i
	
	return(this->PickJ(jPick, jNum));	//call explicit PickI
}; //endmethod


//___________________________________________________
// MATRIX2D::Reset
//	Re-initialize pointer to new size
//___________________________________________________
void MATRIX2D::Reset(unsigned int r, unsigned int c, double val)
{
	MATRIX2D newMatx(r, c, val);
	*this = newMatx;
}; //endmethod


//___________________________________________________
// MATRIX2D::Rows
//	Returns number of rows
//___________________________________________________
int MATRIX2D::Rows()
{
	return (this->Valid() ? this->numI : 0);

}; //endmethod


//___________________________________________________
// MATRIX2D::Size
//	Returns number of elements in the matrix
//___________________________________________________
int MATRIX2D::Size()
{
	if (! Valid() )
		return(0);
	
	return(numI * numJ);
}; //endmethod


//___________________________________________________
// MATRIX2D::SubMatrix
//	Returns a submatrix of the parent at {{r0,c0},{r1,c1}}
//___________________________________________________
MATRIX2D &MATRIX2D::SubMatrix(unsigned int r0, unsigned int c0, 
							  unsigned int r1, unsigned int c1)
{
	static MATRIX2D newMatx;

	unsigned int r, c;
	double **tempData = this->data;

	if (!this->Valid() )	//bad matrix
		return(NULL_MATX2D);

	r1 = (r1 > this->numI) ? (this->numI) : r1;	//adjust ranges
	c1 = (c1 > this->numJ) ? (this->numJ) : c1;
	
	if ((r0 >= r1) || (c0 >= c1))	//nothing to copy
			return(NULL_MATX2D);
		
	newMatx.Reset(r1 - r0, c1 - c0);

	for (r=r0; r < r1; r++)
		for (c=c0; c < c1; c++)
			newMatx.data[r-r0][c-c0] = this->data[r][c];

	newMatx.nameStr = this->nameStr;	//copy name

	return(newMatx);   //reference to static so okay
};


//___________________________________________________
// MATRIX2D::Transpose
//	Copy one memory block to another, with a shift
//___________________________________________________
MATRIX2D &MATRIX2D::Transpose()
{
	unsigned int r, c;
	double **tempData = this->data;

	static MATRIX2D newMatx;
	static MATRIX2D nullMatx;

	if (!this->Valid())			//bad source matrix
		return(nullMatx);

	this->data = 0;	// hide source data for shallow copy
	newMatx = *(this);	
	this->data = tempData;  //restore original

	newMatx.numJ = this->numI;
	newMatx.numI = this->numJ;
	newMatx.data = BlockCreate(newMatx.numI, newMatx.numJ);

	for (r=0; r < newMatx.numI; r++)
		for (c=0; c < newMatx.numJ; c++)
			newMatx.data[r][c] = this->data[c][r];

	return(newMatx);   //reference to static so okay
};


//___________________________________________________
// MATRIX2D::Valid
//	Check that matrix is valid
//___________________________________________________
bool MATRIX2D::Valid()
{
	if (this == NULL)
		return(false);

	return ((this->data) &&
		(this->numI > 0) && (this->numJ > 0) );
};


//_________________________________________
// MATRIX3D Constructor
//	Explicit size specifications
//_________________________________________
MATRIX3D::MATRIX3D(unsigned int p) //constructor
{
	memset(this, 0, sizeof(MATRIX3D));			//zero values

	numK = p;

	if (p)
		page = new MATRIX2D [p];
}; //endmethod


//_________________________________________
// MATRIX3D Destructor
//_________________________________________
MATRIX3D::~MATRIX3D()
{
	int i;

	if (!Valid())
		return;
	
	for (i = 0; i < MAX_MAPS; i++)
		mapArray[i].Reset(0,0);		//zero out maps
	
	if (page) delete [] (this->page);	//zero out pages

}; //endmethod


//_________________________________________
// MATRIX3D Assignment operator
//	Deep copy with smart memory conservation
//_________________________________________
MATRIX3D &MATRIX3D::operator=(const MATRIX3D &rhs)
{
	unsigned int p;

	if (this == &rhs) return *this;		//for self-assignment

	delete [] (this->page);		//remove old 2D information
	this->page = 0;

	this->maxI = rhs.maxI;		//assumes rhs is good (or zero) data
	this->maxJ = rhs.maxJ;
	this->numK = rhs.numK;

	if (rhs.page)					//if something to copy
	{
		this->page = new MATRIX2D [rhs.numK];

		for (p=0; p<rhs.numK; p++)
			*((this->page)+p) = *((rhs.page)+p);	//copy 2D matrices
	};
	
	for (p=0; p < MAX_MAPS; p++)					//copy mapping arrays
		*((this->mapArray)+p) = *((rhs.mapArray)+p); 

	return *this;
		
}; //endmethod


//_________________________________________
// MATRIX3D Copy constructor
//	Uses assignment operator
//_________________________________________
MATRIX3D::MATRIX3D(const MATRIX3D &rhs)
{
	delete (this);				//destroy everything
	operator=(rhs);
}; //endmethod


//___________________________________________________
// MATRIX3D::operator[]		Index operator
//	Smart return of a specific page
//___________________________________________________
MATRIX2D& MATRIX3D::operator[] (unsigned int p)
{
	static MATRIX2D nullMatx;

	if (!this->page)
		return(nullMatx);
	
	if (p >= this->numK) 
		return(nullMatx);

	return(*((this->page)+p));
}; //endmethod


//___________________________________________________
// MATRIX3D::GetMap
//   Returns a mapping array from the 3d matrix
//___________________________________________________
MATRIX2D& MATRIX3D::GetMap(string &mapName)
{
	int i;

	if (this->Valid())	
		for (i = 0; i < MAX_MAPS; i++)
			if (this->mapArray[i].Name() == mapName)
				return(mapArray[i]);	//found it

	return(NULL_MATX2D);	//didnt find anything
}; //endfunc


//___________________________________________________
// MATRIX3D::Overlay
//	Copies a page into an existing 3D matrix
//___________________________________________________
void MATRIX3D::Overlay(MATRIX2D &newPage, unsigned int p)
{
	if (!Valid() || !newPage.Valid() )
		return;

	if (p >= this->numK)
		return;

	*((this->page)+p) = newPage;

}; //endmethod


//___________________________________________________
// MATRIX3D::Page
//	 Public reference to a MATRIX2D page
//___________________________________________________
MATRIX2D& MATRIX3D::Page(unsigned int pageNum)
{
	if (!Valid())
		return NULL_MATX2D;

	return (this->page[pageNum]);
}; //endmethod


//___________________________________________________
// MATRIX3D::Reset
//	Re-initialize pointer to new size
//___________________________________________________
void MATRIX3D::Reset(unsigned int p)
{
	MATRIX3D newMatx(p);
	*this = newMatx;	  //dumps old data, copies new into this
}; //endmethod


//___________________________________________________
// Sheets
//	Returns number of pages in the matrix
//___________________________________________________
int MATRIX3D::Sheets()
{
	return(numK);
}; //endmethod


//___________________________________________________
// MATRIX3D::SubMatrix
//	Returns a submatrix of the parent from pages p0 thru p1}
//___________________________________________________
MATRIX3D &MATRIX3D::SubMatrix(unsigned int p0, unsigned int p1)
{
	unsigned int p;
	static MATRIX3D newMatx;
	static MATRIX3D nullMatx;

	int newMaxI = 0;
	int newMaxJ = 0;

	if(!this->Valid() )
		return(nullMatx);

	p1 = (p1 >= this->numK) ? (this->numK)-1 : p1;	//adjust ranges
	
	if (p0 > p1)				//nothing to copy
		return(nullMatx);
		
	newMatx.Reset(p1 - p0);  //allocate pages
	
	for(p=p0; p < p1; p++)
	{
		*((newMatx.page)+p) = *((this->page)+p);		//now full copy
		newMaxI = max(newMaxI, newMatx.page[p].Rows());
		newMaxJ = max(newMaxJ, newMatx.page[p].Columns());
	};

	newMatx.numK = p1 - p0;

	for (p=0; p < MAX_MAPS; p++)		//copy map arrays
		newMatx.mapArray[p] = this->mapArray[p];

	return(newMatx);   //reference to static so okay
};


//___________________________________________________
// Valid
//	Check that matrix is valid
//___________________________________________________
bool MATRIX3D::Valid()
{
	return ((this->page) && (this->numK > 0));
};


//______________________________________________________
//
// Ceiling
//	 Ceiling function for 2D matrix
//______________________________________________________
bool Ceiling(MATRIX2D &data, double yThreshold, double yCeiling)
{
	int i, j;
		
	if ( !data.Valid() ) 
		return (false);
		
	for (i = 0; i < data.Rows(); i++)
	{
		for (j = 0; j < data.Columns(); j++)
		{
			if (data[i][j] >= yThreshold)
				data[i][j] = yCeiling;
		}; //endloop j
	}; //endloop i

	return(true);
}; //endfunc	

//______________________________________________________
//
// Floor
//	 Floor function for 2D matrix
//______________________________________________________
bool Floor(MATRIX2D &data, double yThreshold, double yFloor)
{
	int i, j;
		
	if ( !data.Valid() ) 
		return (false);
		
	for (i = 0; i < data.Rows(); i++)
	{
		for (j = 0; j < data.Columns(); j++)
		{
			if (data[i][j] < yThreshold)
				data[i][j] = yFloor;
		}; //endloop j
	}; //endloop i

	return(true);
}; //endfunc	


//______________________________________________________
//
// Normalize
//	 Normalizes y-data from range y0 to y1
//______________________________________________________
bool Normalize(MATRIX2D &data, double y0, double y1)
{
	int i, j;
	double yMin=DOUBLE_INFINITY;
	double yMax=-DOUBLE_INFINITY;
	double yScale;
		
	if ( !data.Valid() ) 
		return (false);

	for (i = 0; i < data.Rows(); i++)		//find max and minimum values
	{
		for (j = 0; j < data.Columns(); j++)
		{
			if (data[i][j] < yMin)		
				yMin = data[i][j];

			if (data[i][j] > yMax)
				yMax = data[i][j];
		}; //endloop j
	}; //endloop i

	if (yMax != yMin)
	{
		yScale = (y1 - y0) / (yMax - yMin);
		
		for (i = 0; i < data.Rows(); i++)		//do normalization
			for (j = 0; j < data.Columns(); j++)
				data[i][j] = y0 + (data[i][j] - yMin) * yScale;
	}

	else if ((yMax == yMin) && (yMax > 0.0))	//constant positive -> goto y1
	{
		for (i = 0; i < data.Rows(); i++)		//do normalization
			for (j = 0; j < data.Columns(); j++)
				data[i][j] = y1;
	}

	else //all zero values so undefined
		return(false);

	return(true);
}; //endfunc	


//______________________________________________________
//
// Resample
//	 Resamples data on the new domain xNew
//	 Assumes that 'data' has the form	ROW0: {x0,x1,..xn}
//										ROW1: {y0,y1,..yn}
//______________________________________________________
bool Resample(MATRIX2D &target, MATRIX2D &data, MATRIX2D &xNew)
{
	int i, iW, iE, iMax;
	double *yNew;
	double *x, *y;

	static MATRIX2D newMatx;	
	
	if ( !data.Valid() || !xNew.Valid() ) 
		return (false);
	
	if ( (data.Rows() != 2) || (xNew.Rows() != 1))
		return (false);

	yNew = new double[xNew.Columns()];
	
	x = data[XMAP];		//for speed
	y = data[YMAP];
	iMax = xNew.Columns();

	for(i=iW=iE=0; i<iMax; i++)
	{
		while((x[iE] < xNew[XMAP][i]) && (iE < iMax-1))
			iW = iE++;
		
		if (x[iE] != x[iW])
			yNew[i] = y[iW] + (xNew[XMAP][i] - x[iW]) * (y[iE] - y[iW]) / (x[iE] - x[iW]);

		else
			yNew[i] = y[iE];	//out of bounds or at singularity
	
	}; //endloop i
	
	target.Reset(2, iMax);   //allocate new data
	
	target.Overlay(xNew[XMAP], iMax, XMAP, 0);
	target.Overlay(yNew, iMax, YMAP, 0);

	return(true);
}; // endfunc


// ______________________________________________________________
// string::string
//		Default constructor
// ______________________________________________________________
string::string(unsigned int n)
{
	memset(this, 0, sizeof(string));

	if (n)
		this->data = new char [n];
}; //endmethod


// ______________________________________________________________
// string::string
//		Constructor for typecast from char *
// ______________________________________________________________
string::string(const char *str)
{
	memset(this, 0, sizeof(string));

	if (str)
	{
		this->data = new char [strlen(str)+1];
		strcpy(this->data, str);
	};
}; //endmethod


// ______________________________________________________________
// string::string
//		Constructor for typecast from integer
// ______________________________________________________________
string::string(const int i)
{
	char buffer[80];

	this->Zero();
	
	*this = _itoa(i, buffer, 10);
}; //endmethod


// ______________________________________________________________
// string::string
//		Default destructor
//			called implicitly when exiting scope of a string variable
//			or can be called explicitly using 'delete strPtr'
//			can also be used to fully remove an array of strings,
//			use 'delete[] strPtr'
// ______________________________________________________________
string::~string()
{
	this->Zero();
}; //endmethod


//_________________________________________
// string::operator bool
//	 boolean typecast from string
//_________________________________________
string::operator bool()
{
	if (this)
		if ( this->data )
			return (true);
	
	return (false);
}; //endmethod

string::operator bool() const
{
	if (this)
		if ( this->data )
			return (true);
	
	return (false);
}; //endmethod


//_________________________________________
// string::operator char *
//	Assignment operator with deep copy and new
//_________________________________________
string::operator char*()
{
	if (!this)
		return(0);

	if ((this->data) && *(this->data))
		return (this->data);

	else
		return (0);
}; //endmethod

string::operator const char*() const
{
	if (!this)
		return(0);

	if ((this->data) && *(this->data))
		return (this->data);

	else
		return (0);
}; //endmethod



//_________________________________________
// string::operator=
//	Assignment operator with deep copy and new
//_________________________________________
string& string::operator=(const string &rhs)	//assignment operator
{
	if (this == &rhs) return *this;		//for self-assignment
	
	this->Zero();					//clear out any old data

	if (rhs.data)
	{
		this->data = new char [strlen(rhs.data)+1];
		strcpy(this->data, rhs.data);
	}; //endif

	return *this;
}; //endmethod


string* string::operator=(const string *rhs)	//assignment operator for array
{
	if (this == rhs) return this;
	
	return((string *)rhs);

}; //endmethod


//_________________________________________
// string::operator+
//	Concatenation of strings
//_________________________________________
string string::operator+(const string &rhs)	//assignment operator
{
	char temp[MAX_CHARS*2];
	int i1, i2;
	
	i1 = this->Length();
	i2 = rhs.Length();
	
	if ((!(bool)rhs) || (!(bool)*this))
		return (*this);

	*temp = (char) 0;

	strcat(temp, this->data);
	strcat(temp, rhs.data);

	return ((string) temp);
}; //endfunc


//_________________________________________
// string::Length()
//		Returns number of characters in string
//_________________________________________
int string::Length() const		// <<const after means that 'this' is constant
{
	if (!this)
		return(0);

	if (data)
		return (strlen(data));

	else
		return(0);
}; //endmethod
	

//_________________________________________
// string::string	
//	Copy constructor -- uses assignment operator
//_________________________________________
string::string(const string &rhs)
{
	memset(this, 0, sizeof(*this));		//zero everything
	operator=(rhs);
}; //endmethod


//_________________________________________
// string::operator=
//	test for equality of strings
//_________________________________________
bool string::operator==(const string &rhs)
{
	if (this == &rhs)
		return(true);

	else if (!this)
		return(false);

	else if (data && rhs.data)
		return (!_stricmp(data, rhs.data));

	else 
		return(false);

}; //endmethod

bool string::operator==(const char *rhs)
{
	return(*this == (string)rhs);
};


//_________________________________________
// string::operator new[]
//	 Handles clean initialization of string array
//_________________________________________
void* string::operator new[] (size_t numElements)
{
	size_t memSz = (numElements + 1) * sizeof(string);  //array of string elements
	void *pTemp = malloc(memSz);
	
	if (pTemp != 0)						//set all parts of strings to zero
		memset(pTemp, 0, memSz);

	return(pTemp);
}; //endmethod


void StringListDelete(string **strList)
{
	delete [] *strList;		//zeroes & deallocates each member, then deallocates array
	*strList = 0;
	return;
};


//_________________________________________
// string::Replace
//		Replaces a substring with a new substring
//_________________________________________
char* string::Replace(const string &oldStr, const string &newStr)
{
	int i;
	int lenOld, lenNew;

	char *oldPtr;
	char *newPtr;

	static char rtnChars[MAX_CHARS];
	memset(rtnChars, 0, MAX_CHARS);		//forget old string

	if ((!(bool)*this) || !(bool)oldStr)
		return (*this);

	oldPtr = this->data;		//old string
	newPtr = rtnChars;

	lenOld = oldStr.Length();
	lenNew = newStr.Length();
	
	while ( (unsigned int)(i = strcspn(oldPtr, oldStr.data)) != strlen(oldPtr) )	//distance to next replacement
	{
		strncpy(newPtr, oldPtr, i);		//copy up to replacement point
		newPtr += i;
		
		strncpy(newPtr, newStr, lenNew);
		newPtr += lenNew;
		oldPtr += i + lenOld;
	};

	strcat(newPtr, oldPtr);		//finish string

	return(rtnChars);
}; //endmethod


//_________________________________________
// string::Zero
//	 Zeroes out a string member
//_________________________________________
void string::Zero()
{
	if (this->data)		//assumes that a nonzero address indicates valid data
		delete [] ((char *)this->data);	    //delete array of characters

	memset(this, 0, sizeof(string));
}; //endmethod
		








	 
	
	

