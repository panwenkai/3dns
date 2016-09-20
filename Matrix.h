//  _________________________________________________________
// |
// |   Matrix.h    Matrix Class structures
// |
// |   (C) 1999  Columbia University, MSME
// |   (C) 2000-01	Harvard University, DEAS
// |__________________________________________________________
#ifndef _MATRIXCH
#define _MATRIXCH

#include <malloc.h>
#include <memory.h>		//memory allocation etc

#define MAX_CHARS	((int) 0x100)	 // maximum number of chars in a string
#define MAX_COLUMNS ((int) 0x2000) //maximum columns in 2D matrix
#define MAX_MAPS	((int) 0x10)	//maximum number of mapping arrays
#define MAX_PAGES ((int) 0x1000)	 //maximum pages in 3D matrix
#define MAX_ROWS ((int) 0x8000)    //maximum rows in 2D matrix
#define MAX_2D	((int) 0x80000)	 //maximum elements in 2D matrix
#define MAX_3D  ((int) 0x200000)	 //maximum elements in 3D matrix	
#define MAX_2D_MATRICES	((int) 0x100)  //maximum 2D matrices allocatable
#define MAX_3D_MATRICES ((int) 0x20)	 //maximum 3D matrices allocatable

#define M_DOMAIN ((int) 0)
#define M_RANGE ((int) 1)

#define T_RAW1D		((int) 0x0001)	//a raw 1d matrix of floats  (see RA p.896)
#define T_RAW2D		((int) 0x0002)	//a raw 2d matrix of floats
#define T_MAPPED1D	((int) 0x0004)	//1d data with map index
#define T_MAPPED2D	((int) 0x0008)	//2d data with dual map index

class string
{
	public:
		string(unsigned int n=0);		//constructor
		string(const char *str);        //constructor for typcast from char *
		string(const int i);			//constructor for typecast from int
		string(const string &rhs);	    //copy constructor
		~string();				        //destructor
		
		
		string& operator=(const string &rhs);	//assignment operator
		string* operator=(const string *rhs);	//assignment operator for array
		string operator+(const string &rhs);	//concatenating strings
		bool operator==(const string &rhs);		//equality operator for 'string'
		bool operator==(const char *rhs);		

		string& operator new;
		void* operator new[] (size_t numElements); //allocate clean array
		
		operator bool();		//typecast from 'string' to 'bool'
		operator bool() const;  //typecast from 'const string' to 'bool'
		operator char*();		//typecast from 'string' to 'char *'
		operator const char*() const; //typecast from 'const string' to 'char *'
		
		int Length() const;		//return length of string
		void Zero();
		char* Replace(const string &oldStr, const string &newStr);

	private:
		char *data;
}; //endclass


class MATRIX2D
{
	public:
		MATRIX2D(const MATRIX2D&);				//copy constructor
		MATRIX2D(unsigned int i=0, unsigned int j=0, double val=DOUBLE_INFINITY); //constructor
		~MATRIX2D();							//destructor
		MATRIX2D& operator=(const MATRIX2D&);	//assignment operator

		void Reset(unsigned int i, unsigned int j, double val=DOUBLE_INFINITY);
		
		bool Overlay(MATRIX2D &matx, int rowNum, int colNum);
		bool Overlay(double *array, int numElements, int rowNum, int colNum);
		bool Valid(void);

		int	Columns(void);
		int	Rows(void);
		int Size(void);

		double* operator[](unsigned int r);
		double* CopyRow(unsigned int rowNum);
		double Interpolate(double x, double y=0.0);	//2-dimensional interpolation

		string& Name();

		MATRIX2D& Pick(int *iPick, unsigned int iNum, int *jPick, unsigned int jNum);	//pick i AND/OR j
		MATRIX2D& PickI(int *iPick, unsigned int iNum); //matrix rows
		MATRIX2D& PickJ(int *jPick, unsigned int jNum); //matrix columns
		MATRIX2D& QueryI(unsigned int jQuery, double *dbound, unsigned int dNum);
		MATRIX2D& QueryJ(unsigned int iQuery, double *dbound, unsigned int dNum);

		MATRIX2D& SubMatrix(unsigned int r0, unsigned int c0, unsigned int r1, unsigned int c1);
		MATRIX2D& Transpose(void);	//create a transposed matrix
		
	private:
		string nameStr;
		int dataType;		

		unsigned int numI;			//number of i-coordinates
		unsigned int numJ;			//number of j-coordinates
		unsigned int iLo;			//last referenced i-coordinate
		unsigned int jLo;			//last referenced j-coordinate

		double **data;

		double **BlockCreate(unsigned int r, unsigned int c, double val=DOUBLE_INFINITY);
		void BlockCopy(double *mDest, double *mSource, long int numBytes, long int numShift=0);

}; //endclass


class MATRIX3D
{
	friend bool ReadMatrix3D(MATRIX3D &newMatx, const char *fileName, string *pathSearch);

	public:
		MATRIX3D(unsigned int p=0); //constructor
		~MATRIX3D();							//destructor
		MATRIX3D(const MATRIX3D&);				//copy constructor
		MATRIX3D& operator=(const MATRIX3D&);	//assignment operator
		MATRIX2D& operator[](unsigned int p);  //index operator
		
		MATRIX2D& GetMap(string& mapName);
		MATRIX2D& Page(unsigned int p);
		
		int MaxRows(void);
		int MaxCols(void);
		void Overlay(MATRIX2D &newPage, unsigned int pageNum);
		void Reset(unsigned int p);
		int Sheets(void);		//returns total number of pages in matrix
		MATRIX3D& SubMatrix(unsigned int p0, unsigned int p1);
		bool Valid();		//checks that matrix is valid


	private:
		unsigned int maxI;			//largest i-value		
		unsigned int maxJ;			//largest j-value
		unsigned int numK;			//number of pages (k values)
		unsigned int size;

		MATRIX2D mapArray[MAX_MAPS];
		MATRIX2D *page;
}; //endclass


//_________________________________________________
// Public functions
//_________________________________________________
bool Ceiling(MATRIX2D &data, double yThreshold, double yCeiling);
bool Floor(MATRIX2D &data, double yThreshold, double yFloor);
bool Normalize(MATRIX2D &data, double y0, double y1);
bool Resample(MATRIX2D &target, MATRIX2D &data, MATRIX2D &xNew);


//_________________________________________________
// Public variables
//_________________________________________________

#ifdef EXT_LEVEL
#define EXTERN
#else
#define EXTERN extern
#endif
EXTERN MATRIX2D NULL2D;		// laser energy deposited in each node
EXTERN MATRIX3D NULL3D;
#undef EXTERN

#endif