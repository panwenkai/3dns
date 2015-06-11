//  _________________________________________________________
// |
// |   Parsefunc.h
// |
// |
// |   (C) 1998-99  Columbia University, MSME
// |   (C) 2000-01	Harvard University, DEAS
// |__________________________________________________________
#ifndef _PARSEFUNCH
#define _PARSEFUNCH

//_________________________________________
// typedefs
//_________________________________________
#define T_UNKNOWN		((int)0x0)
#define T_BARRIER		((int)0x0001)	//potential barrier function
#define T_BEGIN 		((int)0x0002)	//beginning of a block
#define T_BOOL			((int)0x0004)	//bool (array or single)
#define T_EDGEWORTH		((int)0x0008)	//edgeworth function (modified gaussian)
#define T_END			((int)0x0010)	//end of a block
#define T_EXPONENTIAL   ((int)0x0020)	//a + b Exp (c x)
#define T_FILE			((int)0x0040)	//general filename
#define T_FLOAT 		((int)0x0080)	//double (array or single)
#define T_GAUSSIAN		((int)0x0100)	//a=fwhm, b=mu, c=E0
#define T_INT 			((int)0x0200)	//int (array or single)
#define T_INTERPOLATION	((int)0x0400)	//interpolated number pairs
#define T_MATRIX1D		((int)0x0800)	//1D interpolation matrix
#define T_MATRIX2D		((int)0x1000)	//2D interpolation matrix
#define T_MATRIX3D		((int)0x2000)	//3D matrix of floats
#define T_POLYNOMIAL	((int)0x4000)	//a + bx + cx^2 + ...
#define T_POWER_SERIES	((int)0x8000)	//a x^b + c x^d + e x^f + ...
#define T_RPN			((int)0x10000)	//reverse polish notation
#define T_RAW			((int)0x20000)	//raw numeric data (no keyword)
#define T_STRING		((int)0x40000)	//general string (array or single)

#define T_CONST (T_INT | T_FLOAT)	//a constant scalar
#define T_FUNCTION (T_BARRIER | T_EDGEWORTH | T_EXPONENTIAL | T_GAUSSIAN | T_INTERPOLATION | T_POLYNOMIAL | T_POWER_SERIES)

#define VAL_NEGATIVE		((int) 0x0001)
#define VAL_ZERO			((int) 0x0002)
#define VAL_POSITIVE		((int) 0x0004)
#define VAL_FLAG_ALL		((int) 0x0008)
#define VAL_FLAG_INFINITY	((int) 0x0010)

#define VAL_ANY				(VAL_NEGATIVE | VAL_ZERO | VAL_POSITIVE)
#define VAL_NON_NEGATIVE	(VAL_ZERO | VAL_POSITIVE)
#define VAL_NON_POSITIVE	(VAL_NEGATIVE | VAL_ZERO)

#define MAX_MEMBERS 	 32		// maximum number of members in a general list
#define MAX_ENVIRONMENT	 32		// maximum number of elements in initialization file
#define NUM_ANY 	0      //any number of elements are allowed

#define INTFLAG_END ((int) -1)
#define INTFLAG_ALL ((int) 0xfffffffe)
#

//_________________________________________
// typedefs
//_________________________________________
typedef double STACK;

typedef struct        // RPN operator
{
   int op_code;
   double value;
} OP;

typedef struct
{
	int dataType;
	int numTerms;
	double *param;
} FUNC;
			
typedef struct			//general RPN function
{
	int nops;
	int max_ops;
	OP *ops;
	int nstack;
	int max_stack;
	STACK *stack;
} RPN;


//_________________________________________
// Typedefs
//_________________________________________
class GENERIC
{
	friend void GenTable(GENERIC &gValue, double **target, int iStart, int iStop);

	public:
		&GENERIC();							//constructor
		GENERIC(const GENERIC &rhs);	    //copy constructor
		~GENERIC();							//destructor
		
		GENERIC& operator=(const GENERIC &rhs);	//assignment operator
		
		double Evaluate(double x, double y = 0.0);

		void Load(bool *b);
		void Load(bool **bArray);
		void Load(int *i);
		void Load(int **iArray);
		void Load(double *d);
		void Load(double **dArray);
		void Load(string *str);
		void Load(string **strList);
		void Load(MATRIX2D *matx2d);
		void Load(FUNC *func);
		void Load(RPN *rpn);

		void Store(const bool b);
		void Store(const bool *bArray, int num);
		void Store(const int i);
		void Store(const int *iArray, int num);
		void Store(const double d);
		void Store(const double *dArray, int num);
		void Store(const string str);
		void Store(const string *strList);
		void Store(MATRIX2D *matx2d);
		void Store(FUNC *func);
		void Store(RPN *rpn);

		bool Valid();
		void Zero();

		int dataType;
	
	private:
		int numElements;

		union		  // GENERIC numerical value
		{
			int *iArray;    
			double *dArray;	
			string *strList;  //pointer to string list
			MATRIX2D *matx2d;	//pointer to matrix2d
			bool *bArray;
			FUNC *func;
			RPN	*rpn;
		} gValue;

}; //endclass


class ENTRY
{
		friend void ParseEntry(ENTRY &newEntry, const string *strList);
		friend void ParseLine(ENTRY &entry, int typeSpec, int numSpec, int valSpec, string *pathSearch);
		friend void ParseLine(ENTRY &entry, int typeSpec, int *numCheck, int valSpec, string *pathSearch);
		friend bool ParseVector(ENTRY &entry, int typeSpec, int numSpec, int valSpec);
		friend double EvalFunction(ENTRY &entry, double x);

	public:
		ENTRY(const ENTRY&);				//copy constructor
		ENTRY(); //constructor

		int typeFound;    	//derived type
		int numElements;	//derived number of elements
		string keyWord;		//keyword string
		GENERIC value;		//generic parsed value

		bool Valid();
		void Zero();
		
	private:
		int dataIndex;		//index to beginning of data
		string *list;   	//raw (temporary) pointer string data 
};


//_________________________________________
// Public functions
//_________________________________________
void EraseComment( char *l);
void GenTable(GENERIC &gValue, double **target, int iStart, int iStop, double factor=1.0);
int IsComment( char *s);
int IsBlank( char *s);

void ParseEntry(ENTRY &newEntry, const string *strList);
void ParseLine(ENTRY &entry, int typeSpec, int numSpec, int valSpec, const string *pathSearch);
void ParseLine(ENTRY &entry, int typeSpec, int *numCheck, int valSpec, const string *pathSearch);
bool ParseVector(ENTRY &entry, int typeSpec, int numSpec, int valSpec);
bool ParseWord(char **strPtr, string &result);

bool StringListAdd(string *strList, string &s, int listNum);
string *StringListCat(const string *str1, const string *str2);
void	StringListCopy(string **strTarget, const string *strSource);
int		StringListCount(const string *strList);
void	StringListDelete(string **strList);
int		StringListIndex(string **strList, const char *str);
string *StringListNew(int n);
void	StringListZero(string *strList);

#undef EXTERN

#endif
