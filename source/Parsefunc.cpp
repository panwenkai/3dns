//  _________________________________________________________
// |
// |   Parsefunc.cpp    Parse and evaluate functions
// |
// |
// |   (C) 1998-99  Columbia University, MSME
// |   (C) 2000-01	Harvard University, DEAS
// |__________________________________________________________

#include "3dns.h"
#include "Files.h"

#define EXT_LEVEL
#include "parsefunc.h"
#undef EXT_LEVEL

#define DELIM " ,;\t\n{}\""

#define RPN_CONSTANT	((int) 2199)
#define	RPN_TEMP	((int) 2200)
#define	RPN_ADD		((int) 2201)
#define	RPN_SUBSTRACT	((int) 2202)
#define	RPN_MULTIPLY	((int) 2203)
#define	RPN_DIVIDE	((int) 2204)
#define	RPN_EXP		((int) 2205)
#define	RPN_LN		((int) 2206)
#define	RPN_LOG		((int) 2207)
#define	RPN_TEN_EXP	((int) 2208)
#define	RPN_SIN		((int) 2209)
#define	RPN_COS		((int) 2210)
#define	RPN_TAN		((int) 2211)
#define	RPN_POWER	((int) 2212)
#define	RPN_SQRT	((int) 2213)
#define	RPN_INV		((int) 2214)
#define	RPN_COSH	((int) 2215)
#define	RPN_SINH	((int) 2216)
#define	RPN_TANH	((int) 2217)
#define	RPN_U		((int) 2218)

#define	COUNT(x)  (sizeof(x) / sizeof(x[0]))
#define NO_STACK ((STACK *) NULL)
#define NO_OP	((OP *) NULL)

typedef struct        // list of strings
{
    char *token;
    int key;
}
SUBCOMMAND;

static SUBCOMMAND RPNOperators[] =
{
  {"T", RPN_TEMP},
  {"+", RPN_ADD},
  {"-", RPN_SUBSTRACT},
  {"*", RPN_MULTIPLY},
  {"/", RPN_DIVIDE},
  {"EXP", RPN_EXP},
  {"LN", RPN_LN},
  {"LOG", RPN_LOG},
  {"TEN_EXP", RPN_TEN_EXP},
  {"SIN", RPN_SIN},
  {"COS", RPN_COS},
  {"TAN", RPN_TAN},
  {"**", RPN_POWER},
  {"^", RPN_POWER},
  {"$$", RPN_POWER},
  {"SQRT", RPN_SQRT},
  {"INV", RPN_INV},
  {"COSH", RPN_COSH},
  {"SINH", RPN_SINH},
  {"TANH", RPN_TANH},
  {"U", RPN_U}
};

//_________________________________________
// Private Variables
//_________________________________________

//_________________________________________
// Private Functions
//_________________________________________
bool BlankLine (char *s);
double EvalRPN (RPN &rpn, double x);

bool ParseBoolean(char **strPtr, bool *rtnVal);
bool ParseFloat(string &strSource, double &result, int testVal);
bool ParseInt(char **strPtr, int *iStart, int *iStop, int testVal);
bool ParseRPN (string *strList, RPN &rpn);
bool ParseType(string &strSource, int &rtnVal);
bool ParseWord(char **strPtr, string &result);

void Pop1 (double * result, RPN &rpn, char *name);
void Pop2 (double * r1, double * r2, RPN &rpn, char *name);
void PushOp (OP * op, RPN &rpn);
void PushStack (double value, RPN &rpn);


//_________________________________________
// ENTRY::ENTRY Constructor
//	Explicit size specifications
//_________________________________________
ENTRY::ENTRY() //constructor
{
	memset(this, 0, sizeof(ENTRY));

}; //endmethod


//___________________________________________________
// ENTRY::Valid
//	Checks for validity of entry
//___________________________________________________
bool ENTRY::Valid()
{
	return ( this->list != NULL );

}; //endmethod


//___________________________________________________
//	ENTRY::Zero
//		Zeroes out entry structure	
//___________________________________________________
void ENTRY::Zero()
{
	value.Zero();
	keyWord.Zero();

	memset(this, 0, sizeof(ENTRY));		//set everything else to zero
}; //endmethod	


//______________________________________________________
//
//      EraseComment
//       Erase comments that appear at the end of a line
//______________________________________________________
void EraseComment( char *l)
{
  	if( l == NULL)
    	return;

  	while( *l != '\0')
    {
      	if ( ( *l == '#') || ((*l=='/') && (*(l+1)=='/')) )
			*l = '\0';
      	else
			l++;
    }; //endwhile
}; //endfunc


// ______________________________________________________________
// EvalRPN
//
// ______________________________________________________________
double EvalRPN (RPN &rpn, double val)
{
	int i;
	OP *op;
	double t1, t2;

	rpn.nstack = 0;

	for (i = 0, op = rpn.ops; i < rpn.nops; i++, op++)
		switch (op->op_code)
      {
      	case RPN_CONSTANT:
				PushStack (op->value, rpn);
				break;

      	case RPN_TEMP:
				PushStack (val, rpn);
				break;

      	case RPN_ADD:
				Pop2 (&t1, &t2, rpn, "ADD");
				PushStack (t1 + t2, rpn);
				break;

      	case RPN_SUBSTRACT:
				Pop2 (&t1, &t2, rpn, "SUBTRACT");
				PushStack (t2 - t1, rpn);
				break;

      	case RPN_MULTIPLY:
				Pop2 (&t1, &t2, rpn, "MULTIPLY");
				PushStack (t1 * t2, rpn);
				break;

         case RPN_DIVIDE:
				Pop2 (&t1, &t2, rpn, "DIVIDE");
				if (t1 == 0.0)
	  			ErrorMsg ("Divide by zero encountered (EvalRPN).\n");
				PushStack (t2 / t1, rpn);
				break;

      	case RPN_EXP:
				Pop1 (&t1, rpn, "EXP");
				t2 = (double) exp ((double) t1);
				PushStack (t2, rpn);
				break;

      	case RPN_LN:
				Pop1 (&t1, rpn, "LN");
				if (t1 <= 0.0)
	  				ErrorMsg ("Argument for LN operator not positive: %g (EvalRPN).\n",
		    		t1);
				t2 = (double) log ((double) t1);
				PushStack (t2, rpn);
				break;

      	case RPN_LOG:
				Pop1 (&t1, rpn, "LOG");
				if (t1 <= 0.0)
	  				ErrorMsg ("Argument for LOG operator not positive: %g (EvalRPN).\n",
		    		t1);
				t2 = (double) log10 ((double) t1);
				PushStack (t2, rpn);
				break;

      	case RPN_TEN_EXP:
				Pop1 (&t1, rpn, "TEN_EXP");
				t2 = (double) pow ((double) 10.0, (double) t1);
				PushStack (t2, rpn);
				break;

      	case RPN_POWER:
				Pop2 (&t1, &t2, rpn, "POWER");
				t2 = (double) pow ((double) t2, (double) t1);
				PushStack (t2, rpn);
				break;

      	case RPN_SIN:
				Pop1 (&t1, rpn, "SIN");
				t2 = (double) sin ((double) t1);
				PushStack (t2, rpn);
         		break;

        case RPN_COS:
				Pop1 (&t1, rpn, "COS");
				t2 = (double) cos ((double) t1);
				PushStack (t2, rpn);
				break;

      	case RPN_TAN:
				Pop1 (&t1, rpn, "TAN");
				t2 = (double) tan ((double) t1);
				PushStack (t2, rpn);
				break;

      	case RPN_SQRT:
				Pop1 (&t1, rpn, "SQRT");
				if (t1 < 0.0)
	  			ErrorMsg ("Negative argument to SQRT function: %g (EvalRPN).\n",
		    		t1);
				t2 = (double) sqrt ((double) t1);
				PushStack (t2, rpn);
				break;

         case RPN_INV:
				Pop1 (&t1, rpn, "INV");
				if (t1 == 0.0)
	  				ErrorMsg ("Zero argument to INV function (EvalRPN).\n");
				t2 = 1.0 / t1;
				PushStack (t2, rpn);
				break;

      	case RPN_COSH:
				Pop1 (&t1, rpn, "COSH");
				t2 = (double) cosh ((double) t1);
				PushStack (t2, rpn);
				break;

      	case RPN_SINH:
				Pop1 (&t1, rpn, "SINH");
				t2 = (double) sinh ((double) t1);
				PushStack (t2, rpn);
				break;

      	case RPN_TANH:
				Pop1 (&t1, rpn, "TANH");
				t2 = (double) tanh ((double) t1);
				PushStack (t2, rpn);
				break;

      	case RPN_U:
				Pop1 (&t1, rpn, "U");
				if (t1 > 1.0E-9)      	//if positive value
	  				PushStack (1.0, rpn);
				else
	  				PushStack (0.0, rpn);
				break;
      	default:
				ErrorMsg ("Internal error: unknown op-code %i (EvalRPN).\n",
		  		op -> op_code);
      } //endswitch

  		if (rpn.nstack == 1)
    		return (rpn.stack[0]);

  		if (rpn.nstack > 1)
    		ErrorMsg ("Extra Operands left on RPN stack (EvalRPN).\n");

  		ErrorMsg ("RPN stack underflow (EvalRPN).\n");

  		return (-1.0);
}; //endfunc


// ______________________________________________________________
// GENERIC::~GENERIC
//		Default destructor
// ______________________________________________________________
GENERIC::~GENERIC()
{
	this->Zero();
}; //endmethod


// ______________________________________________________________
// GENERIC::GENERIC
//		Default constructor
// ______________________________________________________________
GENERIC::GENERIC()
{
	memset(this, 0, sizeof(GENERIC));
}; //endmethod


// ______________________________________________________________
// GENERIC::Evaluate
//     Evaluate an argument regardless of whether
//     it is a polynomial, exponential, etc.
// ______________________________________________________________
double GENERIC::Evaluate(double x, double y)
{
	#define THIS_FUNC this->gValue.func
	
	int i;
	double result=0.0;
	
	if (!Valid())		//no functions to evaluate
		return(0.0);

	switch (dataType)
    {
		case T_BARRIER:			//1-dimensional barrier function
			result = ( (x >= THIS_FUNC->param[0]) && 
						(x < THIS_FUNC->param[1]) ) ? 
						THIS_FUNC->param[2] : 0.0;
			break;

		case T_EXPONENTIAL:		//1-dimensional exponential
			result = THIS_FUNC->param[0] +
					THIS_FUNC->param[1] * exp(THIS_FUNC->param[2] * x);
			break;
    
		case T_FLOAT:			//constant floating point
			this->Load(&result);
			break;

		case T_GAUSSIAN:		//1-dimensional gaussian
			{
				double sigma = THIS_FUNC->param[0] / sqrt(8 * log(2.0));
				
				result = THIS_FUNC->param[2] *
						exp( pow(x - THIS_FUNC->param[1], 2) * (-1.0) /
						(2 * pow(sigma, 2)) );
			}
			break;
						
		case T_INT:				//constant integer		
			this->Load(&i);
			result = (double)i;
			break;
	
		case T_MATRIX1D:		//1-dimensional interpolation matrix
			result = (this->gValue.matx2d)->Interpolate(x);
			break;

		case T_MATRIX2D:		//2-dimensional interpolation matrix
			result = (this->gValue.matx2d)->Interpolate(x, y);
			break;
		
		case T_POLYNOMIAL:		//1-dimensional polynomial
			for (i = THIS_FUNC->numTerms-1; i >= 0; i--)
			{
				result *= x;
				result += THIS_FUNC->param[i];
			}; //endloop

			break;
    
		case T_POWER_SERIES:	//1-dimensional power series
			result = 0.0;

			for (i = 0; i < THIS_FUNC->numTerms; i++)
			{
				if (THIS_FUNC->param[2*i+1] == 0.0)
					result += THIS_FUNC->param[2*i];
				else
					result += THIS_FUNC->param[2*i] * 
						pow (x, THIS_FUNC->param[2*i+1]);
			}; //endloop
	
			break;

		case T_RPN:				//1-dimensional general function
			result = EvalRPN(*(gValue.rpn), x);
			break;

		default:
			return(0.0);
	}; //endswitch

	return(result);
}; //endfunc


//______________________________________________________
//
//   GENERIC::operator=
//		Assignment operator
//______________________________________________________
GENERIC& GENERIC::operator=(const GENERIC &rhs)		//assignment operator
{
	Zero();		//zero this 
	
	numElements = rhs.numElements;
	dataType = rhs.dataType;

	if (numElements == 0)
		return(*this);

	if (dataType & T_CONST)	//constant scalars
	{
		switch (dataType)
		{
			case T_BOOL:
				this->Store(rhs.gValue.bArray, rhs.numElements);
				break;
			case T_INT:
				this->Store(rhs.gValue.iArray, rhs.numElements);
				break;
			case T_FLOAT:
				this->Store(rhs.gValue.dArray, rhs.numElements);
				break;
		};
	} //endif

	else if (dataType & T_FUNCTION)	 //an array
		this->Store(rhs.gValue.func);

	else if ((dataType & T_MATRIX2D) || (dataType & T_MATRIX1D))
		this->Store(rhs.gValue.matx2d);

	else if (dataType & T_STRING)
		this->Store((const string *)rhs.gValue.strList);

	else
		ErrorMsg("GENERIC::operator=: Can't copy this type of generic value");
	
	return(*this);
}; //endmethod


//______________________________________________________
//
//   GENERIC::Load
//		Loads data from generic data structure into
//		various data types.  Allocates new space if
//		necessary, but cannot discard pre-existing space.
//______________________________________________________
void GENERIC::Load(bool *b)
{
	if ((numElements == 1) && (dataType == T_BOOL))
		*b = *(gValue.bArray);
	return;
};

void GENERIC::Load(bool **bPtr)
{
	int i;

	if ((numElements >= 1) && (dataType == T_BOOL))
	{
		*bPtr = new bool[numElements];
	
		for (i = 0; i < numElements; i++)
			(*bPtr)[i] = gValue.bArray[i];
	}; //endif

	return;
}; // endmethod

void GENERIC::Load(int *i)
{
	if ((numElements == 1) && (dataType == T_INT))
		*i = *(gValue.iArray);
	return;
};

void GENERIC::Load(int **iPtr)
{
	int i;

	if ((numElements >= 1) && (dataType == T_INT))
	{
		*iPtr = new int[numElements+1];		//extra space is for -1 termination
	
		for (i = 0; i < numElements+1; i++)
			(*iPtr)[i] = gValue.iArray[i];
	}; //endif

	return;
}; // endmethod

void GENERIC::Load(double *d)
{
	if (numElements == 1)
	{
		switch (dataType)
		{
			case T_FLOAT:
				*d = *(gValue.dArray);
				break;
			case T_INT:	
				*d = (double)*(gValue.iArray);
				break;
		}; //endswitch
	}; //endif

	return;
};

void GENERIC::Load(double **dPtr)
{
	int i;

	if (numElements >= 1)
	{
		switch (dataType)
		{
			case T_FLOAT:
				*dPtr = new double[numElements];
				
				for (i = 0; i < numElements; i++)
					(*dPtr)[i] = gValue.dArray[i];				
				break;

			case T_INT:	
				*dPtr = new double[numElements];
				
				for (i = 0; i < numElements; i++)
					(*dPtr)[i] = (double)gValue.iArray[i];				
				break;
		}; //endswitch
	}; //endif

	return;
}; // endmethod

void GENERIC::Load(string *str)
{
	if ((numElements == 1) && (dataType == T_STRING))
		*str = *(gValue.strList);

	return;
}; //endmethod

void GENERIC::Load(string **strList)
{	
	if ((numElements >= 1) && (dataType == T_STRING))
		StringListCopy(strList, gValue.strList);

	return;
}; // endmethod

void GENERIC::Load(MATRIX2D *matx)
{
	if ((numElements == 1) && (dataType == T_MATRIX2D))	
		*matx = *(gValue.matx2d);

	return;
};

void GENERIC::Load(FUNC *func)
{
	int i;

	if ( this->numElements > 0 ) 
	{
		func->dataType = dataType;
		func->numTerms = gValue.func->numTerms;
		func->param = new double[func->numTerms];
		
		for (i = 0; i < func->numTerms; i++)
			func->param[i] = gValue.func->param[i];
	}; //endif

	return;
};

void GENERIC::Load(RPN *rpnA)		//This is leaky and must fix it later
{
	if ( numElements == 1 )	
		rpnA = gValue.rpn;

	return;
};


//______________________________________________________
//
//   GENERIC::store
//		Deep copy of value into generic storage class
//		clears out pre-existing data in generic, assuming
//		that it doesnt contain garbage
//______________________________________________________
void GENERIC::Store(const bool b)
{
	this->Store(&b, 1);
	return;
}; //endmethod

void GENERIC::Store(const bool *bPtr, int num)
{
	int i;

	Zero();		//zero current values

	if ((bPtr != NULL) && (num > 0) && (num < MAX_COLUMNS))
	{
		gValue.bArray = new bool[num];
		for (i = 0; i < num; i++)
			gValue.bArray[i] = bPtr[i];
	}

	dataType = T_BOOL;
	numElements = num;
	return;
}; // endmethod


void GENERIC::Store(const int i)
{
	this->Store(&i, 1);
	return;
}; //endmethod

void GENERIC::Store(const int *iPtr, int num)
{
	int i;

	Zero();		//zero current values

	if ((iPtr != NULL) && (num > 0) && (num < MAX_COLUMNS))
	{
		gValue.iArray = new int[num+1];
		for (i = 0; i < num; i++)
			gValue.iArray[i] = iPtr[i];
	
		gValue.iArray[i] = INTFLAG_END;			//extra value is termination flag
	}

	dataType = T_INT;
	numElements = num;
	return;
}; // endmethod

void GENERIC::Store(const double d)
{
	this->Store(&d, 1);
	return;
}; //endmethod

void GENERIC::Store(const double *dPtr, int num)
{
	int i;

	Zero();		//zero current values

	if ((dPtr != NULL) && (num > 0) && (num < MAX_COLUMNS))
	{
		gValue.dArray = new double[num];
		for (i = 0; i < num; i++)
			gValue.dArray[i] = dPtr[i];
	}

	dataType = T_FLOAT;
	numElements = num;
	return;
}; // endmethod

void GENERIC::Store(const string str)
{
	Zero();		//zero current values
		
	*(gValue.strList) = str;
	dataType = T_STRING;
	numElements = 1;
	return;
}; //endmethod

void GENERIC::Store(const string *newList)
{
	Zero();		//zero current values
	
	int num = StringListCount(newList);

	if (num <= 0)
		return;
	
	StringListCopy(&gValue.strList, newList);

	dataType = T_STRING;
	numElements = num;
	return;
}; // endmethod

void GENERIC::Store(MATRIX2D *matx)
{
	unsigned int j;

	Zero();		//zero current values
	
	dataType = T_UNKNOWN;		//default is unknown (null)

	if (!(matx->Valid()) )
		return;
	
	if ((matx->Rows() < 2) || (matx->Columns() < 2) )
		return;		//too small so exit
		
	for (j = 2; j < (unsigned int)matx->Columns(); j++)		//check for monotonic map in column 0
	{
		if ((*(matx))[0][j] <= (*(matx))[0][j-1])
				return;
	};

	if (matx->Rows() == 2)		//1D (2 x numJ) interpolation matrix	
	{
		gValue.matx2d = new MATRIX2D;	//allocate memory
		*(gValue.matx2d) = *(matx);	//deep copy into new data

		dataType = T_MATRIX1D;
		numElements = 1;
	}

	else
	{
		for (j = 2; j < (unsigned int)matx->Rows(); j++)		//check for monotonic map
			if ((*(matx))[j][0] <= (*(matx))[j-1][0])
				return;

		gValue.matx2d = new MATRIX2D;	//allocate memory
		*(gValue.matx2d) = *(matx);	//deep copy into new data

		dataType = T_MATRIX2D;
		numElements = 1;
	}; //endif

	return;
};


void GENERIC::Store(FUNC *func)		//deep copy
{
	int i;

	Zero();		//zero current values

	gValue.func = new FUNC;

	if ((func->param != NULL) && (func->numTerms > 0))
	{
		gValue.func->param = new double[func->numTerms];
		
		for (i = 0; i < func->numTerms; i++)
			gValue.func->param[i] = func->param[i];
	}

	dataType = gValue.func->dataType = func->dataType;
	numElements = gValue.func->numTerms = func->numTerms;

	return;
}; // endmethod

void GENERIC::Store(RPN *rpn)		//This is leaky and must fix it later
{
	Zero();		//zero current values
	
	gValue.rpn = rpn;
	dataType = T_RPN;

	numElements = 1;			//single function

	return;
};


//______________________________________________________
//
//   GENERIC::Valid
//		Checks for existing data
//______________________________________________________
bool GENERIC::Valid()
{
	if (this == NULL)
		return (false);

	return ((this->dataType != T_UNKNOWN) && (this->numElements > 0));
}; //endmethod


//______________________________________________________
//
//   GENERIC::Zero
//		Zeroes existing data
//______________________________________________________
void GENERIC::Zero()
{
	if ( this->Valid() )  //something here 
		switch (dataType)
		{
			case T_BOOL:
				delete [] gValue.bArray;
				break;
			case T_INT:
				delete [] gValue.iArray;
				break;
			case T_FLOAT:
				delete [] gValue.dArray;
				break;
			case T_STRING:
				StringListDelete(&gValue.strList);
				break;
			case T_MATRIX2D:
				if (gValue.matx2d->Valid())
				{
					gValue.matx2d->Reset(0,0);	//clear data
					delete gValue.matx2d;		//remove object
				};
				break;
			default:	//everything else just dump, (its leaky right now!)
				break;
		}; //endswitch

	if (this != NULL)
		memset(this, 0, sizeof(this));		//brute force zero

	return;
};
	

// ______________________________________________________________
// GenTable
//     Generate an optimized (integer index) lookup table
// ______________________________________________________________
void GenTable(GENERIC &gValue, double **target, int iStart, int iStop, double factor)
{	
	int i;

	*target = new double[iStop];

	for (i = iStart; i < iStop; i++)
		((*target)[i]) = gValue.Evaluate((double) i) * factor;

	return;
}; //endfunc


//______________________________________________________
//
// IsBlank
//______________________________________________________
int IsBlank( char *s)
{
  	int i;
  	i = strspn( s, DELIM);	       //find position of first non-whitespace char
  	return( (char) *(s + i) == '\0');
}; //endfunc


//______________________________________________________
//
// IsComment
//______________________________________________________
int IsComment( char *s)
{
  	char *t=0;

  	t = s + strspn( s, DELIM);

  	if(( *t == '#') || ((*t=='/') && (*(t+1)=='/')))
    	return( true);
  	else
    	return( false);

}; //endfunc


//______________________________________________________
//
// ParseBoolean
//
//______________________________________________________
bool ParseBoolean(string &strSource, bool &rtnVal)
{
    if( strSource == "TRUE" )
    	rtnVal = true;
  	else if( strSource == "FALSE")
    	rtnVal = false;
    else
    	ErrorMsg("Can't parse BOOLEAN on line %d", Environment.lineNumber);
	
  	return(true);
}; //endfunc


//______________________________________________________
//
//	ParseEntry
//		Copies a string list into the entry class
//		with type checking etc.  The string list must
//		be null-terminated
//______________________________________________________
void ParseEntry(ENTRY &newEntry, const string *strList)
{ 
	int rtnType;	
	double rtnVal;	
	
	newEntry.Zero();		//zero all values in entry
	
	if (!(strList)) return;	//nothing to do

	newEntry.list = (string *)strList;	//set pointer only, no new memory		
	
	newEntry.numElements = StringListCount(newEntry.list);

	if (ParseFloat(newEntry.list[newEntry.dataIndex], rtnVal, VAL_ANY))
	{
		newEntry.typeFound = T_RAW;		//just raw data
	}

	else if (ParseType(newEntry.list[newEntry.dataIndex], rtnType) )  //typed but no keyword
	{
		newEntry.typeFound = rtnType;
		newEntry.dataIndex++;
		newEntry.numElements--;
		return;
	}

	else 
	{
		newEntry.keyWord = newEntry.list[newEntry.dataIndex];		//must be keyword
		newEntry.dataIndex++;
		newEntry.numElements--;

		if (!(bool)(newEntry.list[newEntry.dataIndex]))		//no data
			newEntry.typeFound = T_UNKNOWN;	

		else if (ParseFloat(newEntry.list[newEntry.dataIndex], rtnVal, VAL_ANY))
			newEntry.typeFound = T_FLOAT;

		else if (ParseType(newEntry.list[newEntry.dataIndex], newEntry.typeFound))
		{
			newEntry.dataIndex++;
			newEntry.numElements--;
		}
		
		else	//treat as a list of words
		{
			newEntry.typeFound = T_STRING;
		};
	}; //endif

}; //endfunc


//______________________________________________________
//
//  ParseFloat
//    	Reads the next word from the current string, then
//		parses it as a double
//______________________________________________________
bool ParseFloat(string &strSource, double &result, int testVal)
{
    bool rtnVal=true;
  	char *startPtr=0;
	char *endPtr=0;
  	double x;
    string tempStr;

	if (!(bool)strSource)			//nothing to parse
		return(false);

    if (strSource == "INFINITY")
    	x = DOUBLE_INFINITY;

    else
    {
    	x = strtod((char *)strSource, &endPtr);
        if ((x == 0.0) && ((char *)strSource == endPtr))
        	rtnVal=false;
    }; //endelse

	result = x;		//set up return value

	if ((x == DOUBLE_INFINITY) && ((testVal & VAL_FLAG_INFINITY) != 0))
		return(rtnVal);		//match to 'Infinity' keyword

	if (x < 0)
		rtnVal = rtnVal && ((VAL_NEGATIVE & testVal) != 0);

	if (x == 0)
		rtnVal = rtnVal && ((VAL_ZERO & testVal) != 0);

	if (x > 0) 
		rtnVal = rtnVal && ((VAL_POSITIVE & testVal) != 0);

    return(rtnVal);
}; //endfunc



//______________________________________________________
//
// ParseInt
//    	Reads the next word from the current string, then
//		parses it as an integer.  Will handle range specifications
//		returnint it as
//______________________________________________________
bool ParseInt(string &strSource, int &iStart, int &iStop, int testVal)
{
    bool rtnVal=true;
	char *startPtr=0;
  	char *endPtr=0;
    
	string tempStr;
	
	if (!(bool)strSource)
		return(false);

	if (strSource == "ALL")
    	iStart = iStop = INTFLAG_ALL;

	else
	{
		iStart = strtol(strSource, &endPtr, 10); 		//load base-10 integer

		if ( *(endPtr) == '\0')		// only a single integer
     		iStop = iStart;

		else if ((*(endPtr) == '.') && (*(endPtr+1) == '.'))	  // special ".." range
    		iStop = strtol(endPtr+2, &endPtr, 10);

		else
			return(false);
	}; //endif

	if ((iStart == INTFLAG_ALL) && ((testVal & VAL_FLAG_ALL) != 0))
		return(rtnVal);		//match to 'All' keyword

	if (iStart < 0 || iStop < 0) 
		rtnVal = rtnVal && ((VAL_NEGATIVE & testVal) != 0);

	if ((iStart * iStop) <= 0)
		rtnVal = rtnVal && ((VAL_ZERO & testVal) != 0);

	if (iStart > 0 || iStop > 0)
		rtnVal = rtnVal && ((VAL_POSITIVE & testVal) != 0);

    return(rtnVal);
}; //endfunc


// ______________________________________________________________
// ParseLine
//   Parse a line of data
// ______________________________________________________________
void ParseLine(ENTRY &entry, int typeSpec, int numSpec, int valSpec, string *pathSearch)
{
   	if (entry.typeFound == T_UNKNOWN)  //should have found something
    		ErrorMsg("Unknown data type found on line %d", Environment.lineNumber);

    switch(entry.typeFound)
    {
        case T_BEGIN:
		case T_END:
			break;		//dont do anything
        
		case T_BOOL:
			if ((typeSpec & T_BOOL) == 0)
   				ErrorMsg("Boolean data type is not allowed with %s on line %d", entry.keyWord, Environment.lineNumber);

			if (!ParseVector(entry, T_BOOL, numSpec, VAL_ANY))
				ErrorMsg("Can't get boolean value on line %d", Environment.lineNumber);
			break;
        
		case T_FILE:	//user specifies 'file', true type is inferred from format
		{
			MATRIX2D *data2D = new MATRIX2D;  //allocate space for new matrix
				
			string fileName = entry.list[entry.dataIndex];		//get string directly
			
			if (!((bool)fileName))
				ErrorMsg("Can't get filename on line %d", Environment.lineNumber);
						
			if (typeSpec == T_STRING)		//string only-- load later
			{
				entry.typeFound = T_STRING;

				if (!ParseVector(entry, T_STRING, numSpec, VAL_ANY) )
					ErrorMsg("Can't parse STRING in line %d", Environment.lineNumber);
			}
		
			else	//load and identify file now
			{

				ReadMatrix2D(*data2D, fileName, pathSearch);

				if (!((*data2D).Valid()) )
					ErrorMsg("Cannot load numeric data file");
			
				if ((*data2D).Rows() > 2)
				{	
					entry.typeFound = T_MATRIX2D;

					if ((typeSpec & T_MATRIX2D) == 0)
   						ErrorMsg("2D data is is not allowed with %s on line %d", entry.keyWord, Environment.lineNumber);
				}

				else if ((*data2D).Rows() == 2)
				{
					entry.typeFound = T_MATRIX1D;
					
					if ((typeSpec & T_MATRIX1D) == 0)
   						ErrorMsg("1D data is is not allowed with %s on line %d", entry.keyWord, Environment.lineNumber);
				}

				else
					ErrorMsg("Cannot find proper columns in %s");

				entry.value.Store(data2D);

				if (entry.value.dataType == T_UNKNOWN)
					ErrorMsg("Bad 2D matrix, check mapping column 0 and row 0");
			}; //endif
		};
			break;

		case T_MATRIX3D:
		{
			double pageSpec;
			MATRIX3D data3D;
			MATRIX3D maps3D;
			string fileName;

			if (!((bool)entry.list[entry.dataIndex++]))
				ErrorMsg("Can't get filename on line %d", Environment.lineNumber);
			
			if (!ParseFloat(entry.list[entry.dataIndex], pageSpec, VAL_POSITIVE))
				ErrorMsg("Can't get page specification on line %d", Environment.lineNumber);
				
			ReadMatrix3D(data3D, fileName, pathSearch);

			if ( !(data3D.Valid()) || !(maps3D.Valid()) )
				ErrorMsg("Cannot load 3D file");
				
			//entry.value = (void *)Submatrix(data3D, pageNum)
		}; 				
			break;

		case T_BARRIER:
		{
			if ((typeSpec & T_BARRIER) == 0)
   				ErrorMsg("Barrier/Square function is not allowed with %s on line %d", entry.keyWord, Environment.lineNumber);

			FUNC *func = new FUNC;  //allocate new argument list			
			
			int i, argI;
			double *argP;

			argI = func->numTerms = 3;
			argP = func->param = new double[argI];

			for (i = 0; i < 3; i++)
      			if (!ParseFloat (entry.list[entry.dataIndex++], argP[i], VAL_ANY))
					ErrorMsg ("Unable to parse parameter %d for BARRIER function (%s) in line %d",
		  				i+1, entry.keyWord, Environment.lineNumber);
			
			if (argP[2] < argP[1])
				ErrorMsg ("Unable to resolve BARRIER function (%s)", entry.keyWord);

			func->dataType = T_BARRIER;
			entry.value.Store(func);
		}
      		break;

    	case T_EXPONENTIAL:
		{
			if ((typeSpec & T_EXPONENTIAL) == 0)
   				ErrorMsg("Exponential function is not allowed with %s on line %d", entry.keyWord, Environment.lineNumber);

			FUNC *func = new FUNC;  //allocate new argument list			
			
			int i, argI;
			double *argP;

			argI = func->numTerms = 3;
			argP = func->param = new double[argI];

			for (i = 0; i < 3; i++)
      			if (!ParseFloat (entry.list[entry.dataIndex++], argP[i], VAL_ANY))
					ErrorMsg ("Unable to parse parameter %d for EXPONENTIAL function (%s) in line %d",
		  				i+1, entry.keyWord, Environment.lineNumber);
			
			func->dataType = T_EXPONENTIAL;
			entry.value.Store(func);
		}
      		break;

        case T_FLOAT:
		{
			if ((typeSpec & T_FLOAT) == 0)
   				ErrorMsg("Floating point value is not allowed with %s on line %d", entry.keyWord, Environment.lineNumber);

			if (!ParseVector(entry, T_FLOAT, numSpec, valSpec))
				ErrorMsg("Can't parse numeric values in line %d", Environment.lineNumber);
		};
            break;

		case T_GAUSSIAN:
		{
			if ((typeSpec & T_GAUSSIAN) == 0)
   				ErrorMsg("Gaussian function is not allowed with %s on line %d", entry.keyWord, Environment.lineNumber);

			FUNC *func = new FUNC;  //allocate new argument list			
			
			int i, argI;
			double *argP;

			argI = func->numTerms = 3;
			argP = func->param = new double[argI];

			for (i = 0; i < 3; i++)
      			if (!ParseFloat (entry.list[entry.dataIndex++], argP[i], VAL_NON_NEGATIVE))
					ErrorMsg ("Unable to parse parameter %d for GAUSSIAN function (%s) in line %d",
		  				i+1, entry.keyWord, Environment.lineNumber);			

			func->dataType = T_GAUSSIAN;
			entry.value.Store(func);
		}
      		break;


        case T_INT:
		{			
			if ((typeSpec & T_INT) == 0)
   				ErrorMsg("Integer value is not allowed with %s on line %d", entry.keyWord, Environment.lineNumber);

			if (!ParseVector(entry, T_INT, numSpec, valSpec))
				ErrorMsg("Can't parse numeric values in line %d", Environment.lineNumber);
		};
            break;

		case T_INTERPOLATION:
		{	
			int i;

 			if ((typeSpec & T_INTERPOLATION) == 0)
   				ErrorMsg("Interpolation values are not allowed with %s on line %d", entry.keyWord, Environment.lineNumber);

			if (!ParseVector(entry, T_FLOAT, NUM_ANY, valSpec))
				ErrorMsg("Can't parse numeric values in line %d", Environment.lineNumber);

			if (entry.numElements < 2)
				ErrorMsg("Must have at least one coordinate pair in %s on line %d", entry.keyWord, Environment.lineNumber);

			double *temp;
			entry.value.Load(&temp);	//pull out numeric data

			MATRIX2D *data2D = new MATRIX2D(entry.numElements/(int)2, 2, 0.0);
			
			for (i = 0; i < entry.numElements; i+=2)
				data2D->Overlay( temp+i, 2, i/2, 0);

			delete [] temp;

			if (i != entry.numElements)
				ErrorMsg("Cannot have an odd number of elements in %s on line %d", entry.keyWord, Environment.lineNumber);

			entry.value.Store(&(data2D->Transpose()) );
			entry.typeFound = T_MATRIX1D;  //its become a matrix now

		};
            break;

    	case T_POLYNOMIAL:
		{
			if ((typeSpec & T_POLYNOMIAL) == 0)
   				ErrorMsg("Polynomial function is not allowed with %s on line %d", entry.keyWord, Environment.lineNumber);

			FUNC *func = new FUNC;  //allocate new argument list			
			
			int i, argI;
			double *argP;

			if (!ParseInt (entry.list[entry.dataIndex++], argI, i, VAL_NON_NEGATIVE))
				ErrorMsg ("Unable to get order of %s polynomial in line %d",
		  			entry.keyWord, Environment.lineNumber);

			func->numTerms = argI + 1;		//order + 1

			argP = func->param = new double[argI+1];

      		for (i = 0; i <= argI; i++)
			{
				if (!ParseFloat (entry.list[entry.dataIndex++], argP[i], VAL_ANY))
	 		 		ErrorMsg ("Unable to parse %s, coefficient number %i in line %d",
                		entry.keyWord, i+1, Environment.lineNumber);
			}; //endloop i
			
			func->dataType = T_POLYNOMIAL;
			entry.value.Store(func);
		};
      		break;

    	case T_POWER_SERIES:
		{
			if ((typeSpec & T_POWER_SERIES) == 0)
   				ErrorMsg("Power series function is not allowed with %s on line %d", entry.keyWord, Environment.lineNumber);

			FUNC *func = new FUNC;  //allocate new argument list			
			
			int i, argI;
			double *argP;

      		if (!ParseInt (entry.list[entry.dataIndex++], argI, i, VAL_POSITIVE))
				ErrorMsg ("Unable to get number of terms for %s in line %d",
                	entry.keyWord, Environment.lineNumber);

			func->numTerms = argI;
			argP = func->param = new double[(argI+1)*2];

     		for (i = 0; i <= argI; i++)
			{
	  			if (!ParseFloat (entry.list[entry.dataIndex++], argP[2*i], VAL_ANY))
	    			ErrorMsg ("Unable to parse %s coefficient number %d in line %d",
                    	entry.keyWord, i+1, Environment.lineNumber);
			  	if (!ParseFloat (entry.list[entry.dataIndex++], argP[2*i+1], VAL_ANY))
	    			ErrorMsg ("Unable to parse %s exponent number %d in line %d",
                    	entry.keyWord, i+1, Environment.lineNumber);
			}; //endloop i
			
			func->dataType = T_POWER_SERIES;
			entry.value.Store(func);
		};
      		break;

		case T_RAW:			//int or float raw data
			if (!ParseVector(entry, T_FLOAT, numSpec, valSpec))
				ErrorMsg("Can't parse numeric values in line %d", Environment.lineNumber);
			break;

    	case T_RPN:
		{
			if ((typeSpec & T_RPN) == 0)
   				ErrorMsg("RPN function is not allowed with %s on line %d", entry.keyWord, Environment.lineNumber);

			RPN *rpn = new RPN;  //allocate new argument list
			string *rpnList=0;

			if (!ParseVector(entry, T_STRING, NUM_ANY, VAL_ANY))
				ErrorMsg ("Unable to parse %s RPN expression in line %d",
                		entry.keyWord, Environment.lineNumber);
			
			entry.value.Load(&rpnList);  //load into stringlist
			ParseRPN(rpnList, *rpn);

			entry.value.Store(rpn);   //store back in entry
		};
      		break;

		case T_STRING:
			if ((typeSpec & T_STRING) == 0)
   				ErrorMsg("String value is not allowed with %s on line %d", entry.keyWord, Environment.lineNumber);

			if (!ParseVector(entry, T_STRING, numSpec, VAL_ANY) )
				ErrorMsg("Can't parse STRING in line %d", Environment.lineNumber);
            break;


        default:
        	ErrorMsg("Bad data on line %d", Environment.lineNumber);
    }; //endswitch

//	   	if (entry.typeFound == T_UNKNOWN)  //should have found something
//   		ErrorMsg("Unknown data type found on line %d", Environment.lineNumber);

  	return;
}; //endfunc


void ParseLine(ENTRY &entry, int typeSpec, int *numCheck, int valSpec, string *pathSearch)
{
	ParseLine(entry, typeSpec, NUM_ANY, valSpec, pathSearch); //get line
	
	if (!entry.numElements)
		ErrorMsg("No elements found in line %d", Environment.lineNumber);

	else if (!(*numCheck))						//zero value means unassigned
		*(numCheck) = entry.numElements;	//so use this value
	
	else if (*(numCheck) != entry.numElements)
		ErrorMsg("Must have %d elements in <%s> line %d", *(numCheck), 
			entry.keyWord, Environment.lineNumber);

	return;
}; //endfunc


// ______________________________________________________________
// ParseRPN
//   Parses a 'Reverse-Polish-Notation' expression
// ______________________________________________________________
bool ParseRPN (string *rpnList, RPN &rpn)
{
  	int i=0,j;
  	OP op;

    char *rpnPtr=0;

  	SUBCOMMAND *opList;		//list of operators
    double rpnConst;

  	rpn.nops = 0;
  	rpn.max_ops = 0;
  	rpn.ops = NO_OP;
  	rpn.nstack = 0;
  	rpn.max_stack = 0;
  	rpn.stack = NO_STACK;

	while (rpnList[i])
	{
    	if (ParseFloat(rpnList[i], rpnConst, VAL_ANY))		//if its a numeric value
        {
        	op.op_code = RPN_CONSTANT;
            op.value = rpnConst;
            PushOp(&op, rpn);
        }

        else
        {
            j = 0;
            opList = RPNOperators;
            
			do
        	{
				if (strcmp(rpnList[i],opList[j].token)==0)	//found a match
					break;

				if (j++ == (int)COUNT(RPNOperators) )
            		ErrorMsg("Unknown operator '%s' in line %d", *(rpnList+i), Environment.lineNumber);

        	} while (j<(int)COUNT(RPNOperators));

            op.op_code = (*(opList + j)).key;
            PushOp(&op, rpn);
        }; //endelse

		i++;
    }; //endloop i

  	return (true);
}; //endfunc



//______________________________________________________
//
//   ParseType
//		Returns the integer type specification of a string
//______________________________________________________
bool ParseType(string &strSource, int &rtnVal)
{
	rtnVal = T_UNKNOWN;

	if( (strSource == "BOOL") || (strSource == "BOOLEAN") )
    	rtnVal = T_BOOL;
  	else if( (strSource == "INT") || (strSource == "INTEGER") )
    	rtnVal = T_INT;
	else if( strSource == "INTERPOLATION")
		rtnVal = T_INTERPOLATION;
  	else if( (strSource == "FLOAT") || (strSource == "DOUBLE") )
    	rtnVal = T_FLOAT;
  	else if( strSource == "FILE")
    	rtnVal = T_FILE;
	else if( (strSource == "BARRIER") || (strSource == "SQUARE") )
		rtnVal = T_BARRIER;
  	else if( strSource == "EXPONENTIAL")
    	rtnVal = T_EXPONENTIAL;
	else if( strSource == "GAUSSIAN")
		rtnVal = T_GAUSSIAN;
  	else if( strSource == "POLYNOMIAL")
    	rtnVal = T_POLYNOMIAL;
  	else if( strSource == "POWER_SERIES")
    	rtnVal = T_POWER_SERIES;
  	else if( strSource == "RPN")
    	rtnVal = T_RPN;
  	else if( (strSource == "STRING") || (strSource == "WORD") )
    	rtnVal = T_STRING;
  	else if( strSource == "BEGIN")
    	rtnVal = T_BEGIN;
  	else if( strSource == "END")
    	rtnVal = T_END;

	return (rtnVal != T_UNKNOWN);
}; //endfunc


//______________________________________________________
//
//   ParseVector
//   Parses a vector of any type of data, specified by typespec
//	 returns the number of elements found
//______________________________________________________
bool ParseVector(ENTRY &entry, int typeSpec, int numSpec, int valSpec)
{
    int elemNum=0;
	int elemList= entry.dataIndex;
    int iRange, rangeEnd;				// special parameter for integer ranges
	
	bool *xB;		//temporary values
    int *xI;
	double *xD;
	string *xS; 

	switch (typeSpec)
	{
		case T_BOOL:
			xB = new bool[entry.numElements];				

			while (entry.list[elemList])
				elemNum += ParseBoolean(entry.list[elemList++], xB[elemNum]);
			
			if (elemNum > 0)
				entry.value.Store(xB, elemNum);
				
			delete [] xB;
			break;

		case T_FLOAT:
			xD = new double[entry.numElements];
			while (entry.list[elemList])
				elemNum += ParseFloat(entry.list[elemList++], xD[elemNum], valSpec);
			
			if (elemNum > 0)
				entry.value.Store(xD, elemNum);
				
			delete [] xD;
			break;

		case T_INT:
			xI = new int[MAX_COLUMNS];
			while (entry.list[elemList])
			{
				ParseInt(entry.list[elemList++], iRange, rangeEnd, valSpec);

				do
				{
                	xI[elemNum++] = iRange++;
                    
					if (elemNum >= MAX_COLUMNS)
					{
                    	InfoMsg("Vector is too long in line %d", Environment.lineNumber);
						return(false);
					}; //endif
				} while (iRange <= rangeEnd);
			}; //endwhile

			if (elemNum > 0)
				entry.value.Store(xI, elemNum);
			
			delete [] xI;
	  		break;

		case T_STRING:
		{
			xS = StringListNew(MAX_COLUMNS);

			StringListCopy(&xS, (string *)(entry.list + elemList));

			entry.value.Store(xS);
			elemNum = StringListCount(xS);

			StringListDelete(&xS);
			
		}
			break;
			
		default:
			ErrorMsg("Invalid data type");
	  		break;
	}; //endswitch

	if (elemNum >= MAX_COLUMNS)
	{
        InfoMsg("Vector is too long in line %d", Environment.lineNumber);
		return(false);
	}; //endif

    if (elemNum == 0)
		return(false);   //vector was empty

    if ((numSpec != NUM_ANY) && (elemNum != numSpec))
	{
    	InfoMsg ("Vector needs %d elements, but has %d, in line number %d",
        	numSpec, elemNum, Environment.lineNumber);
		return(false);
	};  //endif

	if (elemNum == 0)
	{
		InfoMsg ("Need at least 1 element in line %d", Environment.lineNumber);
		return(false);
	}; //endif

	entry.numElements = elemNum;		// number of elements read

  	return(true);
}; //endfunc


//____________________________________________________________
//
//  ParseWord
//   	Gets the next delimited word as a string of characters
//		(creates new string memory)
//____________________________________________________________
bool ParseWord(char **strPtr, string &result)
{
  	char *strNew=0;
    int lenStr;

	*(strPtr) += strspn(*(strPtr), DELIM);		// advance past whitespace (if any)

    if (**(strPtr) == '\0')    				// if no more data
    	return(false);

    lenStr = strcspn(*(strPtr), DELIM);     	// get length of word

    
	strNew = new char[lenStr + 1];          // allocate new string memory
	strncpy(strNew, _strupr(*(strPtr)), lenStr),      //all uppercase
    strNew[lenStr] = '\0';

    *(strPtr) += lenStr;                	// move ptr ahead to new position

    result = strNew;
	delete [] strNew;		//must have existing data or been properly initialized

  	return(true);
}; //endfunc


// ______________________________________________________________
// Pop1
//
// ______________________________________________________________
void Pop1 (double * result, RPN &rpn, char *name)
{
  if (rpn.nstack > 0)
    {
      rpn.nstack--;
      *result = rpn.stack[rpn.nstack];
      return;
    }
  ErrorMsg ("No argument for %s function (EvalRPN).\n", name);
}

// ______________________________________________________________
// Pop2
//
// ______________________________________________________________
void Pop2 (double * r1, double * r2, RPN &rpn, char *name)
{
  if (rpn.nstack > 1)
    {
      rpn.nstack--;
      *r1 = rpn.stack[rpn.nstack];
      rpn.nstack--;
      *r2 = rpn.stack[rpn.nstack];
      return;
    }
  ErrorMsg ("Not enough arguments for %s function (EvalRPN).\n", name);
};

// ______________________________________________________________
// PushOp
//
// ______________________________________________________________
void PushOp (OP * op, RPN &rpn)
{
	if (rpn.nops >= rpn.max_ops)
    {
      rpn.max_ops += 8;
      rpn.ops = (OP *) realloc (rpn.ops, rpn.max_ops * sizeof (OP));
    }
	rpn.ops[rpn.nops++] = *op;
};

// ______________________________________________________________
// PushStack
//
// ______________________________________________________________
void PushStack (double value, RPN &rpn)
{
	if (rpn.nstack >= rpn.max_stack)
    {
      rpn.max_stack += 8;
      rpn.stack = (double *) realloc (rpn.stack,
						 rpn.max_stack * sizeof (double));
    }
  
	rpn.stack[rpn.nstack++] = value;
}; //endfunc


// ______________________________________________________________
// StringCompare
//		Returns true if strings are identical, false otherwise
// ______________________________________________________________
bool StringCompare(const char *lhs, const char *rhs)
{
	if (lhs && rhs)
		return (!_stricmp(lhs, rhs));

	else
		return(false);
}; //endfunc


//______________________________________________________
//
// StringListAdd
//______________________________________________________
bool StringListAdd(string *strList, string &s, int listNum)
{
	int i=0;

	if (strList == 0)
		return(false);

	while (i < listNum)
	{
		if ((bool)strList[i])
			i++;

		else
		{
			strList[i] = s;
			return(true);
		}; //endif

	}; //endwhile

	return(false);	
}; //endfunc


//______________________________________________________
//
// StringListCat
//		Combines new strings into a copy string
//______________________________________________________
string *StringListCat(const string *str1, const string *str2)
{
	string *rtnStr;
	int numStr;
	int i=0;
	int j=0;

	numStr = StringListCount(str1) + StringListCount(str2);

	if (!numStr) return (0);

	rtnStr = new string[numStr+1];  //allocate and zero new string

	if (str1)
		while(*(str1+j))
		{
			*(rtnStr + i) = *(str1 + j);
			i++, j++;
		};

	j=0;

	if (str2)
		while(*(str2+j))
		{
			*(rtnStr + i) = *(str2 + j);
			i++, j++;
		};

	return(rtnStr);
}; //endfunc


//______________________________________________________
//
// StringListCopy
//		Returns the number of strings in a string list
//______________________________________________________
void StringListCopy(string **strTarget, const string *strSource)
{
	int i=0;
//	static string *tempStr;

	if (strTarget == 0)
	{
		strTarget = new (string *);
		*strTarget = StringListNew(0);
	}; //endif

	if (**strTarget)		//if target has some data
		delete [] (*strTarget);

	if (*strSource)		//if some data to copy
	{
		*strTarget = new string [StringListCount(strSource) + 1];
		while( *(strSource + i) )
		{
			*((*strTarget) + i) = *(strSource + i);
			i++;
		};
	};

	return;
}; //endfunc


//______________________________________________________
//
// StringListCount
//		Returns the number of strings in a string list
//______________________________________________________
int StringListCount(const string *strList)
{
	int i=0;
	
	if (strList == 0) return(0);

	while( *(strList + i++) );

	return(i-1);
}; //endfunc


//______________________________________________________
//
//  StringListZero
//
//______________________________________________________
void StringListZero(string *strList)
{
	int i=0;
	
	if (!strList) return;

	while ( (bool)*(strList+i) )	//while data exists
	{
		(*(strList+i)).Zero();			//zero individual strings in array
		i++;
	}; //endwhile

	return;
}; //endfunc


//______________________________________________________
//
//      StringListIndex
//
//______________________________________________________
int StringListIndex(string **strList, const char *str)
{
  	int i;
	
	int n = StringListCount(*strList);
	
	for (i=0; i<n; i++)
		if ( *(*(strList)+i) == str) 
			return(i);

  	return(-1);
}; //endfunc


//______________________________________________________
//
//      StringListGet
//
//______________________________________________________
/*char *StringListGet(char **strList, unsigned int i)
{
  	if (i < StringListCound(strList));
		return( *(strList + i));

	else
		return(0);
}; //endfunc
*/

//______________________________________________________
//
//      StringListNew
//
//______________________________________________________
string *StringListNew(int n)
{
    string *newList=0;
	
	newList = new string[n+1];
	memset(newList, 0, (n+1) * sizeof(string));

	return((string *)newList);
};


//______________________________________________________
//
//      StringListZero
//
//______________________________________________________
void StringListZero(string **strList)
{
	int i=0;
	
	while ((*strList)[i].Length())
	{
		delete [] (strList[i]);
		strList[i] = 0;

	}; //endwhile
}; //endfunc



