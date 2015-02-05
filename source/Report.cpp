//  _________________________________________________________
// |
// |   Report.cpp    Output reporting routines
// |
// |
// |  (C) 1998-99  Columbia University, MSME
// |  (C) 2000-01	Harvard University, DEAS
// |__________________________________________________________
#include "3dns.h"
#include "energy.h"
#include "phase.h"
#include "nucleation.h"
#include "thermal.h"
#include "files.h"

#define EXT_LEVEL
#include "report.h"
#undef EXT_LEVEL

//_________________________________________
// Private Functions
//_________________________________________
int PickRange(int **pickList, int rangeLo, int rangeHi, int rtnFlag, char *dimMsg);
void PrintArray(FILE *fPtr, DMATRIX dMatx, IMATRIX iMatx, int *iPick, const int j, const int k, double factor, int coord=1);
void PrintArray(FILE *fPtr, DMATRIX dMatx, IMATRIX iMatx, const int i, int *jPick, const int k, double factor=1);
void PrintArray(FILE *fPtr, DMATRIX dMatx, IMATRIX iMatx, const int i, const int j, int *kPick, double factor=1);
void PrintArray(FILE *fPtr, double *dArray, int *iPick, double factor);
void PrintArray(FILE *fPtr, int *iArray, int *iPick, int factor);
void PrintArray(FILE *fPtr, int *iArray, int numElements);
inline void PrintDouble(FILE *fPtr, double dVal);
void PrintHeader(FILE *fPtr, char *xtraString);
inline void PrintInt(FILE *fPtr, int iVal);
void PrintInterface(FILE *fPtr, int i, int j, int k, 
					double iFrac, double jFrac, double kFrac, 
					int direction, int pSol, int pLiq);
void PrintMap(FILE *fPtr, int **pickV);
void PrintMatrix(FILE *fPtr, DMATRIX dMatx, IMATRIX iMatx, int thickFlag, int **pickV, double factor=1.0);
void PrintPageHeader(FILE *fPtr, double curTime);
void ReportCleanup();
void ReportData();
bool ReportFileOpen(FILE **fPtr, int pickFlag, int **pickV, string suffixStr, string commentStr);
void ReportHistory();			//phase history matrices
void ReportInit();
void ReportInterface();			//instantaneous interface information
void ReportLaser();				//laser energy deposition matrices
void ReportNucleation();		
void ReportProbability();		//nucleation probability matrices
void ReportTemperature();		//temperature matrices
 
//_________________________________________
// Private variables
//_________________________________________
FILE *FpHeterogeneous = NULL;	//heterogeneous nucleation file
FILE *FpHistory = NULL;			//history flags file
FILE *FpHomogeneous = NULL;		//homogeneous nucleation file
FILE *FpInterface = NULL;		//interface file pointer
FILE *FpLaser = NULL;			//laser file
FILE *FpLog = NULL;				//log file
FILE *FpNucleation = NULL;		//instantaneous nucleation file
FILE *FpProperties = NULL;		//properties file
FILE *FpTemperature = NULL;		//temperature file pointer

int DEFAULT_PICK[2] = {0, INTFLAG_END};		//default pick coordinate is 0
 
inline int LENGTH(int *iPick)
{
	int i=0;

	if (iPick == NULL)
		return(0);

	while (*(iPick + i) >= 0 && (i < MAX_I))
		i++;

	return(i);
} //endfunc


// ______________________________________________________________
// LogInit
//
// ______________________________________________________________
void LogInit(char *fileName)
{
  	if (!FileOpenWrite(&FpLog, 
			ParseFileName(Report.dirOutput, Report.baseName, SUFFIX_LOG)))
		ErrorMsg("Cant open log file %s", fileName);

	Environment.fpLog = FpLog;		//set into global environment pointer

	return;
} //endfunc


//_________________________________________________
// PickRange
//
//    Adjusts and checks pick range parameter
//_________________________________________________
int PickRange(int **pickList, int rangeLo, int rangeHi, int rtnFlag, char *dimMsg)
{
	int i;

	if (*(pickList) == 0)			//nothing was provided by user
		*(pickList) = DEFAULT_PICK;
	
	else if (**(pickList) == INTFLAG_ALL)  //user gave 'All' keyword
	{
		*(pickList) = new int[rangeHi - rangeLo + 2];
		
		for (i = rangeLo; i < rangeHi; i++)
			*(*(pickList) + i - rangeLo) = i;

		*(*(pickList) + i - rangeLo) = INTFLAG_END;
	}

	else if(!InRange(*(pickList), 0, -1, rangeLo, rangeHi))
   			ErrorMsg("WATCH_%s out of range", dimMsg);

	for (i = 0; *(*(pickList) + i) != INTFLAG_END; i++);	//length of i
	
	return((i > 1) ? rtnFlag : 0);    //its thick
} //endfunc


// ______________________________________________________________
// PrintArray
//		Prints an array of values
// ______________________________________________________________
void PrintArray(FILE *fPtr, DMATRIX dMatx, IMATRIX iMatx, int *iPick, const int j, const int k, double factor, int coord)
{
	int i=0;

	while (*(iPick + i) >= 0)
	{
		if ((iMatx == 0) && (dMatx != 0))	//a double array
		{
			if (coord == 1)
				PrintDouble(fPtr, dMatx[*(iPick + i++)][j][k] * factor);
			else if (coord == 2)
				PrintDouble(fPtr, dMatx[j][*(iPick + i++)][k] * factor);
			else //coord == 3
				PrintDouble(fPtr, dMatx[j][k][*(iPick + i++)] * factor);
		}

		else if ((iMatx != 0) && (dMatx == 0))
		{
			if (coord == 1)
				PrintInt(fPtr, iMatx[*(iPick + i++)][j][k] * (int)factor);
			else if (coord == 2)
				PrintInt(fPtr, iMatx[j][*(iPick + i++)][k] * (int)factor);
			else //coord == 3
				PrintInt(fPtr, iMatx[j][k][*(iPick + i++)] * (int)factor);
		} //endif

	} //endwhile
	
	fprintf(fPtr, "\n");
} //endfunc

void PrintArray(FILE *fPtr, DMATRIX dMatx, IMATRIX iMatx, const int i, int *jPick, const int k, double factor)
{
	PrintArray(fPtr, dMatx, iMatx, jPick, i, k, factor, 2);
	return;
} //endfunc

void PrintArray(FILE *fPtr, DMATRIX dMatx, IMATRIX iMatx, const int i, const int j, int *kPick, double factor)
{
	PrintArray(fPtr, dMatx, iMatx, kPick, i, j, factor, 3);
	return;
} //endfunc

void PrintArray(FILE *fPtr, double *dArray, int *iPick, double factor)
{
	int i;
	for (i = 0; i < LENGTH(iPick); i++)
		PrintDouble(fPtr, dArray[*(iPick + i)] * factor);

	fprintf(fPtr, "\n");
	return;
} //endfunc

void PrintArray(FILE *fPtr, int *iArray, int *iPick, int factor)
{
	int i;
	for (i = 0; i < LENGTH(iPick); i++)
		PrintInt(fPtr, iArray[*(iPick + i)] * factor);

	fprintf(fPtr, "\n");
	return;
} //endfunc

void PrintArray(FILE *fPtr, int *iArray, int numElements)
{
	int i;

	for (i = 0; i < numElements; i++)
		PrintInt(fPtr, iArray[i]);

	fprintf(fPtr, "\n");
	return;
} //endfunc


// ______________________________________________________________
// PrintDouble
//		Prints a double value
// ______________________________________________________________
inline void PrintDouble(FILE *fPtr, double dVal)
{
	fprintf(fPtr, "%11.9g ", dVal);
	return;
} //endfunc


// ______________________________________________________________
// PrintHeader
//		Routine to open a file given mode and path
// ______________________________________________________________
void PrintHeader(FILE *fPtr, char *xtraString)
{
	time_t t;

	time(&t);			       //get current time

	if(fPtr != NULL)
	{
		fprintf(fPtr, "//____________________________________________________________\n");
		fprintf(fPtr, "// 3DNS - Numerical Simulation v%s\n", VERSION);
		fprintf(fPtr, "// J.P. Leonard, A.B. Limanov, J.S. Im\n");
		fprintf(fPtr, "// Program in Materials Science and Solid State Technology\n");
		fprintf(fPtr, "// Columbia University, New York, NY  10027\n");
		fprintf(fPtr, "//\n");
		fprintf(fPtr, "// COPYRIGHT(C) 1994-2001. All Rights Reserved.\n");
		fprintf(fPtr, "// Run at: %s", ctime(&t));
		fprintf(fPtr, "//\n");
		fprintf(fPtr, "// %s\n", xtraString);
		fprintf(fPtr, "//____________________________________________________________\n");
		fprintf(fPtr, "\n");
	}

  	else
    	ErrorMsg("'%s' File is not open", xtraString);

  	return;
} //endfunc


// ______________________________________________________________
// PrintInt
//		Prints an integer value
// ______________________________________________________________
inline void PrintInt(FILE *fPtr, int iVal)
{
	fprintf(fPtr, "%8d ", iVal);
	return;
} //endfunc


//_________________________________________________
// PrintInterface
//
//		Writes an interface entry to a file
//_________________________________________________
void PrintInterface(FILE *fPtr, int i, int j, int k, 
					double iFrac, double jFrac, double kFrac, 
					int direction, int nodeA, int nodeB)
{
	double tInterface;
	double velocity;
	tInterface = InterpolateT(i, j, k, iFrac, jFrac, kFrac);
	velocity = (A(Velocity) == 0.0) ? B(Velocity) : A(Velocity);

	PrintDouble(fPtr, i + iFrac);
	PrintDouble(fPtr, j + jFrac);
	PrintDouble(fPtr, k + kFrac);
	//test only
	PrintDouble(fPtr, iFrac);
	PrintDouble(fPtr, jFrac);
	PrintDouble(fPtr, kFrac);
	//
	PrintDouble(fPtr, (Geometry.xMap[i] + ((iFrac == 1.0) ? A(DelX) : 0)) * CM_TO_UM);	
	PrintDouble(fPtr, (Geometry.yMap[j] + ((jFrac == 1.0) ? A(DelY) : 0)) * CM_TO_UM);		
	PrintDouble(fPtr, (Geometry.zMap[k] + ((kFrac == 1.0) ? A(DelZ) : 0)) * CM_TO_UM);
	PrintDouble(fPtr, A(DelX) * CM_TO_UM);				//dimensions of node
	PrintDouble(fPtr, A(DelY) * CM_TO_UM);
	PrintDouble(fPtr, A(DelZ) * CM_TO_UM);
	PrintInt(fPtr, B(PhaseCodeSolid));				//solid phase
	PrintInt(fPtr, A(PhaseCodeLiquid));				//liquid phase
	PrintInt(fPtr, direction);
	PrintDouble(fPtr, tInterface);
	PrintDouble(fPtr, velocity * CM_TO_M);
	PrintDouble(fPtr, B(FractionSolid));
	//for testing only
	PrintDouble(fPtr,A(State));
	PrintDouble(fPtr,B(State));
	//end test section
	fprintf(fPtr, "\n");

	return;
} //endfunc


// ______________________________________________________________
// PrintMap
//       Write node mapping information into file
// ______________________________________________________________
void PrintMap(FILE *fPtr, int **pickV)
{
  	if(fPtr == NULL)
    	ErrorMsg("Output File is not open");

	fprintf(fPtr, "I_WATCH  INT ");
	PrintArray(fPtr, pickV[DIM_I], LENGTH(pickV[DIM_I]));      

	fprintf(fPtr, "J_WATCH  INT ");
	PrintArray(fPtr, pickV[DIM_J], LENGTH(pickV[DIM_J]));

	fprintf(fPtr, "K_WATCH  INT ");
	PrintArray(fPtr, pickV[DIM_K], LENGTH(pickV[DIM_K]));
 
	fprintf(fPtr, "DEL_X  FLOAT ");  //node widths [um]
	PrintArray(fPtr, DelX, 0, pickV[DIM_I], 0, 0, CM_TO_UM);

	fprintf(fPtr, "DEL_Y  FLOAT ");  //node heights [um]
	PrintArray(fPtr, DelY, 0, 0, pickV[DIM_J], 0, CM_TO_UM);

	fprintf(fPtr, "DEL_Z  FLOAT ");   //node depths [um]
	PrintArray(fPtr, DelZ, 0, 0, 0, pickV[DIM_K], CM_TO_UM);

	fprintf(fPtr, "MAP_X  FLOAT "); //node positions (left) [um]
	PrintArray(fPtr, Geometry.xMap, pickV[DIM_I], CM_TO_UM);

	fprintf(fPtr, "MAP_Y  FLOAT "); //node positions (top) [um]
	PrintArray(fPtr, Geometry.yMap, pickV[DIM_J], CM_TO_UM);

	fprintf(fPtr, "MAP_Z  FLOAT "); //node positions (front) [um]
	PrintArray(fPtr, Geometry.zMap, pickV[DIM_K], CM_TO_UM);

	fprintf(fPtr, "\n");
} //endfunc


// ______________________________________________________________
// PrintMatrix
//       Write a 2D or 3D matrix to an output file
// ______________________________________________________________
void PrintMatrix(FILE *fPtr, DMATRIX dMatx, IMATRIX iMatx, int thickFlag, int **pickV, double factor)
{
	int *iV, *jV, *kV;
	int **iPtr, **jPtr, **kPtr;	//pointers to pick vectors
	int *iStart, *jStart, *kStart;
	
	iV = pickV[DIM_I];
	jV = pickV[DIM_J];
	kV = pickV[DIM_K];

	switch(thickFlag)		//see RA p.837 for geometry
	{
		case (THICK_I):						//1D (i)	
		case (THICK_I | THICK_J):			//2D (i,j)
		case (THICK_I | THICK_J | THICK_K):	//3D (i,j,k)
			iPtr = &iV; iStart = iV;		//i across
			jPtr = &jV; jStart = jV;		//j down
			kPtr = &kV; kStart = kV;		//k sheet(s)
			break;

		case (THICK_J):						//1D (j)
		case (THICK_J | THICK_K):			//2D (j,k)
			iPtr = &jV; iStart = jV;		//j across
			jPtr = &kV; jStart = kV;		//k down
			kPtr = &iV; kStart = iV;		//i sheet
			break;

		case (THICK_I | THICK_K):			//2D (i,k)
			iPtr = &iV; iStart = iV;		//i across
			jPtr = &kV; jStart = kV;		//k down
			kPtr = &jV; kStart = jV;		//j sheet
			break;

		case (THICK_K):
			iPtr = &kV; iStart = kV;		//k across
			jPtr = &iV; jStart = iV;		//i down
			kPtr = &jV; kStart = jV;		//j sheet
			break;

		default:
			ErrorMsg("Reporting matrix must have more than 1 element");
	} //endswitch

	PrintPageHeader(fPtr, Sim.sClock.curTime);

	for (*kPtr = kStart; **kPtr >= 0; (*kPtr)++)
	{	
		for (*jPtr = jStart; **jPtr >= 0; (*jPtr)++)
		{
			for (*iPtr = iStart; **iPtr >= 0; (*iPtr)++)
			{
				if ((iMatx == 0) && (dMatx != 0))	//a double matrix
					PrintDouble(fPtr, dMatx[*iV][*jV][*kV] * factor);
				
				else if ((iMatx != 0) && (dMatx == 0))	//a double matrix
					PrintInt(fPtr, iMatx[*iV][*jV][*kV] * (int)factor);
			
			} // endloop iV			
			fprintf(fPtr, "\n");    //new line

		} //endloop jV	
		
		if (LENGTH(*kPtr) > 1)		//a 3D matrix
		{
			fprintf(fPtr, "//end sheet ");
			PrintInt(fPtr, (**kPtr));
			fprintf(fPtr, "\n");
		} //endif
	
	} //endloop kV

	return;
} //endfunc


// ______________________________________________________________
// PrintPageHeader
//		Print page header
// ______________________________________________________________
void PrintPageHeader(FILE *fPtr, double curTime)
{
	fprintf(fPtr, "\nPAGE FLOAT ");
	PrintDouble(fPtr, curTime);
	fprintf(fPtr, "\n");
	
	return;
} //endfunc

// ______________________________________________________________
// ReportCleanup
//              Closes all report files,
//              brutal but we only have these 4 files
// ______________________________________________________________
void ReportCleanup()
{
  	if(FpHeterogeneous)
   	{
      	fflush(FpHeterogeneous);
      	FileClose(FpHeterogeneous);
   	}

    if(FpHistory)
    {
    	fflush(FpHistory);
        FileClose(FpHistory);
    }

  	if(FpHomogeneous)
   	{
      	fflush(FpHomogeneous);
      	FileClose(FpHomogeneous);
   	}

	if(FpInterface)
   	{
      	fflush(FpInterface);
      	FileClose(FpInterface);
   	}

  	if(FpLaser)
   	{
      	fflush(FpLaser);
      	FileClose(FpLaser);
   	}

	if(FpLog)
   	{
      	fflush(FpLog);
      	FileClose(FpLog);
   	}

    if(FpNucleation)
    {
    	fflush(FpNucleation);
        FileClose(FpNucleation);
    }

  	if(FpTemperature)
   	{
      	fflush(FpTemperature);
      	FileClose(FpTemperature);
   	}

  	return;
} //endfunc


//_________________________________________________
// ReportData
//
//    Sends new information to each open report file
//______________________________________________
void ReportData(void)
{
	if (Sim.sClock == Report.sClock)
	{
/*		InfoMsg("Time=%.7gns dt=%.4gps dT=%.4g ", 
			Sim.sClock.curTime * S_TO_NS, Sim.sClock.curDTime * S_TO_PS,
			Sim.dTmax);

  		if (Sim.numberSlush > 0)
			InfoMsg("dI=%.4g ", Sim.dImax);
*/		
		InfoMsg("Time=%.7gns(%.7gps) ", Sim.sClock.curTime * S_TO_NS,
			Report.prevStep * S_TO_PS);

		ReportHistory();						//report history file
		ReportLaser();							//laser input data
		ReportProbability();  		 			//nucleation density or rates in specific nodes
		ReportTemperature();					//watch temperatures in specific nodes
		ReportNucleation();						//report number of nucleation events
		ReportInterface();						//interface data file
		InfoMsg("\n");

		Report.sClock.Increment();		//set up next reporting time
	}

	return;
} //endfunc


//_________________________________________________
// ReportFileOpen
//
//    Opens a single reporting file
//______________________________________________
bool ReportFileOpen(FILE **fPtr, int pickFlag, int **pickV, string suffixStr, string commentStr)
{
	string fileName;
	
	fileName = ParseFileName(Report.dirOutput, Report.baseName, suffixStr);

  	if (!FileOpenWrite(fPtr, fileName))
		ErrorMsg("Cant open file %s", fileName);

  	PrintHeader(*fPtr, commentStr);
	PrintMap(*fPtr, pickV);

	return(true);
} //endfunc
		

//_________________________________________________
// ReportHistory
//
//		Reports history output codes
//_________________________________________________
void ReportHistory()
{
 	if(!Report.watchHistory)
		return;
	
	if(FpHistory == NULL)
    	ErrorMsg("History file not open");

	PrintMatrix(FpHistory, 0, PhaseHistory, Report.phaseThick, Report.phaseWatch, 1);

   return;
} //endfunc


//_________________________________________________
// ReportInit
//
//    Opens and initializes all reporting files
//______________________________________________
void ReportInit()
{
	string fileName;
	
	if (Report.watchTemperature || Report.watchLaser)	
	{
		Report.tempThick = 0;
		Report.tempThick |= PickRange(Report.tempWatch + DIM_I, I_FIRST, I_LAST, THICK_I, "TEMPERATURE_I");
		Report.tempThick |= PickRange(Report.tempWatch + DIM_J, J_FIRST, J_LAST, THICK_J, "TEMPERATURE_J");
		Report.tempThick |= PickRange(Report.tempWatch + DIM_K, K_FIRST, K_LAST, THICK_K, "TEMPERATURE_K");
	} //endif

	if (Report.watchHistory || Report.watchProbability)
	{
		Report.phaseThick = 0;
		Report.phaseThick |= PickRange(Report.phaseWatch + DIM_I, I_FIRST, I_LAST, THICK_I, "PHASE_I");
		Report.phaseThick |= PickRange(Report.phaseWatch + DIM_J, J_FIRST, J_LAST, THICK_J, "PHASE_J");
		Report.phaseThick |= PickRange(Report.phaseWatch + DIM_K, K_FIRST, K_LAST, THICK_K, "PHASE_K");
	} //endif

	if (Report.watchInterface)
	{
		Report.interfaceThick = 0;		
		Report.interfaceThick |= PickRange(Report.interfaceWatch + DIM_I, I_FIRST, I_LAST, THICK_I, "INTERFACE_I");
		Report.interfaceThick |= PickRange(Report.interfaceWatch + DIM_J, J_FIRST, J_LAST, THICK_J, "INTERFACE_J");
		Report.interfaceThick |= PickRange(Report.interfaceWatch + DIM_K, K_FIRST, K_LAST, THICK_K, "INTERFACE_K");
	} //endif

	//___________ History report initialization _______________
	if (Report.watchHistory)
		ReportFileOpen(&FpHistory, Report.phaseThick, Report.phaseWatch, 
			SUFFIX_HISTORY, (string)"History output file binary flags ");

	//___________ Interface report initialization _______________
	if (Report.watchInterface)
	{
  		fileName = ParseFileName(Report.dirOutput, Report.baseName, SUFFIX_INTERFACE);

  		FileOpenWrite(&(FpInterface), fileName);
  		PrintHeader(FpInterface, "Interface output file");
  		fprintf(FpInterface, "// Slush nodes output as  i, j, k, X(um), Y, Z, dx(um), dy, dz, Phase0, Phase1, DirSlush(W=1,E=2,N=4,S=8), T(K), v(m/sec), FractionSolid\n\n");

		Report.prevSlush = Sim.numberSlush;
		Report.prevLiquid = Sim.numberLiquid;
	} //endif

	//___________ Laser report initialization _______________
	if (Report.watchLaser)
		ReportFileOpen(&FpLaser, Report.tempThick, Report.tempWatch, 
			SUFFIX_LASER, "Laser energy absorption file");

	//___________ Event report initialization _______________
	if (Report.watchNucleation)
	{
		Report.nucThisClock = Report.nucThisRept = 0;
    										//always produce instantaneous file
		fileName = ParseFileName(Report.dirOutput, Report.baseName, SUFFIX_NUCLEATION);
		FileOpenWrite(&(FpNucleation), fileName);
		PrintHeader(FpNucleation, "Nucleation event output file");
  		fprintf(FpNucleation, "// Nucleation events output as  i, j, k, X(um), Y, Z, dx(um), dy, dz, Phase0, Phase1, DirSlush(W=1,E=2,N=4,S=8), T(K), v(m/sec), FractionSolid\n\n");
	} //endif

	//___________ Probability report initialization _______________
	if (Report.watchProbability)
	{
		ReportFileOpen(&FpHeterogeneous, Report.phaseThick, Report.phaseWatch, 
			SUFFIX_HETEROGENEOUS, "Heterogeneous Nucleation Probabilities");

		ReportFileOpen(&FpHomogeneous, Report.phaseThick, Report.phaseWatch, 
			SUFFIX_HOMOGENEOUS, "Homogeneous Nucleation Probabilities");
	} //endif

	//___________ Temperature report initialization _______________
	if (Report.watchTemperature)	
		ReportFileOpen(&FpTemperature, Report.tempThick, Report.tempWatch, 
			SUFFIX_TEMPERATURE, "Temperature Output File");

	if (!InRange( Report.tInterval, 0, Sim.eraMax, 0, MAX_TIME))
    	ErrorMsg("REPORT_INTERVAL is out of range {0, %g}, or not specified", MAX_TIME);	

	Report.sClock.Reset(Sim.eraMax, Sim.eraDuration, Report.tInterval, Report.tInterval);

	return;
} //endfunc


//_________________________________________________
// ReportInterface
//
//  Reports snapshot of interface at a given time,
//  in the form
//      {x[m],y[m],z[um],T[K],v[cm/s],Direction,FractionSolid}
//_________________________________________________
void ReportInterface()
{
	int i, j, k;
	int iD, jD, kD;
	unsigned int nodeA;
	double smallFrac, bigFrac;
	smallFrac = 1e-12;
	bigFrac = 1. - 1e-12;

	double iPos, jPos, kPos;

	bool nowMelting = ((Sim.numberSlush + Sim.numberLiquid) > 0);
	bool previousMelting = ((Report.prevSlush + Report.prevLiquid) > 0);

	if ( nowMelting && !previousMelting )
    	InfoMsg("Melting ");
  	else if( !nowMelting && previousMelting )
    	InfoMsg("Resolidified ");

  	if(Sim.numberSlush)
    	InfoMsg("Slush:%d ", Sim.numberSlush);

	if(Sim.numberLiquid)
		InfoMsg("Liquid:%d ", Sim.numberLiquid);

  	Report.prevSlush = Sim.numberSlush;		//update stats
	Report.prevLiquid = Sim.numberLiquid;

	if (!nowMelting && previousMelting && Sim.modeResolidify)
		throw(MSG_END_PROGRAM);


  	if(!Report.watchInterface)
    	return;

    	PrintPageHeader(FpInterface, Sim.sClock.curTime);

   	//___________ SLUSH NODES _______________
	for(kD = 0; (k = Report.interfaceWatch[DIM_K][kD]) != INTFLAG_END; kD++)
	{
		for(jD = 0; (j = Report.interfaceWatch[DIM_J][jD]) != INTFLAG_END; jD++)	// check all j-elements
   		{
   			for(iD = 0; (i = Report.interfaceWatch[DIM_I][iD]) != INTFLAG_END; iD++)	// check all i-elements
			{
				nodeA = IJKToIndex(i, j, k);
				
				if (A(State) == SLUSH)		//standard slush node
				{	
					//calculate instantaneous interface position
					iPos = A(IPos); jPos = A(JPos); kPos = A(KPos);  //make sure we have up-to-date values
					
					//check if the the "slush node" contains an internode interface 					
					if (A(JPos) >= bigFrac)
						PrintInterface(FpInterface,i,j,k,0.5,1.,0.5,A(SlushDirection),nodeA,IJKToIndex(i,J_SOUTH,k));
					else if (A(IPos) >= bigFrac)
						PrintInterface(FpInterface,i,j,k,1.,0.5,0.5,A(SlushDirection),nodeA,IJKToIndex(I_EAST,j,k));
					else if (A(KPos) >= bigFrac)
						PrintInterface(FpInterface,i,j,k,0.5,0.5,1.,A(SlushDirection),nodeA,IJKToIndex(i,j,K_BACK));
					else //a real slush node!
					{
						PositionInterface (A(SlushDirection), A(FractionSolid), iPos, jPos, kPos);
						PrintInterface(FpInterface, i, j, k, iPos, jPos, kPos,A(SlushDirection), nodeA, nodeA);
					}
				} 
            		
				else if (A(State) == LIQUID)		//check for internode interfaces
				{
					if (State[I_WEST][j][k] == SOLID)
						PrintInterface(FpInterface, i, j, k, 0, 0.5, 0.5, EAST, nodeA, IJKToIndex(I_WEST, j, k));

					if (State[I_EAST][j][k] == SOLID)
						PrintInterface(FpInterface, i, j, k, 1, 0.5, 0.5, WEST, nodeA, IJKToIndex(I_EAST, j, k));

					if (State[i][J_NORTH][k] == SOLID)
						PrintInterface(FpInterface, i, j, k, 0.5, 0, 0.5, SOUTH, nodeA, IJKToIndex(i, J_NORTH, k));

					if (State[i][J_SOUTH][k] == SOLID)
						PrintInterface(FpInterface, i, j, k, 0.5, 1, 0.5, NORTH, nodeA, IJKToIndex(i, J_SOUTH, k));
				
					if (State[i][j][K_FRONT] == SOLID)
						PrintInterface(FpInterface, i, j, k, 0.5, 0.5, 0, BACK, nodeA, IJKToIndex(i, j, K_FRONT));
					
					if (State[i][j][K_BACK] == SOLID)
						PrintInterface(FpInterface, i, j, k, 0.5, 0.5, 1, FRONT, nodeA, IJKToIndex(i, j, K_BACK));
				} //endif
			} //endloop iD
		} //endloop jD
   	} //endloop kD

	return;
} //endfunc


//_________________________________________________
// ReportLaser
//
//    Output of laser energy delivered to individual nodes
//_________________________________________________
void ReportLaser()
{
	if (!Report.watchLaser)
		return;
	
	if (FpLaser == NULL)		//need to open output file
		ErrorMsg("Laser file is not open");

	PrintMatrix(FpLaser, QLaser, 0, Report.tempThick, Report.tempWatch);

	return;
} //endfunc


//_________________________________________________
// ReportNucleation
//
//		Reports instantaneous nucleation events
//_________________________________________________
void ReportNucleation()
{

	if(Report.watchNucleation && (Report.nucThisRept>0))
   		InfoMsg("Nucleated:%d  ", Report.nucThisRept);

	Report.nucThisRept = 0;	       //reset nucleation counter

	return;
} 


//_________________________________________________
// ReportNucleationEvent
//
//		Reports instantaneous nucleation events
//_________________________________________________
void ReportNucleationEvent(int nodeA, int dirToLiq)
{
	int i, j, k;
	double iPos, jPos, kPos;

	(Report.nucThisClock)++;	//number of nucleations in this clock
	(Report.nucThisRept)++;		//number of nucleations in this reporting session

	if ( !Report.watchNucleation )
    	return;		//dont make .nuc file if not watching nucleation
	
	PositionInterface(dirToLiq, 0, iPos, jPos, kPos);	//force position update
	IndexToIJK(nodeA, i, j, k);		//force recalculation of coordinates

  	if(Report.nucThisClock == 1)    		// first nucleation this clock
		PrintPageHeader(FpNucleation, Sim.sClock.curTime);

	PrintInterface(FpNucleation, i, j, k, iPos, jPos, kPos, 
						dirToLiq, nodeA, nodeA);		

	return;
} //endfunc


//_________________________________________________
// ReportProbability
//
//    	Reports a snapshot of nucleation densities
//		or probabilities
//_________________________________________________
void ReportProbability()
{
	if (!Report.watchProbability)
		return;
		
	if ((FpHeterogeneous == NULL) || (FpHomogeneous == NULL))
		ErrorMsg("Nucleation probability files are not open");

	PrintMatrix(FpHeterogeneous, DensityHetNuc, 0, Report.phaseThick, Report.phaseWatch);
	PrintMatrix(FpHomogeneous, DensityHomNuc, 0, Report.phaseThick, Report.phaseWatch);

	return;
} //endfunc


//_________________________________________________
// ReportProperties
//
//    Reports materials properties from lookup tables
//	  file is closed after write
//_________________________________________________
void ReportProperties(char *baseName)
{
	int t, phaseCode;
    string fileName;
	static double dummy[MAX_DEGREES];
	
	#define THIS_PHASE Phase[phaseCode]
	#define CHKNULL(q) if(q == NULL) q=dummy

    if(!Report.watchProperties)
    	return;

    fileName = ParseFileName(Report.dirOutput, baseName, SUFFIX_PROPERTIES);
    FileOpenWrite(&(FpProperties), fileName);
    PrintHeader(FpProperties, "Temperature dependent materials properties");

    for (phaseCode = 0; phaseCode < Sim.numberPhases; phaseCode++)
    {
        if(THIS_PHASE == NULL)		//ignore bad layer defs
			continue;

		CHKNULL(THIS_PHASE->heatCapacity);         //check that pointers are non-null
		CHKNULL(THIS_PHASE->thermConductivity);
		CHKNULL(THIS_PHASE->enthalpyH);
		CHKNULL(THIS_PHASE->indexN);
		CHKNULL(THIS_PHASE->indexK);
		CHKNULL(THIS_PHASE->emissivity);
		CHKNULL(THIS_PHASE->interfaceResponse);
		CHKNULL(THIS_PHASE->hetNucleation);
		CHKNULL(THIS_PHASE->homNucleation);

		fprintf(FpProperties, "//PHASE %s %s\n", THIS_PHASE->matlClassName, THIS_PHASE->phaseName);
		fprintf(FpProperties, "PAGE INT %d  //phase code\n", phaseCode);

		fprintf(FpProperties, "//T[K]\t Cp[J/cm3K]\t Kt[W/cmK]\t H[J/cm3]\t n\t k\t e\t IRF[cm/s]\t Ihet[#/cm2s]\t Ihom[#/cm3s]\n");

		for(t=MIN_DEGREES; t<MAX_DEGREES; t+=2)
		{
			PrintInt(FpProperties, t);
			
			PrintDouble(FpProperties, THIS_PHASE->heatCapacity[t]);
			PrintDouble(FpProperties, THIS_PHASE->thermConductivity[t]);
			PrintDouble(FpProperties, THIS_PHASE->enthalpyH[t]);
			PrintDouble(FpProperties, THIS_PHASE->indexN[t]);
			PrintDouble(FpProperties, THIS_PHASE->indexK[t]);
			PrintDouble(FpProperties, THIS_PHASE->emissivity[t]);
			PrintDouble(FpProperties, THIS_PHASE->interfaceResponse[t]);
			PrintDouble(FpProperties, THIS_PHASE->hetNucleation[t]);
			PrintDouble(FpProperties, THIS_PHASE->homNucleation[t]);

			fprintf(FpProperties, "\n");
		} //endloop t
	} //endloop phaseCode
    
	FileClose(FpProperties);

	return;
} //endfunc


//_________________________________________________
// ReportTemperature
//
//    Reports a snapshot of node temperatures.
//_________________________________________________
void ReportTemperature()
{
	if (!Report.watchTemperature)
		return;

	InfoMsg("T[%d][%d][%d]=%.4fK ", *(Report.tempWatch[DIM_I]), 
			*(Report.tempWatch[DIM_J]), *(Report.tempWatch[DIM_K]), 
			T[*(Report.tempWatch[DIM_I])][*(Report.tempWatch[DIM_J])][*(Report.tempWatch[DIM_K])] );

	if (FpTemperature == NULL)		//need to open output file
		ErrorMsg("Temperauture output file is not open");

	PrintMatrix(FpTemperature, T, 0, Report.tempThick, Report.tempWatch);

	return;
} //endfunc

