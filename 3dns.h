//  _________________________________________________________
// |
// |   3dns.h
// |
// |   (C) 1998-99  Columbia University, MSME
// |   (C) 2000-01	Harvard University, DEAS
// |_________________________________________________________
#ifndef THREEDNSH
#define THREEDNSH
#ifndef _TWODPCH
#define _TWODPCH

#include <ctype.h>
#include <sys/types.h>
#include <malloc.h>
#include <math.h>
#include <limits.h>

#include <time.h>

#endif

#include "compat.h"
#include "memory.h"
#include "Matrix.h"
#include "Parsefunc.h"
#include "Clock.h"

#define VERSION  "3.7.1"
#define BUILD	"20131104"

//_____ node states __________
#define LIQUID ((int)0x001)
#define SOLID ((int)0x002)
#define SLUSH ((int)0x004)

#define LAYER_FIRST 0
#define LAYER_LAST 	Geometry.jZones

#define I_FIRST  	Geometry.iBound[0]
#define J_FIRST  	Geometry.jBound[0]
#define K_FIRST		Geometry.kBound[0]

#define I_LAST	  	Geometry.iBound[Geometry.iZones]
#define J_LAST      Geometry.jBound[Geometry.jZones]
#define K_LAST		Geometry.kBound[Geometry.kZones]

#define I_PERIODIC false
#define J_PERIODIC false
#define K_PERIODIC false

#define LOOP_I	for(i=I_FIRST; i<I_LAST; i++)
#define LOOP_J	for(j=J_FIRST; j<J_LAST; j++)
#define LOOP_K	for(k=K_FIRST; k<K_LAST; k++)

#define MIN_DEGREES  ((int)77)		       //minimum temperature
#define MAX_DEGREES	 ((int)4000)	       //maximum degrees in lookup table

#define MAX_DIM		 ((int)(3))			 //maximum number of dimensions
#define MAX_I        ((int)(0x400))      //maximum number of i-nodes
#define MAX_J        ((int)(0x200))      //maximum number of j-nodes
#define MAX_K		 ((int)(0x100))		 //maximum number of k-nodes
#define MAX_NODES	 ((int)(0x1000000))	 //maximum total number of nodes
#define MAX_JREGIONS 16		       //maximum vertical regions (layers)
#define MAX_IREGIONS 16		       //maximum horizontal regions
#define MAX_KREGIONS 16			   //maximum front-back regions
#define MAX_REGIONS	 256		   //maximum total regions (including overlays)
#define MIN_NODESIZE ((double)5E-7)	//minimum node size [cm]
#define MAX_NODESIZE ((double)1.0001)    //maximum node size [cm]

#define MAX_PHASES		1024		//maximum number of unique phases in simulation
#define MAX_MATERIAL_CLASS	16	    //maximum number of unique material classes

#define MAX_TCHANGE	 ((double)1000.0)	//maximum temperature change
#define MAX_ICHANGE	 ((double)0.98)		//maximum fraction solid change
#define MAX_RETRIES	 ((int) 6)			//maximum number of retries for clock reduction

#define MIN_LAMBDA	 ((double)1E-7) //minimum laser wavelength [m]
#define MAX_LAMBDA	 ((double)1E-5)	//maximum laser wavelength [m]
#define MIN_ENERGY	 ((double)0.0)  //minimum laser energy density [J/cm2]
#define MAX_ENERGY	 ((double)100.0)	//maximum laser energy density [J/cm2]
#define MAX_BZ		 ((int) 3)		//maximum number of boundary zones

#define MSG_END_PROGRAM ((int) 0x01)	//signal to end program

#define WEST ((int)0x001)
#define EAST ((int)0x002)
#define NORTH ((int)0x004)
#define SOUTH ((int)0x008)
#define FRONT ((int)0x010)
#define BACK ((int)0x020)
#define HETEROGENEOUS ((int)0x040)	//special hemispherical interface
#define HOMOGENEOUS ((int)0x080)	//special spherical interface
#define SURFACEMELT	((int)0x100)
#define ALL_DIRECTIONS ((int)0x1FF)   // NESWUDHHS = 9 filled bits
#define ALL_DIRORDINARY	((int)0x03F)	// NESWUD = 6 filled bits

#define CM_TO_UM     ((double)(10000.0))
#define M_TO_UM		 ((double)(1E-6))
#define CM_TO_M		 ((double)(0.01))
#define M_TO_CM		 ((double)(100.0))

#define DIM_I ((int) 0x000)	//packed bit positions
#define DIM_J ((int) 0x001)
#define DIM_K ((int) 0x002)

#define THICK_I ((int) 0x001)
#define THICK_J ((int) 0x002)
#define THICK_K ((int) 0x004)

#define THIS_ERA Sim.sClock.eCurrent

#define CHANGE_OK ((int) 0x000)
#define CHANGE_HI ((int) 0x001)
#define CTRL_PROP ((double) 0.9)    //proportional control factor

#define LOOK_UP(table, T) (*(table+(int)T) + (T-(int)T) * (*(table+(int)T+1)-*(table+(int)T))) //lookup function for velocity T = temperature?
#define CHANGED(a,b,c)  (a.b != c.b)
#define A(x) VAL(x, nodeA) //VAL retreives matrix value using single-index notation
#define B(x) VAL(x, nodeB) //nodeA and nodeB should be an integer index of x, where x is an (I/D/B)MATRIX

#define NEW(x) VAL(x, nodeNew)
#define OLD(x) VAL(x, nodeOld)

//_________________________________________
// Typedefs
//_________________________________________
typedef double *TABLE;

typedef struct CELL
{
	IMATRIX cHistory;
	IMATRIX cState;
	IMATRIX cDirection;
	IMATRIX cPhaseLiquid;
	IMATRIX cPhaseSolid;

	DMATRIX cFraction;
	DMATRIX cTimeX;
} CELL;

typedef struct
{
	bool *canChangeI;		// if column can change
	bool *canChangeJ;		// if row can change
	bool *canChangeK;		// if page can change
	bool *modePeriodic;		//periodic boundary conditions for ijk
	bool *isThick;			// if dimension is thick

    int dimThick;		// packed integer info on dimensionality
	int iZones;			// number of i-zones
    int jZones;			// number of j-zones (layers)
	int kZones;			// number of k-zones
	int numberRegions;	// total number of regions including overlays
    int *iBound;		// i-coordinates of zones
    int *jBound;		// j-coordinates of zones
	int *kBound;		// k-coordinates of zones
	int *nodesI;		// number of horizontal columns in each i-zone
	int *nodesJ;		// number of vertical rows in each j-zone
	int *nodesK;		// number of pages in each k-zone
	
    double *delX;		// width of nodes in each i-zone
    double *delY;		// height of nodes in each j-zone
	double *delZ;		// depth of noded in each k-zone
	double *sizeI;		// total width of each i-zone
	double *sizeJ;		// total height of each j-zone
	double *sizeK;		// total depth of each k-zone

	double *xMap;		// WEST edge of each i-column
	double *yMap;		// NORTH edge of each j-row
	double *zMap;		// UP edge of each k-sheet

} GEOMETRY;

typedef struct
{
	int bZone[ALL_DIRECTIONS+1];

	int numNeighbors[MAX_DIM];			//number of neighbors in each bZone
	
	int dirKeep[MAX_BZ][ALL_DIRECTIONS+1];		//directions to clear in higher bZones
	int dirInverse[ALL_DIRECTIONS+1];		//inverse direction

	int iCrystal[ALL_DIRECTIONS+1];		//pseudo-crystallographic components
	int jCrystal[ALL_DIRECTIONS+1];
	int kCrystal[ALL_DIRECTIONS+1];

	int crystalType[ALL_DIRECTIONS+1];	//crystal class
	
} INTERFACE;	

typedef struct			// Laser parameters
{
    int **beginDielectric;		// first j pos where absorption begins
    int **beginAbsorption;	// first j pos where vacuum end and material begins
    
	double energyDensity;	// peak laser energy density [J/cm2]
	double peakTemporal;	// peak temporal flux [J/(cm2 sec)]
	double normTemporal;	// normalization factor for temporal curve
    double waveLength;		// wavelength (Meters)

	MATRIX2D egySpatial;		// indexed in i-space
	
	GENERIC dataSpatialX;		//input spatial x-profile
	GENERIC dataSpatialZ;		//input spatial z-profile
    GENERIC dataTemporal;	//input temporal coordinates
}
LASER;

typedef struct
{
	string phaseName;		//"LIQUID" or any other solid phase name
	int phaseType;			// 0=liquid, 1=primary solid, 2..n=alternate solids

    string matlClassName;	//material name used to identify class of materials
	int matlClass;			//integer associated with a unique class
    
    double *enthalpyH;   	//lookup table heat of fusion
    double *heatCapacity;	//lookup table [MAX_DEGREES]
    double *hetNucleation;	//heterogeneous nucleation
    double *homNucleation;	//homogeneous nucleation
    double *indexN;
    double *indexK;
    double *interfaceResponse;	//lookup table
    double *optAbsorption;		//lookup table
    double *thermConductivity;	//lookup table
    double *emissivity;

	int numTransitions;		//number of phases accessible from this phase
	int *transitionTemp;	//single (melting&nucleation) temp for solids, interface temp for liquid
	int *transitionTo;		//solids: liquid phase number, liquid: solid phase number(s)
} PHASE;

typedef struct
{
	int phaseLiquid;	//index to liquid phase
	int phaseStart;		//index to starting phase
	int numberI;			//number of i locations
	int numberJ;			//number of j locations
	int numberK;			//number of k locations
	int *iLocations;	//rows that this region spans
	int *jLocations;	//columns that this region spans
	int *kLocations;	//sheets that this region spans
	int *phaseIndex;	//list of phase numbers associated with this region

	bool canChange;		//can form liquid phase
    bool catalyzeFreezing;	//true if it can catalyze freezing in neighboring regions
    bool catalyzeMelting;	//true if it can catalyze melting in neighboring regions
    bool canHetNucleate;
    bool canHomNucleate;

	double thresholdHet;	//heterogeneous nucleation threshold
    double thresholdHom;	//homogeneous nucleation threshold
} REGION;

typedef struct			// Reporting information
{
	bool watchHistory;		//true if history report
    bool watchInterface;	//true if interface report
	bool watchLaser;		//true if laser report
    bool watchNucleation;	//true if report on instantaneous nucleation events
	bool watchProbability;  //true if report on nucleation probabilities
	bool watchProperties;	//true if dump lookup tables to a file
    bool watchTemperature;	//true if temperature output
	
	int nucThisClock;		//nucleations in this clock
	int nucThisRept;		//nucleations in this reporting interval
	int interfaceThick;		//packed integer for ijk interface watch dimensions
	int phaseThick;			//packed integer for ijk phase watch dimensions
	int prevSlush;			//previous slush nodes
	int prevLiquid;			//previous liquid nodes
	int tempThick;			//packed integer for ijk temperature watch dimensions
    int *tempWatch[3];		//array of elements to report thermal info
    int *phaseWatch[3];		//array of elements to report phase info
	int *interfaceWatch[3];	//array of elements to report interface info
    
	string baseName;		//base name
    string dirRoot;		//root directory name-- executable location
    string dirOutput;		//output directory
	string *pathSearch;		//list of directories to search for include files

	double prevStep;		//previous time step     
	double *tInterval;		//reporting range and interval

	CLOCK sClock;			//reporting clock
} REPORT;

typedef struct			// Simulation environment
{
    bool modeRadiation;		// true if radiation effects are modeled
	bool modeResolidify;	// true if stop upon resolidification
	bool modeStochastic;	// true if stochastic effects are active
	bool calcIntraNode;		// true if intra-node calculations were made
	bool calcHeatFlow;		// true if heat-flow calculations were made

	int dTnode;				// maximum temperature change
	int dInode;				// maximum node change
    int curLine;			// current line number for error reporting
    int modeError;			// selects actions to perform on error condition
	int eraMax;				// number of eras during simulation
    int numberNucleated;	// total number of nucleations in simulations
	int numberPhases;		// number of unique phases
	int numberSlush;		// current number of slush nodes
	int numberLiquid;		// current number of liquid nodes
	int seedStochastic;		// seed for rand() generator

	double dTmax;			// maximum temperature change
	double dImax;			// maximum interface change
	double nucThreshold;	// threshold probability for nucleation
	double tEnvironment;	// environmental temperature for radiation

	double *dTemperatureLo;		// lower bound temperature change
	double *dTemperatureHi;		// upper bound temperature change
	double *dInterfaceLo;		// lower bound interface change
	double *dInterfaceHi;		// upper bound interface change
	double *dTimeLo;			// lower bound for clock
	double *dTimeHi;			// upper bound for clock

	double *eraDuration;	// duration of each era

    GENERIC tSubstrate;		// temperature of substrate (file or value)	
	CLOCK sClock;

} SIM;

//_________________________________________________
// public functions
//_________________________________________________

#include "error.h"

//_________________________________________________
// Public variables
//  
//_________________________________________________
#ifdef EXT_LEVEL
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN GEOMETRY Geometry;				// sample geometry
EXTERN INTERFACE Interface;
EXTERN LASER Laser;
EXTERN REPORT Report;
EXTERN SIM Sim;			// time and miscellaneous parameters

EXTERN DMATRIX T;		// Temperature [K] at each node
EXTERN DMATRIX DelX;		//width (in cm) of each node
EXTERN DMATRIX DelY;		//height (in cm) of each node
EXTERN DMATRIX DelZ;		//depth (in cm) of each node
EXTERN DMATRIX DelXY;	//diagonal distance XY
EXTERN DMATRIX DelXZ;	//diagonal distance XZ
EXTERN DMATRIX DelYZ;	//diagonal distance YZ
EXTERN DMATRIX DelXYZ; //diagonal distance XYZ
	
EXTERN DMATRIX AreaXY;	//area in xy plane
EXTERN DMATRIX AreaXZ;	//area in xz plane
EXTERN DMATRIX AreaYZ;	//area in yz plane
EXTERN DMATRIX Volume;	//node volume

EXTERN IMATRIX LiquidCode;	//liquid phase code identifier
EXTERN IMATRIX SolidCode;	//solid phase code identifier
EXTERN PHASE **Phase;  //all phase structures are stored in single array
EXTERN REGION **Region;	 //phase and geometry information on a region of nodes 

#undef EXTERN
#endif