//  _________________________________________________________
// |
// |   Phase.h
// |
// |
// |   (C) 1998  Columbia University, MSME
// |   (C) 2000-01	Harvard University, DEAS
// |__________________________________________________________
#include "3dns.h"

#ifndef _PHASEH
#define _PHASEH

#define ISO_THRESHOLD_SOLID ((double) 2.00)  //prevent 110/111 change if interfaces too close
#define ISO_THRESHOLD_MELT	((double) 2.00)
#define ISO_VELOCITY_110	  ((double) 2.40)  //tuning factor for 110 motion
#define ISO_VELOCITY_111	  ((double) 3.60)  //tuning factor for 111 motion
#define	ISO_LINEAR			//if defined use linear models for 110/111 position/fraction

#define ENTHALPY_MIN  ((double) 0.0)
#define ENTHALPY_MAX  ((double) 1E6)
#define VINTERFACE_MAX ((double) 1E6)		//[cm/sec] maximum solidification velocity
#define VINTERFACE_MIN ((double) -1E6)		//[cm/sec] maximum melting velocity

#define P_FIRST		((int) 1)			//first phase index
#define P_LAST		Sim.numberPhases	//last phase index

#define FULL_LIQUID ((double)0.0)
#define FULL_SOLID ((double)1.0)

#define TMELT transitionTemp[0]

#define LIQ(x) VAL(x, nodeLiquid) 
#define SOL(x) VAL(x, nodeSolid)
#define OTHER(x) VAL(x, nodeOther)
#define SLH(x) VAL(x, nodeSlush)

#define I_EAST	(i < (I_LAST-1) ? i+1 : (I_PERIODIC ? I_FIRST : i))
#define I_WEST	(i > I_FIRST ? i-1 : (I_PERIODIC ? I_LAST-1 : i))
#define J_SOUTH (j < (J_LAST-1) ? j+1 : (J_PERIODIC ? J_FIRST : j))
#define J_NORTH (j > J_FIRST ? j-1 : (J_PERIODIC ? J_LAST-1 : j))
#define K_BACK	(k < (K_LAST-1) ? k+1 : (K_PERIODIC ? K_FIRST : k))
#define K_FRONT	(k > K_FIRST ? k-1 : (K_PERIODIC ? K_LAST-1 : k))

#define BZONE(dir) Interface.bZone[dir]
#define IS_ORDINARY(dir) ((dir & (ALL_DIRORDINARY)) == dir)		//a new direction

#define H_STATE ((int) 0x1)			//instantaneous state of node
#define H_PHASESOL	((int) 0x4)		//instantaneous solid phase code
#define H_GRAINID	((int) 0x1000)		//historical grain code
#define H_DIRSTART	((int) 0x400000)	//historical solidification direction

//_________________________________________________
// Public functions
//_________________________________________________
void InterfaceInit (void);
void InterNode();		//interface motion across boundaries, melting, nucleation
int  IntraNode();		//interface motion inside nodes
void IntraNodeDump();
void IntraNodeFinalize();
void NodeSolidify(CELL &cellNew, CELL &cellOld, const int nodeLiquid, const int dirToLiq, 
				const double tInterface, const int solidIndex, const double timeLeft, const int grainCode);
void NodeMelt(CELL &cellNew, CELL &cellOld, const int nodeA, const int dirToLiq, const double tInterface, const double timeLeft);
void PhaseCleanup (void);
int PhaseFind(int *phaseList, string phaseName);
void PhaseInit (void);
PHASE *PhaseNew(int i);
void PositionInterface (const int dirSlush, const double fracSolid, double &iFrac, double &jFrac, double &kFrac);
int StateSample(void);


//_________________________________________________
// Public variables
//_________________________________________________

#ifdef EXT_LEVEL
#define EXTERN
#else
#define EXTERN extern
#endif

#undef EXT_LEVEL
EXTERN IMATRIX CanChange;		// if node can change
EXTERN BMATRIX CatalyzeMelting;		//if node catalyzes melting of neighbors
EXTERN BMATRIX CatalyzeFreezing;	//if node catalyzes freezing of neighbors
EXTERN DMATRIX FractionSolid;	// fraction of node that is solid
EXTERN IMATRIX GrainCode;		// parent grain code
EXTERN DMATRIX IPos;			// fractional position within node
EXTERN DMATRIX JPos;			// fractional position within node
EXTERN DMATRIX KPos;			// fractional position within node			
EXTERN IMATRIX MaterialClass;	// index for common material classes
EXTERN IMATRIX PhaseCodeLiquid;	// index of liquid phase in each node
EXTERN IMATRIX PhaseCodeSolid;	// index of solid phase in each node
EXTERN IMATRIX PhaseHistory;  	// flags for history of phase evolution
EXTERN DMATRIX QInterface;		// energy change associated with interfacial motion
EXTERN IMATRIX SlushDirection;	// N&S&E&W direction of slush
EXTERN IMATRIX State;		// state of node (SOLID/LIQUID/SLUSH)
EXTERN DMATRIX TInterface;	// temperature of the interface
EXTERN DMATRIX Velocity;	// Velocity (cm/sec)

#undef EXTERN

#endif
