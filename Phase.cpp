//  _________________________________________________________
// |
// |   Phase.cpp    Phase change / interface routines
// |
// |
// |   (C) 1998-99  Columbia University, MSME
// |   (C) 2000-01	Harvard University, DEAS
// |__________________________________________________________

#include "3dns.h"
#include "nucleation.h"
#include "parsefunc.h"		       //for LOOK_UP
#include "report.h"
#include "thermal.h"

#define EXT_LEVEL
#include "phase.h"
#undef EXT_LEVEL

#define DX A(DelX)
#define DY A(DelY)
#define DZ A(DelZ)
#define PHASE_NULL ((int) 0)
 
//___________________________________
// Private functions
//___________________________________
void CalcInterNode(CELL &cellNew, CELL &cellOld);
double CalcIntraNode(CELL &cellNew, CELL &cellOld);
void CalcTInterface(DMATRIX tInterface);
void CheckNodeMelt(CELL &cellNew, CELL &cellOld, const int nodeSolid, const int nodeOther, const int dirSolToLiq, const double tInterface);
void CheckNodeMeltAlien(CELL &cellNew, CELL &cellOld, const int nodeSolid, const int nodeOther, const int dirSolToLiq, const double tInterface);
void CheckNodeMeltSurface(CELL &cellNew, CELL &cellOld, int nodeA, int nodeB, int dirAtoSurf);
void CheckNodeSolidify(CELL &cellNew, CELL &cellOld, const int nodeSolid, const int nodeOther, const int dirSoltoLiq, const double tInterface);
double Fraction100(const double posRaw, const int comp, bool &isComplete);
double Fraction110(const double posRaw, const int comp, bool &isComplete);
double Fraction111(const double posRaw, const int comp, bool &isComplete);
double FractionHigh(const double posRaw, const int comp, bool &isComplete);
void HistorySet(IMATRIX histNew, const int nodeA, const int state, const int phaseSolid, int grainCode, int dirToLiq);
void HistorySet(IMATRIX histNew, const int nodeA, const int state, const int phaseSolid);
double MotionNode (CELL &cellNew, CELL &cellOld, const int i, const int j, const int k);
void NeighborInteract(CELL &cellNew, CELL &cellOld, const int nodeA, const int nodeB, const int dirAtoB);
void NodeMelt(CELL &cellNew, CELL &cellOld, const int nodeA, const int dirToLiq, const double tInterface, const double timeLeft);
void NodeMeltComplete(CELL &cellNew, CELL &cellOld, const int nodeA, double &fracUnused);
void NodeSolidifyComplete(CELL &cellNew, CELL &cellOld, const int nodeA, double &fracUnused);
void NormalAdjust(CELL &cellNew, CELL &cellOld, const int nodeSlush, const int nodeOther,  const int dirPositive, const double tInterface);
inline bool PhaseWillSolidify(const int nodeA, const double tInterface, int &phaseCode);
double Position100(const double fraction, const int component);
double Position110(const double frac, const int comp);
double Position111(const double frac, const int comp);
double PositionHigh(const double frac);

inline double TMelt(const IMATRIX phaseSolid, const int nodeA);
inline int TSolidify(const int nodeA, const int phaseNum);

//___________________________________
// Private variables
//___________________________________
IMATRIX canInteract;		//packed bit directions for interaction
IMATRIX historyNew;
IMATRIX stateNew;			//new state integer
IMATRIX slushDirectionNew;
IMATRIX slushDir[MAX_BZ];	//new slush direction for neighbors in zones (0-2)
IMATRIX pLiquidNew;
IMATRIX pSolidNew;
	
DMATRIX distSlush;		//distance to the nearest interface
DMATRIX distSlushNew;	//distance to the nearest interface (this clock)
DMATRIX fractionNew;	//new fraction solid
DMATRIX timeXInter;		//extra time from inter-node calculations
DMATRIX timeXIntra;		//extra time from intra-node calculations

int oldSlush, oldLiquid;
	
inline double TIMELEFT(double iNew, double iOld)
{	return( (iNew>1) ? (iNew-1.0)/(iNew-iOld) : ((iNew<0) ? iNew/(iNew-iOld) : 0) ); }



//______________________________________________________
//
//	CheckNodeMelt
//		Checks for simple melting in a node with
//		a LIQUID neighbor in the SAME material class
//______________________________________________________
void CheckNodeMelt(CELL &cellNew, CELL &cellOld, const int nodeA, 
		const int nodeOther, const int dirToLiq, const double tInterface)
{
	if (!A(CanChange))			//if node can't melt
		return;

	if ( BZONE(dirToLiq) > 0)		//higher order solidification
		if (A(distSlush) < ISO_THRESHOLD_MELT)	//another interface is too close, or no interfaces around
			return;

	if (tInterface > TMelt(cellOld.cPhaseSolid, nodeA))		 //liquid always will consume solid node of
	{
		NodeMelt(cellNew, cellOld, nodeA, dirToLiq, tInterface, OTHER(cellOld.cTimeX));
		OTHER(Velocity) = 0.0;		//reset velocity in source node
	};

	return;
}; //endfunc


//______________________________________________________
//
//	CheckNodeMeltAlien
//		Checks for simple melting in a solid phase with
//		a neighbor in an ALIEN material class
//______________________________________________________
void CheckNodeMeltAlien(CELL &cellNew, CELL &cellOld, const int nodeSolid, const int nodeOther, const int dirSolToLiq, const double tInterface)
{
	if (!SOL(CanChange))	//if node can't melt
		return;							
							//same material class (RA p.755)
	if (OTHER(CatalyzeMelting) && (tInterface > TMelt(cellOld.cPhaseSolid, nodeSolid)) )
	{
		NodeMelt(cellNew, cellOld, nodeSolid, dirSolToLiq | SURFACEMELT, tInterface, 0.0);
		OTHER(Velocity) = 0.0;
	};

	return;
}; //endfunc


//______________________________________________________
//
//  CheckNodeMeltSurface
//      Moves slush in a node, returns the new fraction
//      solidified
//______________________________________________________
void CheckNodeMeltSurface(CELL &cellNew, CELL &cellOld, int nodeA, int nodeB, int dirAtoSurf)
{
	double tSurf;	//extrapolated surface temperature
	
	if (((A(cellOld.cState) == SOLID) && A(CanChange)) ||
		((A(cellOld.cState) == SLUSH) && IS_ORDINARY(A(cellOld.cDirection)) ))
	{
		tSurf = 1.5 * A(T) - 0.5 * B(T);	//assumes equal size nodes ! change this
	
		if (tSurf > TMelt(cellOld.cPhaseSolid, nodeA))
			NodeMelt(cellNew, cellOld, nodeA, dirAtoSurf | SURFACEMELT, tSurf, 0.0);
	}; //endif

	return;
}; //endfunc


//______________________________________________________
//
//	CheckNodeSolidify
//		Checks for simple solidification in a liquid node
//		with a neighbor in the SAME material class
//______________________________________________________
void CheckNodeSolidify(CELL &cellNew, CELL &cellOld, const int nodeA, const int nodeOther, const int dirToLiq, const double tInterface)
{
	int i;
	int solidIndex = -1;			//index into the table of available solid phases

	double tCritMin = DOUBLE_INFINITY;
	double tCrit;		
	
	if (BZONE(dirToLiq) > 0)		//higher order solidification
		if (A(distSlush) < ISO_THRESHOLD_SOLID)	//another interface is too close, or no interfaces around
			return;

	if (OTHER(CatalyzeFreezing) &&		 // a solid neighbor in same material class
		A(cellOld.cTimeX) <= 0.0)		 // and this liquid node didnt just completely melt
									     // so try find best solid phase to form (RA p.754)
	{	
		for (i = 0; i < Phase[A(cellOld.cPhaseLiquid)]->numTransitions; i++)
		{
			tCrit = Phase[A(cellOld.cPhaseLiquid)]->transitionTemp[i];

			if ((tInterface < tCrit) && (tCrit < tCritMin))	//can solidify to this phase
			{
				tCritMin = tCrit;
				solidIndex = i;
			}; //endif
		}; //endloop i

	} //endif

	else //could still solidify if same solid phase is available and temperature is low enough
	{
		for (i = 0; i < Phase[A(cellOld.cPhaseLiquid)]->numTransitions; i++)
		{
			if (Phase[A(cellOld.cPhaseLiquid)]->transitionTo[i] == OTHER(cellOld.cPhaseSolid) )
			{
				if (tInterface < Phase[A(cellOld.cPhaseLiquid)]->transitionTemp[i])
				{
					solidIndex = i;
					break;
				}; //endif
			}; //endif
		}; //endloop i
	
	}; //endif

	if (solidIndex >= 0)	//found a new phase to form
	{
		NodeSolidify(cellNew, cellOld, nodeA, dirToLiq, tInterface, solidIndex, 
			OTHER(cellOld.cTimeX), OTHER(GrainCode));

		OTHER(Velocity) = 0.0;		//reset velocity in other node
	};

	return;
};


//_________________________________________________
// HistorySet
//		Records historical information on node in packed integer
//
//		Bits [01][23456789AB][CDEF012345][6789ABCD][EF]
//			state  phaseSolid  grainCode  dirStart  unused
//_________________________________________________
void HistorySet(IMATRIX histNew, const int nodeA, const int state, const int phaseSolid, int grainCode, int dirToLiq)
{
	A(histNew) = ((state-1) * H_STATE) + (phaseSolid * H_PHASESOL) +
					(grainCode * H_GRAINID) + (dirToLiq * H_DIRSTART);

	return;
}

void HistorySet(IMATRIX histNew, const int nodeA, const int state, const int phaseSolid)
{
	A(histNew) -= (A(histNew) & (H_GRAINID-1));  //clear instantaneous bits
	A(histNew) += ((state-1) * H_STATE) + (phaseSolid * H_PHASESOL);

	return;
}; 


//__________________________________________________________
// InterNode
//		Inter-node calculations-- always permanent
//__________________________________________________________
void InterNode()
{
	int i, iE, iW;
	int j, jS, jN;
	int k, kB, kF;
	int nodeA;
	int bz, dir;

	CELL cellNew, cellOld;
	
	cellOld.cHistory = PhaseHistory;
	cellOld.cState = State;
	cellOld.cDirection = SlushDirection;
	cellOld.cFraction = FractionSolid;
	cellOld.cTimeX = timeXIntra;	//time left over from intra-node motion
	cellOld.cPhaseLiquid = PhaseCodeLiquid;
	cellOld.cPhaseSolid = PhaseCodeSolid;
	
	MatrixNew(&stateNew);
	MatrixCopy(stateNew, State); //nn: copies the newstate to the current state

	MatrixNew(&slushDirectionNew);
	MatrixCopy(slushDirectionNew, SlushDirection);

	cellNew.cHistory = PhaseHistory;
	cellNew.cState = stateNew;
	cellNew.cDirection = slushDirectionNew;
	cellNew.cFraction = FractionSolid;
	cellNew.cTimeX = timeXInter;	//calculate time to put into neighbors
	cellNew.cPhaseLiquid = PhaseCodeLiquid;
	cellNew.cPhaseSolid = PhaseCodeSolid;

	for (i = 0; i < MAX_BZ; i++)			//new slush directions
		MatrixZero(MatrixNew(slushDir+i));

	MatrixZero(MatrixNew(&distSlushNew));
	MatrixZero(timeXInter);		//clear any residual values

	Report.nucThisClock = 0;	//reset nucleation counter

	//___________ STANDARD NEIGHBOR NODE INTERACTIONS ____________
	for (k=K_FIRST; kB=K_BACK, kF=K_FRONT, k < K_LAST; k++)	
	{
		if ( !(Geometry.canChangeK[k] || Geometry.canChangeK[kB] ||
				Geometry.canChangeK[kF]))
			continue;		//these 3 pages cant change

		for (j=J_FIRST; jS = J_SOUTH, jN = J_NORTH, j < J_LAST; j++)
		{
			if ( !(Geometry.canChangeJ[j] || Geometry.canChangeJ[jS] ||
				Geometry.canChangeJ[jN]))
				continue;		//these 3 columns cant change so skip
			
			for (i=I_FIRST; iE=I_EAST, iW=I_WEST, i < I_LAST; i++)
			{	
				if ( !(Geometry.canChangeI[i] || Geometry.canChangeI[iE] ||
					Geometry.canChangeI[iW]))
					continue;		//these 3 rows cant change so skip

				nodeA = IJKToIndex(i, j, k);
			
				//_____ (111) interactions ________
				NeighborInteract(cellNew, cellOld, nodeA, IJKToIndex(iE, jS, kB), EAST|SOUTH|BACK);
				NeighborInteract(cellNew, cellOld, nodeA, IJKToIndex(iE, jS, kF), EAST|SOUTH|FRONT);
				NeighborInteract(cellNew, cellOld, nodeA, IJKToIndex(iE, jN, kB), EAST|NORTH|BACK);
				NeighborInteract(cellNew, cellOld, nodeA, IJKToIndex(iW, jS, kB), WEST|SOUTH|BACK);

				//_____ (110) interactions ________
				NeighborInteract(cellNew, cellOld, nodeA, IJKToIndex(iE, jS, k), EAST|SOUTH);
				NeighborInteract(cellNew, cellOld, nodeA, IJKToIndex(iE, j, kB), EAST|BACK);
				NeighborInteract(cellNew, cellOld, nodeA, IJKToIndex(iE, jN, k), EAST|NORTH);
				NeighborInteract(cellNew, cellOld, nodeA, IJKToIndex(iW, j, kB), WEST|BACK);
				NeighborInteract(cellNew, cellOld, nodeA, IJKToIndex(i, jS, kB), SOUTH|BACK);
				NeighborInteract(cellNew, cellOld, nodeA, IJKToIndex(i, jS, kF), SOUTH|FRONT);
				
				//_____ (100) interactions ________ 
				NeighborInteract(cellNew, cellOld, nodeA, IJKToIndex(iE, j, k), EAST);
				NeighborInteract(cellNew, cellOld, nodeA, IJKToIndex(i, jS, k), SOUTH);
				NeighborInteract(cellNew, cellOld, nodeA, IJKToIndex(i, j, kB), BACK);

			}; //endloop i
		}; //endloop j
	}; //endloop k
	
	NucleateHomogeneous(cellNew, cellOld);		//Homogeneous nucleation

	//___________ SPECIAL EDGE CONDITIONS ____________
	if (Geometry.canChangeJ[J_FIRST])	
	{	
		j = J_FIRST;

		for (i = I_FIRST; i < I_LAST; i++)			//top surface melting
			for (k = K_FIRST; k < K_LAST; k++)
				CheckNodeMeltSurface(cellNew, cellOld, IJKToIndex(i, j, k), 
						IJKToIndex(i, J_SOUTH, k), NORTH);
	}; //endif

	//___________ RECALCULATE SLUSH POSITIONS ____________
	for (i = I_FIRST; i < I_LAST; i++)
	{
		if (! Geometry.canChangeI[i])
			continue;

		for (j = J_FIRST; j < J_LAST; j++)
		{			
			if (! Geometry.canChangeJ[j])
				continue;

			for (k = K_FIRST; k < K_LAST; k++)
			{
				nodeA = IJKToIndex(i, j, k);
					
				if (A(cellNew.cState) == SLUSH)
				{
					for (bz = MAX_BZ-1, dir=0; bz >= 0; bz--)	//combine all zone contributions
					{
						dir &= Interface.dirKeep[bz][A(slushDir[bz])];	//clear some lesser zone bits
						dir |= A(slushDir[bz]);		//add directions from this zone
					}; //endloop bz
  					
					if (!IS_ORDINARY(A(cellNew.cDirection))) //previously catalyzed
						A(cellNew.cDirection) |= dir;		//only add directions
					
					else      //recalculate all slush nodes
						A(cellNew.cDirection) = dir;	//update direction


					if ((A(cellOld.cState) != SLUSH) ||	//new slush or changed direction
						(cellNew.cDirection != cellOld.cDirection) )

						PositionInterface(A(cellNew.cDirection), A(cellNew.cFraction), 
							A(IPos), A(JPos), A(KPos) );	//update positions

					A(TInterface) = InterpolateT (i, j, k, A(IPos), A(JPos), A(KPos));
												//always will be new temperatures
				}; //endif

			}; //endloop k
		}; //endloop j
	}; //endloop i

	for (i = 0; i < MAX_BZ; i++)
		MatrixRecycle(slushDir[i]);			//dump old slush directions

	MatrixRecycle(distSlush);
	distSlush = distSlushNew;

	memset(&cellNew, 0, sizeof(CELL)); 
	memset(&cellOld, 0, sizeof(CELL));
		//workaround-- optimizing compiler seems to mess with FractionSolid
		//if we let it destruct cellNew naturally.								

	MatrixRecycle(State);				//make state changes permanent	
	State = stateNew;

	MatrixRecycle(SlushDirection);		//slush direction changes are permanent
	SlushDirection = slushDirectionNew;

	return;
}; //endfunc


//__________________________________________________________
// IntraNode
//		Calculate node states based on conditions within
//		individual nodes
//__________________________________________________________
int IntraNode()
{	
	static bool zeroedQT=false;  //true if matrices are already zeroed

	int i,j,k;

	double dF;		//maximum fractional change

	CELL cellOld, cellNew;

	if (Sim.dImax > 0.0)		//must have been some movement previously
	{
		MatrixZero(QInterface);
		MatrixZero(timeXIntra);
	}; //endif

	Sim.dImax = 0;		//reset maximum interface motion counter

	if (Sim.numberSlush == 0)	//nothing to do
		return(CHANGE_OK);

	Sim.calcIntraNode = true;

	oldSlush = Sim.numberSlush;		//save number of slush nodes
	oldLiquid = Sim.numberLiquid;	//save number of liquid nodes

	cellOld.cHistory = PhaseHistory;
	cellOld.cState = State;
	cellOld.cDirection = SlushDirection;
	cellOld.cFraction = FractionSolid;
	cellOld.cTimeX = timeXInter;	//time left over from intra-node motion
	cellOld.cPhaseLiquid = PhaseCodeLiquid;
	cellOld.cPhaseSolid = PhaseCodeSolid;


	MatrixNew(&historyNew);
	MatrixCopy(historyNew, PhaseHistory);

	MatrixNew(&stateNew);
	MatrixCopy(stateNew, State);

	MatrixNew(&slushDirectionNew);
	MatrixCopy(slushDirectionNew, SlushDirection);

	MatrixNew(&fractionNew);
	MatrixCopy(fractionNew, FractionSolid);

	MatrixNew(&pSolidNew);
	MatrixCopy(pSolidNew, PhaseCodeSolid);

	cellNew.cHistory = historyNew;
	cellNew.cState = stateNew;
	cellNew.cDirection = slushDirectionNew;
	cellNew.cFraction = fractionNew;
	cellNew.cTimeX = timeXIntra;		//new extra times from intra-node motion
	cellNew.cPhaseLiquid = PhaseCodeLiquid;
	cellNew.cPhaseSolid = pSolidNew;
	
	//___________ INTRA-NODE MOTION ____________
	for (i = I_FIRST; i < I_LAST; i++)
	{
		if (!Geometry.canChangeI[i])
			continue;

		for (j = J_FIRST; j < J_LAST; j++)
		{
			if (!Geometry.canChangeJ[j])
				continue;

			for (k = K_FIRST; k < K_LAST; k++)
			{
				if ((cellOld.cState)[i][j][k] == SLUSH)
				{
					dF = MotionNode(cellNew, cellOld, i, j, k);
					
					if (dF > Sim.dImax)		//find maximum node motion
					{
						Sim.dImax = dF;
						Sim.dInode = IJKToIndex(i, j, k);
					}; //endif
				}; //endif
			}; //endloop k
		}; //endloop j
	}; //endloop i

	memset(&cellOld, 0, sizeof(CELL));		//workaround !!
	memset(&cellNew, 0, sizeof(CELL));

	if (Sim.dImax > Sim.dInterfaceHi[THIS_ERA])	//too much change, throw away calculations
		throw(CHANGE_HI);
	
	return(CHANGE_OK);

}; //endfunc


// ______________________________________________________________
// IntraNodeDump
//		Dump interface calculations, free up temporary matrices
//		works even if matrices are not allocated
// ______________________________________________________________
void IntraNodeDump()
{
	if (Sim.calcIntraNode)
	{
		MatrixRecycle(historyNew);
		MatrixRecycle(fractionNew);
		MatrixRecycle(slushDirectionNew);
		MatrixZero(timeXIntra);		//dump any intranode motion extra times
		MatrixZero(timeXInter);		//these extra values were from larger clock
		MatrixRecycle(pSolidNew);
		MatrixRecycle(stateNew);

		Sim.numberSlush = oldSlush;
		Sim.numberLiquid = oldLiquid;
	
		Sim.calcIntraNode = false;
	}; //endif

	return;
}; //endfunc


// ______________________________________________________________
// IntraNodeFinalize
//      Finalize interface states
// ______________________________________________________________
void IntraNodeFinalize()
{
	if (Sim.calcIntraNode)
	{
		MatrixRecycle(PhaseHistory);		//change is okay, keep calculations
		PhaseHistory = historyNew;

		MatrixRecycle(State);
		State = stateNew;

		MatrixRecycle(FractionSolid);
		FractionSolid = fractionNew;

		MatrixRecycle(SlushDirection);
		SlushDirection = slushDirectionNew;

		//(keep timeXIntra)

		MatrixRecycle(PhaseCodeSolid);
		PhaseCodeSolid = pSolidNew;

		Sim.calcIntraNode = false;		//reset for more calculations
	}; //endif

	return;
}; //endfunc


// ______________________________________________________________
// InterfaceInit
//      Interface tracking model parameters
// ______________________________________________________________
void InterfaceInit (void)
{
	#define PSEUDO(x,y) ((i&x)>0) ? (((i&y)>0) ? 2 : 1) : (((i&y)>0) ? -1 : 0)
	
	int i, invDir, dir0, dir1;

	const int nonOrdinary = (ALL_DIRECTIONS) ^ (ALL_DIRORDINARY);

	#define DIMS ((int)Geometry.isThick[DIM_I] + (int)Geometry.isThick[DIM_J] +(int)Geometry.isThick[DIM_K])

	Interface.numNeighbors[0] = 2 * DIMS;
	Interface.numNeighbors[1] = (DIMS > 1) ? (DIMS * 8) - 12 : 0;
	Interface.numNeighbors[2] = (DIMS > 2) ? 8 : 0;

	for (i = 0; i <= ALL_DIRECTIONS; i++)
	{
		Interface.iCrystal[i] = PSEUDO(EAST, WEST);		//pseudo-crystallographic components
		Interface.jCrystal[i] = PSEUDO(SOUTH, NORTH);
		Interface.kCrystal[i] = PSEUDO(BACK, FRONT);

		Interface.crystalType[i] = abs(Interface.iCrystal[i] * 256) + 
			abs(Interface.jCrystal[i] * 16) + abs(Interface.kCrystal[i]);

		invDir = 0;
		invDir |= ((i & EAST) ? WEST : 0);		//inverse components
		invDir |= ((i & WEST) ? EAST : 0);
		invDir |= ((i & SOUTH) ? NORTH : 0);
		invDir |= ((i & NORTH) ? SOUTH : 0);
		invDir |= ((i & BACK) ? FRONT : 0);
		invDir |= ((i & FRONT) ? BACK : 0);

		Interface.dirInverse[i] = invDir;

		dir0 = dir1 = 0;		//save nucleation specifiers
		
		dir0 |= ((i & (EAST|WEST)) ? (NORTH|SOUTH|FRONT|BACK) : 0); //bits to dump
		dir0 |= ((i & (SOUTH|NORTH)) ? (EAST|WEST|FRONT|BACK) : 0); //bits to dump
		dir0 |= ((i & (FRONT|BACK)) ? (NORTH|SOUTH|EAST|WEST) : 0); //bits to dump

		Interface.dirKeep[0][i] = (i > 0) ?
			((dir0 ^ (ALL_DIRECTIONS)) | nonOrdinary) : ALL_DIRECTIONS;

		dir1 |= ((i & (EAST|WEST)) ? (EAST|WEST) : 0);		//bits to save
		dir1 |= ((i & (SOUTH|NORTH)) ? (SOUTH|NORTH) : 0);
		dir1 |= ((i & (FRONT|BACK)) ? (FRONT|BACK) : 0);
	
		Interface.dirKeep[1][i] = (i > 0) ? (dir1 | nonOrdinary) : ALL_DIRECTIONS;
		
		Interface.dirKeep[2][i] = 0;	//no zones beyond 2!

		switch (Interface.crystalType[i])
		{
			case ((int)(0x111)):
				Interface.bZone[i] = 2;
				break;

			case (0x110):
			case (0x101):
			case (0x011):
				Interface.bZone[i] = 1;
				break;

			case ((int)(0x100)):
			case ((int)(0x010)):
			case ((int)(0x001)):
				Interface.bZone[i] = 0;
				break;

			default:
				Interface.bZone[i] = 0;
		}; //endswitch
	}; //endloop i

	return;
}; //endfunc


//______________________________________________________
//
//  MotionNode
//      Moves slush in a node, returns the new fraction
//      solidified.  Assumes that SlushDirection, State,
//		and IJK Positions are already resolved for this
//		clock.
//______________________________________________________
double MotionNode(CELL &cellNew, CELL &cellOld, const int i, const int j, const int k)
{
	int nodeA = IJKToIndex(i, j, k);	//index into arrays //note that #def A(x) VAL(x,nodeA) 
	bool isComplete;			//complete melting or solidification

	double clockUnused = 0.0;	//left over clock
	double delPos;				//change in position
	double fractionOld, posOld, posNew;		//old fraction solid and position
	double f=0;
	
								//assumes TInterface, IPos, JPos, KPos are up-to-date from InterNode
    A(Velocity) = LOOK_UP (Phase[A(cellNew.cPhaseSolid)]->interfaceResponse, A(TInterface)); //calculate velocity at i,j,k

	fractionOld = A(cellOld.cFraction);

	switch ( Interface.crystalType[A(cellNew.cDirection)] )		// uses (ijk) crystallographic notation
	{
		#define ICOMP Interface.iCrystal[A(cellNew.cDirection)]
		#define JCOMP Interface.jCrystal[A(cellNew.cDirection)]
		#define KCOMP Interface.kCrystal[A(cellNew.cDirection)]

		case (0x100):
			posOld = A(IPos);
			posNew = posOld + (Sim.sClock.curDTime + A(cellOld.cTimeX)) * A(Velocity) * ICOMP / DX;
			A(cellNew.cFraction) = Fraction100(posNew, ICOMP, isComplete);
			break;

		case (0x010):
			posOld = A(JPos);
			posNew = posOld + (Sim.sClock.curDTime + A(cellOld.cTimeX)) * A(Velocity) * JCOMP / DY;
			A(cellNew.cFraction) = Fraction100(posNew, JCOMP, isComplete);
			break;

		case (0x001):
			posOld = A(KPos);
			posNew = posOld + (Sim.sClock.curDTime + A(cellOld.cTimeX)) * A(Velocity) * KCOMP / DZ;
			A(cellNew.cFraction) = Fraction100(posNew, KCOMP, isComplete);
			break;

		case (0x110):
			posOld = A(IPos);

			delPos = (Sim.sClock.curDTime + A(cellOld.cTimeX)) * A(Velocity) * ISO_VELOCITY_110 * A(DelXY) / 
				(2 * A(AreaXY));
			
			posNew = posOld + delPos * ICOMP;
			A(cellNew.cFraction) = Fraction110(posNew, ICOMP, isComplete);
			break;

		case (0x101):
			posOld = A(IPos);

			delPos = (Sim.sClock.curDTime + A(cellOld.cTimeX)) * A(Velocity) * ISO_VELOCITY_110 * A(DelXZ) / 
				(2 * A(AreaXZ));

			posNew = posOld + delPos * ICOMP;
			A(cellNew.cFraction) = Fraction110(posNew, ICOMP, isComplete);
			break;

		case (0x011):
			posOld = A(JPos);

			delPos = (Sim.sClock.curDTime + A(cellOld.cTimeX)) * A(Velocity) * ISO_VELOCITY_110 * A(DelYZ) / 
				(2 * A(AreaYZ));

			posNew = posOld + delPos * JCOMP;
			A(cellNew.cFraction) = Fraction110(posNew, JCOMP, isComplete);
			break;

		case (0x111):
			posOld = A(IPos);
			
			delPos = (Sim.sClock.curDTime + A(cellOld.cTimeX)) * A(Velocity) * ISO_VELOCITY_111 *
				sqrt (DX*DX*DY*DY + DX*DX*DZ*DZ + DY*DY*DZ*DZ) / 
				(3 * DX*DY*DZ);
			
			posNew = posOld + delPos * ICOMP;
			A(cellNew.cFraction) = Fraction111(posNew, ICOMP, isComplete);
			break;

		default:		//use 3D scaling
			posOld = PositionHigh(A(cellOld.cFraction));	//position=fraction here
			posNew = posOld + 2.0 * (Sim.sClock.curDTime + A(cellOld.cTimeX)) * A(Velocity) *
				(1/DX + 1/DY + 1/DZ);
			A(cellNew.cFraction) = FractionHigh(posNew, ICOMP, isComplete);														//fraction = posNew
			break;

		}; //endswitch

	clockUnused = TIMELEFT(posNew, posOld);		//fraction of clock unused

	if (isComplete && (A(cellNew.cFraction) > (FULL_SOLID - 0.01)))
		NodeSolidifyComplete(cellNew, cellOld, nodeA, clockUnused);

	else if (isComplete && (A(cellNew.cFraction) < (0.01 + FULL_LIQUID)))
		NodeMeltComplete(cellNew, cellOld, nodeA, clockUnused);

    A(QInterface) = (A(cellNew.cFraction) - fractionOld) * A(Volume) *
    	LOOK_UP(Phase[A(cellNew.cPhaseSolid)]->enthalpyH, A(TInterface) );		//energy associated with move

    return( fabs(posNew - posOld) );  //change in fractionSolid + extra
}; //endfunc


//__________________________________________________________
// NeighborInteract
//		Calculate behavior at boundary of neighboring nodes
//		initiate solidification or melting if necessary
//__________________________________________________________
void NeighborInteract(CELL &cellNew, CELL &cellOld, const int nodeA, const int nodeB, const int dirAtoB)
{	
	double tInterface;
	int nodeLiquid, nodeSolid, nodeSlush;	//node indices

	#define MATCH_STATE(x) ((A(State) == x) ? nodeA : ((B(State) == x) ? nodeB : -1 ))
	#define OTHER_NODE(x)  ((x == nodeA) ? nodeB : nodeA)
	#define DIR(x, y)  ((x == nodeA) ? dirAtoB : Interface.dirInverse[dirAtoB])

	if ((A(canInteract) & dirAtoB) != dirAtoB)		//no interaction in this direction
		return;

	//_______________ SAME MATERIAL CLASS _____________________
	if ( A(MaterialClass) == B(MaterialClass) )
	{
		switch (A(cellOld.cState) | B(cellOld.cState))		//see interaction table RA p.750
		{
			case (SOLID | SOLID):
				if (A(cellOld.cPhaseSolid) == B(cellOld.cPhaseSolid))
					return;						//exactly the same phase

				tInterface = InterpolateT(nodeA, dirAtoB);	//interface temperature

				CheckNodeMelt(cellNew, cellOld, nodeA, nodeB, DIR(nodeA, nodeB), tInterface);
				CheckNodeMelt(cellNew, cellOld, nodeB, nodeA, DIR(nodeB, nodeA), tInterface);
	
				break;

			case (LIQUID | LIQUID):
				return;
				
				break;

			case (SOLID | LIQUID): 
				nodeSolid = MATCH_STATE(SOLID);				
				nodeLiquid = OTHER_NODE(nodeSolid);

				tInterface = InterpolateT(nodeA, dirAtoB);	//interface temperature

				CheckNodeMelt(cellNew, cellOld, nodeSolid, nodeLiquid, DIR(nodeSolid, nodeLiquid), tInterface );
				CheckNodeSolidify(cellNew, cellOld, nodeLiquid, nodeSolid, DIR(nodeSolid, nodeLiquid), tInterface);

				break;

			case (SOLID | SLUSH):
				nodeSolid = MATCH_STATE(SOLID);
				nodeSlush = OTHER_NODE(nodeSolid);

				NormalAdjust(cellNew, cellOld, nodeSlush, nodeSolid, DIR(nodeSolid, nodeSlush), SLH(TInterface));
				break;

			case (LIQUID | SLUSH):				
				nodeLiquid = MATCH_STATE(LIQUID);
				nodeSlush = OTHER_NODE(nodeLiquid);

				NormalAdjust(cellNew, cellOld, nodeSlush, nodeLiquid, DIR(nodeSlush, nodeLiquid), SLH(TInterface));
				break;

			case (SLUSH | SLUSH):
				break;

		}; //endswitch
	}

	//_______________ ALIEN MATERIAL CLASS _____________________
	else if ( BZONE(dirAtoB) == 0)   //zone 0 interactions only
	{
		switch (A(State) | B(State))		//see interaction table RA p.750
		{
			case (SOLID | SOLID):
				tInterface = InterpolateT(nodeA, dirAtoB);	//interface temperature

				CheckNodeMeltAlien(cellNew, cellOld, nodeA, nodeB, DIR(nodeA, nodeB), tInterface);
				CheckNodeMeltAlien(cellNew, cellOld, nodeB, nodeA, DIR(nodeB, nodeA), tInterface);

				break;

			case (LIQUID | LIQUID):
				return;
				
				break;

			case (SOLID | LIQUID):
				nodeSolid = MATCH_STATE(SOLID);				
				nodeLiquid = OTHER_NODE(nodeSolid);
				
				tInterface = InterpolateT(nodeA, dirAtoB);	 //interface temperature

				NucleateHeterogeneous(cellNew, cellOld, nodeLiquid, nodeSolid, DIR(nodeSolid, nodeLiquid), tInterface);
				CheckNodeMeltAlien(cellNew, cellOld, nodeSolid, nodeLiquid, DIR(nodeSolid, nodeLiquid), tInterface);
				break;

			case (SOLID | SLUSH):	//no effect from alien interfaces
				break;

			case (LIQUID | SLUSH):	//no effect from alien interfaces			
				break;

			case (SLUSH | SLUSH):	//no effect from alien interfaces
				break;

		}; //endswitch

	}; //endif

	return;
}; //endfunc


//_________________________________________________
// NodeMelt
//		Initiates melting in nodeA from neighbor nodeB.
//_________________________________________________
void NodeMelt(CELL &cellNew, CELL &cellOld, const int nodeA, const int dirToLiq, const double tInterface, const double timeLeft)
{
	if (A(cellNew.cState) == SOLID)		//initial melting
	{
		Sim.numberSlush++;
		HistorySet(cellNew.cHistory, nodeA, SLUSH, A(cellOld.cPhaseSolid));
	}; //endif

	A(cellNew.cState) = SLUSH;
	A(slushDir[ BZONE(dirToLiq) ]) = dirToLiq;
		
	A(cellNew.cPhaseLiquid) = *(Phase[A(cellOld.cPhaseSolid)]->transitionTo);

	if ( BZONE(dirToLiq) == 0)		
		A(cellNew.cTimeX) += timeLeft;		//extra time from zone-0 neighbors	

	return;
}; //endfunc


//_________________________________________________
// NodeMeltComplete
//		Complete solidification of nodeA.
//		Always occurs after (stateNew --> State) update
//_________________________________________________
void NodeMeltComplete(CELL &cellNew, CELL &cellOld, const int nodeA, double &fracUnused)
{
	A(cellNew.cTimeX) = Sim.sClock.curDTime * fracUnused;

    A(cellNew.cState) = LIQUID;
	A(cellNew.cDirection) = 0;
	A(FractionSolid) = 0.;
    Sim.numberSlush--;
	Sim.numberLiquid++;

	HistorySet(cellNew.cHistory, nodeA, LIQUID, 0, 0, 0);	//erase solid information

	return;
}; //endfunc


//_________________________________________________
// NodeSolidify
//		Initiates solidification in nodeA from neighbor nodeB.
//		entry point from nucleation
//_________________________________________________
void NodeSolidify(CELL &cellNew, CELL &cellOld, const int nodeLiquid, const int dirToLiq, 
				const double tInterface, const int solidIndex, const double timeLeft, const int grainCode)
{
	if (LIQ(cellOld.cState) != LIQUID)		//old state MUST be liquid
		ErrorMsg("Cant treat non-liquid node in NodeSolidify at time %g", Sim.sClock.curTime);

	LIQ(cellNew.cPhaseSolid) = *(Phase[LIQ(cellOld.cPhaseLiquid)]->transitionTo + solidIndex);	

	if (LIQ(cellNew.cState) == LIQUID)	//first solidification direction 
	{
		Sim.numberSlush++;
		Sim.numberLiquid--;
		
		LIQ(GrainCode) = grainCode;

		HistorySet(cellNew.cHistory, nodeLiquid, SLUSH, LIQ(cellNew.cPhaseSolid), grainCode, dirToLiq);
 	}; //endif
	
	LIQ(cellNew.cState) = SLUSH;
	LIQ(slushDir[ BZONE(dirToLiq) ]) |= dirToLiq;
	LIQ(TInterface) = tInterface;	//update tInterface

	if ( BZONE(dirToLiq) == 0)		
		LIQ(cellNew.cTimeX) += timeLeft;		//extra time from zone-0 neighbors

	return;
}; //endfunc


//_________________________________________________
// NodeSolidifyComplete
//		Complete solidification of nodeA.
//		Always occurs after (stateNew-->State) update
//_________________________________________________
void NodeSolidifyComplete(CELL &cellNew, CELL &cellOld, const int nodeA, double &fracUnused)
{        
	A(cellNew.cTimeX) = Sim.sClock.curDTime * fracUnused;

	A(cellNew.cState) = SOLID;
	A(cellNew.cDirection) = 0;
	A(FractionSolid) = 1.;
	HistorySet(cellNew.cHistory, nodeA, SOLID, A(cellOld.cPhaseSolid));
	Sim.numberSlush--;

	return;
}; //endfunc


//______________________________________________________
//
//	NormalAdjust
//		Change direction of interface based on nearest
//		neighbor states
//______________________________________________________
void NormalAdjust(CELL &cellNew, CELL &cellOld, const int nodeSlush, const int nodeOther, const int dirPositive, const double tInterface)
{
	static double timeX[3] = {0.25, 0.125, 0.0625};	//diffusion of extra time to neighbors
	static double dist, di, dj, dk;
	static int dirDist;		//direction to calculate distance (in positive slush direction)
	
	dirDist = (OTHER(cellOld.cState) == LIQUID) ? dirPositive : Interface.dirInverse[dirPositive];

	di = (dirDist & WEST) ? SLH(IPos) : (dirDist & EAST) ? (1.0 - SLH(IPos)) : 0;
	dj = (dirDist & NORTH) ? SLH(JPos) : (dirDist & SOUTH) ? (1.0 - SLH(JPos)) : 0;
	dk = (dirDist & FRONT) ? SLH(KPos) : (dirDist & BACK) ? (1.0 - SLH(KPos)) : 0;

	dist = sqrt(di*di + dj*dj + dk*dk);		//distance to neighbor node boundary

	if ((OTHER(distSlushNew) > dist) || (OTHER(distSlushNew) == 0.0))
		OTHER(distSlushNew) = dist;			//keep smallest distance

	if ((OTHER(cellOld.cState) == SOLID) && (tInterface > TMelt(cellOld.cPhaseSolid, nodeSlush)))
		return;
	if ((OTHER(cellOld.cState) == LIQUID) && (tInterface < TMelt(cellOld.cPhaseSolid, nodeSlush)))
		return;

	SLH(slushDir[BZONE(dirPositive)]) |= dirPositive;		//turn on this direction	

	if(BZONE(dirPositive) == 0)
		SLH(cellNew.cTimeX) += timeX[ BZONE(dirPositive) ] * OTHER(cellOld.cTimeX);  
								//get extra time if any
	return;
}; //endfunc


//_________________________________________________
// PhaseCleanup
//_________________________________________________
void PhaseCleanup ()
{
    MatrixFree (canInteract);
	MatrixFree (CanChange);
    MatrixFree (FractionSolid);
    MatrixFree (PhaseCodeSolid);
	MatrixFree (PhaseCodeLiquid);
    MatrixFree (PhaseHistory);
    MatrixFree (QInterface);
    MatrixFree (SlushDirection);
    MatrixFree (TInterface);
    MatrixFree (State);
    MatrixFree (Velocity);

}; //endfunc


//_________________________________________________
// PhaseFind
//_________________________________________________
int PhaseFind(int *phaseList, string phaseName)
{
	int i = 0;
	
	if (phaseList[i] == 0)		//no list to search
		return(0);

	while (phaseList[i] > 0)
	{
		if (Phase[phaseList[i]] == NULL)
			ErrorMsg("Null phase in PhaseFind");

		if (phaseName == Phase[phaseList[i]]->phaseName)
			return(phaseList[i]);

		i++;
	}; //endwhile

	return(0);		//null-phase
}; //endfunc


//_________________________________________________
// PhaseInit
//_________________________________________________
void PhaseInit ()
{
	char *typeMsg;
	int typeVal;

	int p, r;
	int i, j, k;	//pointers to region positions
	int nodeA;

    MatrixNew (&CanChange);
    MatrixNew (&FractionSolid);
	MatrixZero(MatrixNew(&GrainCode));
	MatrixNew (&IPos);
	MatrixNew (&JPos);
	MatrixNew (&KPos);
	MatrixNew (&MaterialClass);
    MatrixNew (&PhaseCodeLiquid);
	MatrixNew (&PhaseCodeSolid);
    MatrixZero(MatrixNew(&PhaseHistory));
    MatrixZero(MatrixNew(&QInterface));
    MatrixZero(MatrixNew(&SlushDirection));
    MatrixNew (&TInterface);
    MatrixNew (&State);
    MatrixNew (&Velocity);

	MatrixZero(MatrixNew(&canInteract));
	MatrixZero(MatrixNew(&distSlush));
	
	MatrixZero(MatrixNew(&timeXInter));	//extra times from internode motion
	MatrixZero(MatrixNew(&timeXIntra));	//extra times from intranode motion

	MatrixNew(&CatalyzeFreezing);
	MatrixNew(&CatalyzeMelting);

	Geometry.canChangeI = new bool[I_LAST];					//if columns can change
	Geometry.canChangeJ = new bool[J_LAST];					//if rows can change
	Geometry.canChangeK = new bool[K_LAST];					//if pages can change

	if ((Geometry.canChangeI == NULL) || (Geometry.canChangeJ == NULL) ||
			(Geometry.canChangeK == NULL) )
		ErrorMsg("Out of memory");		

	memset(Geometry.canChangeI, 0, sizeof(bool)*I_LAST);
	memset(Geometry.canChangeJ, 0, sizeof(bool)*J_LAST);
	memset(Geometry.canChangeK, 0, sizeof(bool)*K_LAST);

   	if (!InRange( Sim.dInterfaceLo, 0, Sim.eraMax, 0, MAX_ICHANGE))
    	ErrorMsg("INTERFACE_MOTION_LO is out of range, or not specified");	

   	if (!InRange( Sim.dInterfaceHi, 0, Sim.eraMax, 0, MAX_ICHANGE))
    	ErrorMsg("INTERFACE_MOTION_HI is out of range, or not specified");	

	PhaseNew(0);	//initialize NULL phase

	for (p=1; p<Sim.numberPhases; p++)	//first phase is '1', while '0' is NULL phase
	{
		if (Phase[p] == NULL)
			ErrorMsg("Phase %d not defined", p);

		if (!InRange(Phase[p]->phaseType, 0, MAX_PHASES))
			ErrorMsg("Bad phase type in phase %d", p);

		if (Phase[p]->numTransitions > 0)	// this phase can change
		{
			if (!InRange(Phase[p]->transitionTo, 0, Phase[p]->numTransitions, 1, Sim.numberPhases))
				ErrorMsg ("Cant find transition phase index in %s %s", 
					Phase[p]->matlClassName, Phase[p]->phaseName);

			if (!InRange(Phase[p]->transitionTemp, 0, Phase[p]->numTransitions, MIN_DEGREES, MAX_DEGREES))
				ErrorMsg ("Transition temperature out of range in %s %s", 
					Phase[p]->matlClassName, Phase[p]->phaseName);
		
			if (Phase[p]->phaseType != LIQUID)		//solid nodes
			{			
				if (!InRange (Phase[p]->enthalpyH, MIN_DEGREES, MAX_DEGREES, ENTHALPY_MIN, ENTHALPY_MAX))
					ErrorMsg ("Enthalpy out of range {%g,%g} in %s %s",
        				ENTHALPY_MIN, ENTHALPY_MAX, Phase[p]->matlClassName, Phase[p]->phaseName);
			
				if (!InRange (Phase[p]->interfaceResponse,
						MIN_DEGREES, (int)(Phase[p]->TMELT), 0.0, VINTERFACE_MAX))
					ErrorMsg ("Interface response out of range for T<Tm in %s %s", 
						Phase[p]->matlClassName, Phase[p]->phaseName);

				if (!InRange (Phase[p]->interfaceResponse,
						(int)(Phase[p]->TMELT)+1, MAX_DEGREES, VINTERFACE_MIN, 0.0))
					ErrorMsg ("Interface response out of range for T>Tm in %s %s", 
						Phase[p]->matlClassName, Phase[p]->phaseName);

				if (!InRange( LOOK_UP(Phase[p]->interfaceResponse, Phase[p]->TMELT), -1E-9, 1E-9))
					ErrorMsg ("Interface response function must be 0 at transition temperature in %s %s",
    					Phase[p]->matlClassName, Phase[p]->phaseName);
			}; //endif	
		}; //endif
	}; //endloop p
		
	//_________ CHECK REGION INFORMATION _______________
	for (r=0; r<Geometry.numberRegions; r++)
	{
		typeMsg = (r < Geometry.jZones) ? "LAYER" : "OVERLAY";
		typeVal = (r < Geometry.jZones) ? r : r - Geometry.jZones;
				
		if (Region[r] == NULL)
    		ErrorMsg("Bad %s number %d, check .dat file", typeMsg, typeVal);

		if (!InRange(Region[r]->phaseStart, 1, Sim.numberPhases))
			ErrorMsg("Bad starting phase in %s %d", typeMsg, typeVal);

		if (Region[r]->canChange)
		{
			if (!InRange(Region[r]->phaseLiquid, 1, Sim.numberPhases) ||
				(Phase[Region[r]->phaseLiquid]->phaseType != LIQUID) )
				ErrorMsg("No liquid phase avaiable in  %s %d", typeMsg, typeVal);
		}; //endif

		for (i = 0; i < Region[r]->numberI; i++)
		{
			for (j = 0; j < Region[r]->numberJ; j++)
			{
				for (k = 0; k < Region[r]->numberK; k++)
				{
					nodeA = IJKToIndex(Region[r]->iLocations[i], 
							Region[r]->jLocations[j], Region[r]->kLocations[k]);

					A(CanChange) = Region[r]->canChange;
					A(CatalyzeFreezing) = Region[r]->catalyzeFreezing;
					A(CatalyzeMelting) = Region[r]->catalyzeMelting;
					A(MaterialClass) = Phase[Region[r]->phaseStart]->matlClass;

					Geometry.canChangeI[Region[r]->iLocations[i]] =
						Geometry.canChangeI[Region[r]->iLocations[i]] || A(CanChange);
					Geometry.canChangeJ[Region[r]->jLocations[j]] =
						Geometry.canChangeJ[Region[r]->jLocations[j]] || A(CanChange);
					Geometry.canChangeK[Region[r]->kLocations[k]] =
						Geometry.canChangeK[Region[r]->kLocations[k]] || A(CanChange);

					if (Phase[Region[r]->phaseStart]->phaseType == LIQUID)
					{
						if (A(State) != LIQUID)		//a new liquid node
							Sim.numberLiquid++;

						A(State) = LIQUID;						
						A(PhaseCodeLiquid) = Region[r]->phaseStart;
						A(PhaseCodeSolid) = PHASE_NULL;
						A(FractionSolid) = FULL_LIQUID;
					}
					
					else //solid nodes
					{
						if (A(State) == LIQUID)
							Sim.numberLiquid--;

						A(State) = SOLID;
						A(PhaseCodeSolid) = Region[r]->phaseStart;
						A(PhaseCodeLiquid) = PHASE_NULL;
						A(FractionSolid) = FULL_SOLID;
					}; //endif
						
					HistorySet(PhaseHistory, nodeA, A(State), A(PhaseCodeSolid));
				}; //endloop k
			}; //endloop j
		}; //endloop i
	}; //endloop r
			
	Sim.numberSlush = 0;	 // initially no slush nodes

	for (i=I_FIRST; i<I_LAST; i++)			//determine boundaries of interactions
		for (j=J_FIRST; j<J_LAST; j++)
			for (k=K_FIRST; k<K_LAST; k++)
			{
			
				nodeA = IJKToIndex(i, j, k);
						
				if ((i != I_EAST) && (A(CanChange) || CanChange[I_EAST][j][k]))	
					A(canInteract) |= EAST;
				if ((i != I_WEST) && (A(CanChange) || CanChange[I_WEST][j][k]))
 					A(canInteract) |= WEST;
				if ((j != J_SOUTH) && (A(CanChange) || CanChange[i][J_SOUTH][k]))
 					A(canInteract) |= SOUTH;
				if ((j != J_NORTH) && (A(CanChange) || CanChange[i][J_NORTH][k]))
 					A(canInteract) |= NORTH;
				if ((k != K_BACK) && (A(CanChange) || CanChange[i][j][K_BACK]))
 					A(canInteract) |= BACK;
				if ((k != K_FRONT) && (A(CanChange) || CanChange[i][j][K_FRONT]))
 					A(canInteract) |= FRONT;
			};
	
	return;

}; //endfunc


//______________________________________________________
//
//	PhaseNew
//      Allocates a new region or connects to existing
//______________________________________________________
PHASE *PhaseNew(int i)
{
	if (Phase[i] == NULL)
	{
		Phase[i] = new PHASE;
		
		if (Phase[i] == NULL)
			ErrorMsg("Out of memory");

		memset(Phase[i], 0, sizeof(PHASE));
		
	}; //endif

	return(Phase[i]);
}; //endfunc


//______________________________________________________
//
//	PhaseWillSolidify
//		Checks all available phases to determine if 
//		solidification is possible
//______________________________________________________
inline bool PhaseWillSolidify(const int nodeA, const double tInterface, int &phaseCode)
{
	int tThisCrit;
	int tMinCrit=MAX_DEGREES;
	int i;

	for (i=0; i< Phase[A(PhaseCodeLiquid)]->numTransitions; i++)  // Loop thru potential solid phases
	{
		tThisCrit = TSolidify(nodeA, i);		//critical temp for growth of this phase

		if ((tInterface < tThisCrit) && (tThisCrit < tMinCrit)) 	//if this is a better choice
		{
			tMinCrit = tThisCrit;
			phaseCode = *(Phase[A(PhaseCodeLiquid)]->transitionTo+i);
		}; //endif
	}; //endloop

	return (tMinCrit < MAX_DEGREES);	//true if found a phase
}; //endfunc


//_________________________________________________
// Fraction100
//		Calculates fraction in 1-dimensional case
//_________________________________________________
double Fraction100(const double posRaw, const int comp, bool &isComplete)
{
	double posn;

	posn = (posRaw < 0.0) ? 0.0 : ((posRaw > 1.0) ? 1.0 : posRaw);
	isComplete = !(posn == posRaw);		//complete melting or solidification

	switch (comp)
	{
		case 1:
			return(posn);
			break;

		case -1:
			return(1.0 - posn);
			break;
		
		default:
			ErrorMsg("Bad component call to Fraction100 at time %g", Sim.sClock.curTime);
	}; //endswitch			

	return(0);
}; //endfunc


//______________________________________________________
//
//	Fraction110
//		Fractional solid within a node for
//		(110) directions.  See RA p.729
//______________________________________________________
double Fraction110(const double posRaw, const int comp, bool &isComplete)
{
	#ifdef ISO_LINEAR
	return(Fraction100(posRaw, comp, isComplete));
	
	#else 
	
	//__________ Non-Linear Model __________
	double posn;

	posn = (posRaw < 0.0) ? 0.0 : ((posRaw > 1.0) ? 1.0 : posRaw);
	isComplete = !(posn == posRaw);		//complete melting or solidification

	switch (comp)
	{
		case (1):	//positive direction (SOUTH, EAST, DOWN)
			return( (posn <= 0.5) ? 
				2 * pow(posn,2) : 1 - 2 * pow((1 - posn),2));
			break;
	
		case (-1):	//negative direction (NORTH, WEST, UP)
			return( (posn <= 0.5) ? 
				1 - 2*pow(posn,2) : 2 * pow((1 - posn),2));
			break;

		case (0):
		case (2):
			return (0);
			break;
	}; //endswitch

	return(0);
	#endif
}; //endfunc


//______________________________________________________
//
//	Fraction111
//		Fractional solid within a node for
//		(111) directions. See RA p.734
//______________________________________________________
double Fraction111(const double posRaw, const int comp, bool &isComplete)
{
	#ifdef ISO_LINEAR
	return(Fraction100(posRaw, comp, isComplete));
	
	#else 

	//__________ Non-Linear Model __________
	double posn;

	posn = (posRaw < 0.0) ? 0.0 : ((posRaw > 1.0) ? 1.0 : posRaw);
	isComplete = !(posn == posRaw);		//complete melting or solidification

	switch (comp)
	{
		case (1):
			if (posn < 0.3333) 
				return (0.16667 * pow(3 * posn, 3));
			else if (posn < 1 - 0.3333)
				return (2 * (posn - 0.25));
			else 
				return (1 - 0.16667 * pow(3 * (1 - posn), 3));
			break;
		case (-1):
			if (posn < 0.3333)
				return (1 - 0.16667 * pow(3 * posn, 3));
			else if (posn < 1 - 0.3333)
				return (2 * (0.75 - posn));
			else
				return (0.16667 * pow(3 * (1 - posn), 3));
			break;

		case (0):
		case (2):
			return (0);
			break;
	}; //endswitch
	return(0);

	#endif
}; //endfunc


//______________________________________________________
//
//	FractionHigh
//		Fraction solid calculation for higher order nodes
//		approximate as sphere equal to node volume
//______________________________________________________
double FractionHigh(const double posRaw, const int comp, bool &isComplete)
{
	double posn;

	posn = (posRaw < 0.0) ? 0.0 : ((posRaw > 1.0) ? 1.0 : posRaw);
	isComplete = !(posn == posRaw);		//complete melting or solidification		

	return(posn);
}; //endfunc
	


//_________________________________________________
// Position100
//		Calculates interface position for 
//		(100) class interfaces
//_________________________________________________
double Position100(const double fraction, const int component)
{
	switch (component)
	{
		case 0:
			return(0.5);
			break;

		case 1:
			return(fraction);
			break;

		case -1:
			return(1.0 - fraction);
			break;
			
		case 2:
			return(0.5);
			break;
	}; //endswitch			

	return(0);
}; //endfunc


//______________________________________________________
//
//	Position110
//		Fractional position within a node for
//		(110) directions.  See RA p.729
//______________________________________________________
double Position110(const double fraction, const int comp)
{
	#ifdef ISO_LINEAR
	return(Position100(fraction, comp));
	
	#else 
	//__________ Non-Linear Model __________
	switch (comp)
	{
		case (1):	//positive direction (SOUTH, EAST, DOWN)
			return( (fraction <= 0.5) ? 
				0.5 * sqrt(2 * fraction) : 1 - 0.5 * sqrt(2 * (1-fraction)) );
			break;
		
		case (-1):	//negative direction (NORTH, WEST, UP)
			return( (fraction <= 0.5) ? 
				1 - 0.5 * sqrt(2 * fraction) : 0.5 * sqrt(2 * (1-fraction)) );
			break;

		case (0):
		case (2):
			return (0.5);
			break;
	}; //endswitch
	
	return(0);
	
	#endif

}; //endfunc


//______________________________________________________
//
//	Position111
//		Fractional position within a node for
//		(111) directions. See RA p.730
//______________________________________________________
double Position111(const double fraction, const int comp)
{
	#ifdef ISO_LINEAR
	return(Position100(fraction, comp));
	
	#else 
	//__________ Non-Linear Model __________
	switch (comp)
	{
		case (1):
			if (fraction < 0.1667) 
				return (0.3333 * pow(6 * fraction, 0.3333));
			else if (fraction < 1 - 0.1667)
				return (0.25 + 0.5 * fraction);
			else 
				return (1 - 0.3333 * pow(6 * (1 - fraction), 0.3333));
			break;
		case (-1):
			if (fraction < 0.1667)
				return (1 - 0.3333 * pow(6 * fraction, 0.3333));
			else if (fraction < 1 - 0.1667)
				return (0.75 - 0.5 * fraction);
			else
				return (0.3333 * pow(6 * (1 - fraction), 0.3333));
			break;

		case (0):
		case (2):
			return (0.5);
			break;
	}; //endswitch

	return(0);

	#endif
}; //endfunc


//______________________________________________________
//
//	PositionHigh
//		Fractional position within a node for
//		higher order directions (see RA p.778)
//______________________________________________________
double PositionHigh(const double frac)
{
	return(frac);
}; //endfunc


//______________________________________________________
//
//	PositionInterface
//		Calculates fractional position of the slush
//		interface within a node
//______________________________________________________
void PositionInterface (const int dirToLiq, const double fracSolid, double &iFrac, double &jFrac, double &kFrac)
{
    if ((fracSolid < FULL_LIQUID) || (fracSolid > FULL_SOLID))
        ErrorMsg ("Fraction solid out of range in CalcSlushPosition");

    if ((dirToLiq < 0) || (dirToLiq > ALL_DIRECTIONS))
        ErrorMsg ("Direction out of range in CalcSlushPosition");

    switch ( Interface.crystalType[dirToLiq] )
	{
		case ((int)(0x100)):
		case ((int)(0x010)):
		case ((int)(0x001)):
			iFrac = Position100(fracSolid, Interface.iCrystal[dirToLiq]);
			jFrac = Position100(fracSolid, Interface.jCrystal[dirToLiq]);
			kFrac = Position100(fracSolid, Interface.kCrystal[dirToLiq]);
			break;

		case (0x110):
		case (0x101):
		case (0x011):
			iFrac = Position110(fracSolid, Interface.iCrystal[dirToLiq]);
			jFrac = Position110(fracSolid, Interface.jCrystal[dirToLiq]);
			kFrac = Position110(fracSolid, Interface.kCrystal[dirToLiq]);

			break;

		case (0x111):
			iFrac = Position111(fracSolid, Interface.iCrystal[dirToLiq]);
			jFrac = Position111(fracSolid, Interface.jCrystal[dirToLiq]);
			kFrac = Position111(fracSolid, Interface.kCrystal[dirToLiq]);
			break;

		default:		// higher orders
			break;	// use last known positions
	}; //endswitch

	return;
}; //endfunc		
	

//______________________________________________________
//
//	TMelt
//		Melting temperature of current solid phase
//______________________________________________________
inline double TMelt(const IMATRIX phaseSolid, const int nodeA)
{
	return ( *(Phase[A(phaseSolid)]->transitionTemp) );
}; //endinline


//______________________________________________________
//
//	TSolidify
//		Solidification temperature of a solid phase (phaseNum)
//		accesible from the current liquid phase
//______________________________________________________
inline int TSolidify(const IMATRIX phaseLiquid, const int nodeA, const int phaseNum)
{
	return ( *((Phase[A(phaseLiquid)]->transitionTemp) + phaseNum) );
}; //endinline
