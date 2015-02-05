//  _________________________________________________________
// |
// |   Phase.cpp    Phase change / interface routines
// |
// |
// |   (C) 1998-99  Columbia University, MSME
// |__________________________________________________________
#include <conio.h>
#include "3dns.h"
#include "nucleation.h"
#include "parsefunc.h"		       //for LOOK_UP
#include "report.h"
#include "thermal.h"
#include "curvature.h"

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
void CalcTInterface(DMATRIX tInterface);
void ChangeSlushDir(const int nodeSlush, const int nodeOther,  const int dirToOther, const double tInterface);
void CheckNodeMelt(const int nodeSolid, const int nodeOther, const int dirSolToLiq, const double tInterface);
void CheckNodeMeltAlien(const int nodeSolid, const int nodeOther, const int dirSolToLiq, const double tInterface);
void CheckNodeSolidify(const int nodeSolid, const int nodeOther, const int dirSoltoLiq, const double tInterface);
void CheckSurfaceMelt(int nodeA, int nodeB, int dirAtoSurf);
void MotionNode (const int i, const int j, const int k);
void MotionNeighbor(const int nodeA, const int nodeB, const int dirAtoB);
void NodeMelt(const int nodeA, const int dirToLiq, const double tInterface, const double timeLeft);
void NodeMeltComplete(const int nodeA, double &fracUnused);
void NodeSolidifyComplete(const int nodeA, double &fracUnused);
inline bool PhaseWillSolidify(const int nodeA, const double tInterface, int &phaseCode);
double Position100(const double fraction, const int component);
double Position110(const double frac, const int comp);
double Position111(const double frac, const int comp);
double PositionHigh(const double frac);

inline double TMelt(int node);
inline int TSolidify(const int nodeA, const int phaseNum);

//___________________________________
// Private variables
//___________________________________
DMATRIX timeExtra;		//extra time for new interface motion
DMATRIX timeUnused;		//time unused during solidification

IMATRIX stateNew;			//new state integer
IMATRIX slushDirNew;		//new slush direction
int numNeighbors;		//number of nearest neighbors
	

inline double TIMELEFT(double iNew, double iOld)
{	return( (iNew>1) ? (iNew-1.0)/(iNew-iOld) : ((iNew<0) ? iNew/(iNew-iOld) : 0) ); }


//______________________________________________________
//
//	ChangeSlushDir
//		Change direction of slush based on nearest
//		neighbor states
//______________________________________________________
void ChangeSlushDir(const int nodeSlush, const int nodeOther,  const int dirPositive, const double tInterface)
{
	if (OTHER(stateNew) == SOLID)
	{
		if (tInterface < TMelt(nodeSlush))			//solidifying slush node	
			SLH(slushDirNew) |= dirPositive;		//turn on this direction
	}

	else if (OTHER(stateNew) == LIQUID)
	{
		if (tInterface > TMelt(nodeSlush))		//melting slush node
			SLH(slushDirNew) |= dirPositive;		//turn on this direction
	}; //endif

	if (OTHER(timeUnused) > 0)
		SLH(timeExtra) += OTHER(timeUnused) / numNeighbors;

	return;
}; //endfunc


//______________________________________________________
//
//	CheckNodeMelt
//		Checks for simple melting in a node with
//		a LIQUID neighbor in the SAME material class
//______________________________________________________
void CheckNodeMelt(const int nodeA, const int nodeOther, const int dirSolToLiq, const double tInterface)
{
	if (!A(CanChange))			//if node can't melt
		return;

	if (tInterface > TMelt(nodeA))		 //liquid always will consume solid node of
	{
		NodeMelt(nodeA, dirSolToLiq, tInterface, OTHER(timeUnused));
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
void CheckNodeMeltAlien(const int nodeSolid, const int nodeOther, const int dirSolToLiq, const double tInterface)
{
	if (!SOL(CanChange))	//if node can't melt
		return;							
							//same material class (RA p.755)
	if (OTHER(CatalyzeMelting) && (tInterface > TMelt(nodeSolid)) )
	{
		NodeMelt(nodeSolid, dirSolToLiq, tInterface, 0.0);
		OTHER(Velocity) = 0.0;
	};

	return;
}; //endfunc


//______________________________________________________
//
//	CheckNodeSolidify
//		Checks for simple solidification in a liquid node
//		with a neighbor in the SAME material class
//______________________________________________________
void CheckNodeSolidify(const int nodeA, const int nodeOther, const int dirSoltoLiq, const double tInterface)
{
	int i;
	int solidIndex = -1;			//index into the table of available solid phases

	double tCritMin = DOUBLE_INFINITY;
	double tCrit;		
	
	if (OTHER(CatalyzeFreezing))	//liquid so try find best solid phase to form (RA p.754)
	{	
		for (i = 0; i < Phase[A(PhaseCodeLiquid)]->numTransitions; i++)
		{
			tCrit = Phase[A(PhaseCodeLiquid)]->transitionTemp[i];

			if ((tInterface < tCrit) && (tCrit < tCritMin))	//can solidify to this phase
			{
				tCritMin = tCrit;
				solidIndex = i;
			}; //endif
		}; //endloop i

	} //endif

	else //could still solidify if same solid phase is available and temperature is low enough
	{
		for (i = 0; i < Phase[A(PhaseCodeLiquid)]->numTransitions; i++)
		{
			if (Phase[A(PhaseCodeLiquid)]->transitionTo[i] == OTHER(PhaseCodeSolid) )
			{
				if (tInterface < Phase[A(PhaseCodeLiquid)]->transitionTemp[i])
				{
					solidIndex = i;
					break;
				}; //endif
			}; //endif
		}; //endloop i
	
	}; //endif

	if (solidIndex >= 0)	//found a new phase to form
	{
		NodeSolidify(nodeA, dirSoltoLiq, tInterface, solidIndex, OTHER(timeUnused), H_FREEZE_INTERFACE);
		OTHER(Velocity) = 0.0;		//reset velocity in other node
	};

	return;
};

		
//______________________________________________________
//
//  CheckSurfaceMelt
//      Moves slush in a node, returns the new fraction
//      solidified
//______________________________________________________
void CheckSurfaceMelt(int nodeA, int nodeB, int dirAtoSurf)
{
	double tSurf;	//extrapolated surface temperature
	
	if ((A(State) == SOLID) && A(CanChange))
	{
		tSurf = 1.5 * A(T) - 0.5 * B(T);	//assumes equal size nodes ! change this
	
		if (tSurf > TMelt(nodeA))
			NodeMelt(nodeA, dirAtoSurf, tSurf, 0.0);
	}; //endif

	return;
}; //endfunc
		

//_________________________________________________
// HistoryClear
//		Clears flags in the history integer 
//_________________________________________________
void HistoryClear(const int nodeA, const int clrVal)
{
	A(PhaseHistory) &= H_COMPLEMENT(clrVal);
    return;
};


//_________________________________________________
// HistorySet
//		Sets flags in the history integer
//_________________________________________________
void HistorySet(const int nodeA, const int setVal)
{
	A(PhaseHistory) |= setVal;
    return;
}; //endfunc


//__________________________________________________________
// InterfaceAutomata
//		Calculate interface behavior using
//		the cellular automata model
//__________________________________________________________
void InterfaceAutomata()
{
	int i, j, k;
	int iNext, jNext, kNext;

	int nodeA;
		
	MatrixCopy(MatrixNew(&stateNew), State, 0);				//new state based on old
	MatrixCopy(MatrixNew(&slushDirNew), SlushDirection, 0);	//recalculate all slush directions anew
	MatrixZero(MatrixNew(&timeExtra));		//extra (unused) time for fully solidified/melted nodes

	//___________ STANDARD NEIGHBOR NODE INTERACTIONS ____________
	for (k=K_FIRST; kNext = (k + 1) % K_LAST, k < K_LAST; k++)	
	{
		if ( !( *(Geometry.canChangeK + k) ||  *(Geometry.canChangeK + kNext) ))
			continue;		//these 2 pages cant change

		for (j=J_FIRST; jNext = (j + 1) % J_LAST, j < J_LAST; j++)
		{
			if ( !( *(Geometry.canChangeJ + j) ||  *(Geometry.canChangeJ + jNext) ))
				continue;		//these 2 rows cant change
			
			for (i=I_FIRST; iNext = (i + 1) % I_LAST, i < I_LAST; i++)
			{	
				if ( !( *(Geometry.canChangeI + i) ||  *(Geometry.canChangeI + iNext) ))
					continue;		//these 2 columns cant change

				nodeA = IJKToIndex(i, j, k);

				if ((i < I_LAST-1) || (I_PERIODIC && Geometry.isThick[DIM_I]))
					MotionNeighbor(nodeA, IJKToIndex(iNext, j, k), EAST);
					
				if ((j < J_LAST-1) || (J_PERIODIC && Geometry.isThick[DIM_J]))
					MotionNeighbor(nodeA, IJKToIndex(i, jNext, k), SOUTH);

				if ((k < K_LAST-1) || (K_PERIODIC && Geometry.isThick[DIM_K]))
					MotionNeighbor(nodeA, IJKToIndex(i, j, kNext), BACK);
			}; //endloop i
		}; //endloop j
	}; //endloop k

	//___________ SPECIAL EDGE CONDITIONS ____________
	if (Geometry.canChangeJ[J_FIRST] && Geometry.isThick[DIM_J])	
	{											
		for (i = I_FIRST; i < I_LAST; i++)			//top surface melting
			for (k = K_FIRST; k < K_LAST; k++)
				CheckSurfaceMelt(IJKToIndex(i, J_FIRST, k), 
						IJKToIndex(i, J_FIRST+1, k), NORTH);
	}; //endif

	//___________ RECALCULATE SLUSH TEMPERATURES & POSITIONS ____________
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
					
				if (A(stateNew) == SLUSH)
				{
					if (A(slushDirNew) != A(SlushDirection))
						PositionInterface(A(slushDirNew), A(FractionSolid), 
							A(IPos), A(JPos), A(KPos) );	//update positions

					A(TInterface) = InterpolateT (i, j, k, A(IPos), A(JPos), A(KPos));
				}; //endif

			}; //endloop k
		}; //endloop j
	}; //endloop i

	MatrixZero(timeUnused);				//zero unused times from previous clock
	MatrixZero(QInterface);					//energy due to interface motion
	MatrixRecycle(State);				//dump old state
	MatrixRecycle(SlushDirection);		//dump old slush direction

	State = stateNew;				//use newly resolved states
	SlushDirection = slushDirNew;	//use new slush directions
	
	//___________ CURVATUE UNIT ____________
	for (i = I_FIRST; i < I_LAST; i++)
	{
		if (!Geometry.canChangeI)
			continue;

		for (j = J_FIRST; j < J_LAST; j++)
		{
			if (!Geometry.canChangeJ)
				continue;

			for (k = K_FIRST; k < K_LAST; k++)
			{
				if (State[i][j][k] == SLUSH)
				{
					Curvature(i, j, k, S_Curv[i][j][k], DT_GT [i][j][k]);
				};

			}; //endloop k
		}; //endloop j
	}; //endloop i
	//___________ INTRA-NODE MOTION ____________
	for (i = I_FIRST; i < I_LAST; i++)
	{
		if (!Geometry.canChangeI)
			continue;

		for (j = J_FIRST; j < J_LAST; j++)
		{
			if (!Geometry.canChangeJ)
				continue;

			for (k = K_FIRST; k < K_LAST; k++)
			{
				if (State[i][j][k] == SLUSH)
					MotionNode(i, j, k);
			}; //endloop k
		}; //endloop j
	}; //endloop i
	
	MatrixRecycle(timeExtra);			//dump residual times from this clock
										//but keep timeUnused	
	return;
}; //endfunc


//__________________________________________________________
// MotionNeighbor
//		Calculate behavior at boundary of neighboring nodes
//		initiate solidification or melting if necessary
//__________________________________________________________
void MotionNeighbor(const int nodeA, const int nodeB, const int dirAtoB)
{	
	double tInterface;
	int nodeLiquid, nodeSolid, nodeSlush;	//node indices

	#define MATCH_STATE(x) ((A(State) == x) ? nodeA : ((B(State) == x) ? nodeB : -1 ))
	#define OTHER_NODE(x)  ((x == nodeA) ? nodeB : nodeA)
	#define DIR(x, y)  ((x == nodeA) ? dirAtoB : Geometry.invertDir[dirAtoB])

	//_______________ SAME MATERIAL CLASS _____________________
	if ( A(MaterialClass) == B(MaterialClass) )
	{
		switch (A(State) | B(State))		//see interaction table RA p.750
		{
			case (SOLID | SOLID):
				if (A(PhaseCodeSolid) == B(PhaseCodeSolid))
					return;						//exactly the same phase

				tInterface = 0.5 * ( A(T) + B(T) );	//interface temperature

				CheckNodeMelt(nodeA, nodeB, DIR(nodeA, nodeB), tInterface);
				CheckNodeMelt(nodeB, nodeA, DIR(nodeB, nodeA), tInterface);
	
				break;

			case (LIQUID | LIQUID):
				return;
				
				break;

			case (SOLID | LIQUID):
				nodeSolid = MATCH_STATE(SOLID);				
				nodeLiquid = OTHER_NODE(nodeSolid);

				tInterface = 0.5 * ( A(T) + B(T) );	//interface temperature

				CheckNodeMelt(nodeSolid, nodeLiquid, DIR(nodeSolid, nodeLiquid), tInterface );
				CheckNodeSolidify(nodeLiquid, nodeSolid, DIR(nodeSolid, nodeLiquid), tInterface);

				break;

			case (SOLID | SLUSH):
				nodeSolid = MATCH_STATE(SOLID);
				nodeSlush = OTHER_NODE(nodeSolid);

				ChangeSlushDir(nodeSlush, nodeSolid, DIR(nodeSolid, nodeSlush), SLH(TInterface));
				break;

			case (LIQUID | SLUSH):				
				nodeLiquid = MATCH_STATE(SOLID);
				nodeSlush = OTHER_NODE(nodeLiquid);

				ChangeSlushDir(nodeSlush, nodeLiquid, DIR(nodeSlush, nodeLiquid), SLH(TInterface));
				break;

			case (SLUSH | SLUSH):
				break;

		}; //endswitch
	}

	//_______________ ALIEN MATERIAL CLASS _____________________
	else
	{
		switch (A(State) | B(State))		//see interaction table RA p.750
		{
			case (SOLID | SOLID):
				tInterface = 0.5 * ( A(T) + B(T) );	//interface temperature

				CheckNodeMeltAlien(nodeA, nodeB, DIR(nodeA, nodeB), tInterface);
				CheckNodeMeltAlien(nodeB, nodeA, DIR(nodeB, nodeA), tInterface);

				break;

			case (LIQUID | LIQUID):
				return;
				
				break;

			case (SOLID | LIQUID):
				nodeSolid = MATCH_STATE(SOLID);				
				nodeLiquid = OTHER_NODE(nodeSolid);
				
				tInterface = 0.5 * ( A(T) + B(T) );	//interface temperature

				NucleateHeterogeneous(nodeLiquid, nodeSolid, DIR(nodeSolid, nodeLiquid), tInterface);
				CheckNodeMeltAlien(nodeSolid, nodeLiquid, DIR(nodeSolid, nodeLiquid), tInterface);
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


//______________________________________________________
//
//  MotionNode
//      Moves slush in a node, returns the new fraction
//      solidified.  Assumes that SlushDirection, State,
//		and IJK Positions are already resolved for this
//		clock.
//______________________________________________________
void MotionNode (const int i, const int j, const int k)
{
	int nodeA = IJKToIndex(i, j, k);	//index into arrays
	bool isComplete;			//complete melting or solidification

	double clockUnused = 0.0;	//left over clock
	double timeXtra = A(timeExtra) * numNeighbors;
	double fractionOld;		//old fraction solid and position

/*	if ((Sim.curTime > 0E-9) && 
	(((abs(crystalI[A(SlushDirection)]) == 1) && (abs(crystalJ[A(SlushDirection)]) == 1)) ||
	((abs(crystalJ[A(SlushDirection)]) == 1) && (abs(crystalK[A(SlushDirection)]) == 1)) ||
	((abs(crystalK[A(SlushDirection)]) == 1) && (abs(crystalI[A(SlushDirection)]) == 1)))) 
	{
		printf("i=%3d j=%2d k=%2d   %2d %2d %2d   DT_GT=%6.3f K  S_Curv=%6.6f nm2\n", i, j, k, 
		crystalI[A(SlushDirection)], crystalJ[A(SlushDirection)], crystalK[A(SlushDirection)],
		A(DT_GT), A(S_Curv)*1E14); 
		getche();
	}; 
*/
	A(Velocity) = 1.2013E7 * exp( -1.65E-19/(1.381E-23 * A(TInterface))) * (1-exp(-4206 * ((1685 - 
		A(DT_GT) ) - A(TInterface) )/(5E22 * (1685 - A(DT_GT) ) * 1.381E-23 * A(TInterface) )));

	if (((abs(crystalI[A(SlushDirection)]) == 1) && (abs(crystalJ[A(SlushDirection)]) == 1) && (crystalK[A(SlushDirection)] == 0)) ||
		((abs(crystalJ[A(SlushDirection)]) == 1) && (abs(crystalK[A(SlushDirection)]) == 1) && (crystalI[A(SlushDirection)] == 0)) ||
		((abs(crystalK[A(SlushDirection)]) == 1) && (abs(crystalI[A(SlushDirection)]) == 1) && (crystalJ[A(SlushDirection)] == 0))) 
		A(Velocity) = 1.05 * A(Velocity); //1.07

	fractionOld = A(FractionSolid);
	A(FractionSolid) += (Sim.curClock + timeXtra) * A(Velocity) * A(S_Curv)/A(Volume);
	clockUnused = TIMELEFT(A(FractionSolid), fractionOld);
	isComplete = !(clockUnused == 0.0);
	A(FractionSolid) = (A(FractionSolid) < 0.0) ? 0.0 : ((A(FractionSolid) > 1.0) ? 1.0 : A(FractionSolid));
	PositionInterface(A(slushDirNew), A(FractionSolid), A(IPos), A(JPos), A(KPos) );	


	if (isComplete && (A(FractionSolid) > (0.99 * FULL_SOLID)))
		NodeSolidifyComplete(nodeA, clockUnused);

	else if (isComplete && (A(FractionSolid) < (0.01 + FULL_LIQUID)))
		NodeMeltComplete(nodeA, clockUnused);

	if (fabs(A(FractionSolid) - fractionOld)  > Sim.maxErrPhase)
		ErrorMsg("Interface %2d %2d %2d moved too far (FSolNew=%g FSolOld=%g) in node [%d][%d][%d] at time %g  dT_GT=%g Surf=%g[nm2] (Vol)^2/3=%g[nm2]",
			crystalI[A(SlushDirection)], crystalJ[A(SlushDirection)], crystalK[A(SlushDirection)],
			A(FractionSolid), fractionOld, i, j, k, Sim.curTime, A(DT_GT), A(S_Curv) * 1E14, (pow(A(Volume),0.666667)) * 1e14 );

    A(QInterface) = (A(FractionSolid) - fractionOld) * A(Volume) *
    	LOOK_UP(Phase[A(PhaseCodeSolid)]->enthalpyH, A(TInterface) );		//energy associated with move

    return;
}; //endfunc

//_________________________________________________
// NodeMelt
//		Initiates melting in nodeA from neighbor nodeB.
//		If nodeA is already slush then just changes direction
//_________________________________________________
void NodeMelt(const int nodeA, const int dirToLiq, const double tInterface, const double timeLeft)
{
	A(stateNew) = SLUSH;

	A(slushDirNew) = dirToLiq;
	
	PositionInterface(A(slushDirNew), A(FractionSolid),	//update fractional positions
							A(IPos), A(JPos), A(KPos) );

	A(TInterface) = tInterface;

	HistorySet(nodeA, (H_MELT_INTERFACE | (H_MELT_N * dirToLiq)));
		
	Sim.numberSlush++;

	A(PhaseCodeLiquid) = *(Phase[A(PhaseCodeSolid)]->transitionTo);

	A(timeExtra) = timeLeft / (double)numNeighbors;

	return;
}; //endfunc


//_________________________________________________
// NodeMeltComplete
//		Complete solidification of nodeA.
//		Always occurs after (stateNew --> State) update
//_________________________________________________
void NodeMeltComplete(const int nodeA, double &fracUnused)
{
	A(timeUnused) = Sim.curClock * fracUnused;

    A(State) = LIQUID;
	A(SlushDirection) = 0;

    Sim.numberSlush--;
	Sim.numberLiquid++;

	HistorySet(nodeA, H_MELT_COMPLETE);  //flag complete melting
	HistoryClear(nodeA, (H_FREEZE_HETEROGENEOUS | H_FREEZE_HOMOGENEOUS |
		H_FREEZE_INTERFACE | H_FREEZE_N | H_FREEZE_E | H_FREEZE_S | H_FREEZE_W));

	return;
}; //endfunc


//_________________________________________________
// NodeSolidify
//		Initiates solidification in nodeA from neighbor nodeB.
//		If nodeA is already slush, just changes direction
//_________________________________________________
void NodeSolidify(IMATRIX thisState, IMATRIX thisDir, const int nodeLiquid, const int dirToLiq, 
				const double tInterface, const int solidIndex, const double timeLeft, const int hCode)
{
	LIQ(thisState) = SLUSH;
	LIQ(thisDir) = dirToLiq;
	LIQ(PhaseCodeSolid) = *(Phase[LIQ(PhaseCodeLiquid)]->transitionTo + solidIndex);
	
	PositionInterface(LIQ(thisDir), LIQ(FractionSolid),	//update fractional positions
							LIQ(IPos), LIQ(JPos), LIQ(KPos) );

	LIQ(TInterface) = tInterface;	//update tInterface

	HistorySet(nodeLiquid, (hCode | (H_FREEZE_N * dirToLiq)));
	Sim.numberSlush++;
	Sim.numberLiquid--;

	LIQ(timeExtra) += timeLeft / (double)numNeighbors;		//extra time from neighbor


	return;
}; //endfunc

void NodeSolidify(const int nodeLiquid, const int dirToLiq, 
					const double tInterface, const int solidIndex, 
					const double timeLeft, const int hCode)
{	
	NodeSolidify(stateNew, slushDirNew, nodeLiquid, dirToLiq, 
				tInterface, solidIndex, timeLeft, hCode);

	return;
}; //endfunc


//_________________________________________________
// NodeSolidifyComplete
//		Complete solidification of nodeA.
//		Always occurs after (stateNew-->State) update
//_________________________________________________
void NodeSolidifyComplete(const int nodeA, double &fracUnused)
{        
	A(timeUnused) = Sim.curClock * fracUnused;

	A(State) = SOLID;
	A(SlushDirection) = 0;

	Sim.numberSlush--;

	HistorySet(nodeA, H_FREEZE_COMPLETE);

	return;
}; //endfunc


//_________________________________________________
// PhaseCleanup
//_________________________________________________
void PhaseCleanup ()
{
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
    MatrixFree (DT_GT);     // @
    MatrixFree (S_Curv);    // @

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
	MatrixNew (&IPos);
	MatrixNew (&JPos);
	MatrixNew (&KPos);
	MatrixNew (&MaterialClass);
    MatrixNew (&PhaseCodeLiquid);
	MatrixNew (&PhaseCodeSolid);
    MatrixZero(MatrixNew(&PhaseHistory));
    MatrixNew (&QInterface);
    MatrixNew (&SlushDirection);
    MatrixNew (&TInterface);
    MatrixNew (&State);
    MatrixNew (&Velocity);
    MatrixNew (&DT_GT);    // @
    MatrixNew (&S_Curv);   // @

	MatrixZero(MatrixNew(&timeUnused));
	MatrixNew(&CatalyzeFreezing);
	MatrixNew(&CatalyzeMelting);

	Geometry.canChangeI = new bool[I_LAST];					//if columns can change
	memset(Geometry.canChangeI, 0, sizeof(bool)*I_LAST);

	Geometry.canChangeJ = new bool[J_LAST];					//if rows can change
	memset(Geometry.canChangeJ, 0, sizeof(bool)*J_LAST);

	Geometry.canChangeK = new bool[K_LAST];					//if pages can change
	memset(Geometry.canChangeK, 0, sizeof(bool)*K_LAST);

	PhaseNew(0);	//initialize NULL phase

	for (p=1; p<Sim.numberPhases; p++)	//first phase is '1', while '0' is NULL phase
	{
		if (Phase[p] == NULL)
			ErrorMsg("Phase %d not defined", p);

		if (!InRange(Phase[p]->phaseType, 0, MAX_PHASES))
			ErrorMsg("Bad phase type in phase %d", p);

		if (Phase[p]->transitionTo[0] != PHASE_NULL)	// this phase can change
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
						A(State) = LIQUID;
						A(PhaseCodeLiquid) = Region[r]->phaseStart;
						A(PhaseCodeSolid) = PHASE_NULL;
						A(FractionSolid) = FULL_LIQUID;
						Sim.numberLiquid++;
					}
					
					else //solid nodes
					{
						A(State) = SOLID;
						A(PhaseCodeSolid) = Region[r]->phaseStart;
						A(PhaseCodeLiquid) = PHASE_NULL;
						A(FractionSolid) = FULL_SOLID;
					}; //endif
						
				}; //endloop k
			}; //endloop j
		}; //endloop i
	}; //endloop r
	
	#define PSEUDO(x,y) ((i&x)>0) ? (((i&y)>0) ? 2 : 1) : (((i&y)>0) ? -1 : 0)
	
	for (i = 0; j = 0, i < MAX_DIRECTIONS; i++)		//construct inverse directions table
	{
		crystalI[i] = PSEUDO(EAST, WEST);
		crystalJ[i] = PSEUDO(SOUTH, NORTH);
		crystalK[i] = PSEUDO(BACK, FRONT);

		crystalClass[i] = abs(crystalI[i] * 256) + abs(crystalJ[i] * 16) +
							abs(crystalK[i]);

	}; //endloop i
		
	for (i = 0, numNeighbors=0; i < 3; i++)
		numNeighbors += 2 * (Geometry.isThick[i] ? 1 : 0);

	Sim.numberSlush = 0;	 // initially no slush nodes
	
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
double Position110(const double frac, const int comp)
{
	switch (comp)
	{
		case (1):	//positive direction (SOUTH, EAST, DOWN)
			return( (frac <= 0.5) ? 
				0.5 * sqrt(2 * frac) : 1 - 0.5 * sqrt(2 * (1-frac)) );
			break;
		
		case (-1):	//negative direction (NORTH, WEST, UP)
			return( (frac <= 0.5) ? 
				1 - 0.5 * sqrt(2 * frac) : 0.5 * sqrt(2 * (1-frac)) );
			break;

		case (0):
		case (2):
			return (0.5);
			break;
	}; //endswitch
	return(0);
}; //endfunc


//______________________________________________________
//
//	Position111
//		Fractional position within a node for
//		(111) directions. See RA p.730
//______________________________________________________
double Position111(const double frac, const int comp)
{
	switch (comp)
	{
		case (1):
			if (frac < 0.1667) 
				return (0.3333 * pow(6 * frac, 0.3333));
			else if (frac < 1 - 0.1667)
				return (0.25 + 0.5 * frac);
			else 
				return (1 - 0.3333 * pow(6 * (1 - frac), 0.3333));
			break;
		case (-1):
			if (frac < 0.1667)
				return (1 - 0.3333 * pow(6 * frac, 0.3333));
			else if (frac < 1 - 0.1667)
				return (0.75 - 0.5 * frac);
			else
				return (0.3333 * pow(6 * (1 - frac), 0.3333));
			break;

		case (0):
		case (2):
			return (0.5);
			break;
	}; //endswitch
	return(0);
}; //endfunc


//______________________________________________________
//
//	PositionHigh
//		Fractional position within a node for
//		higher order directions (see RA p.778)
//______________________________________________________
double PositionHigh(const double frac)
{
	return(pow(frac, 0.333));
}; //endfunc


//______________________________________________________
//
//	PositionInterface
//		Calculates fractional position of the slush
//		interface within a node
//______________________________________________________
void PositionInterface (const int dirSlush, const double fracSolid, double &iFrac, double &jFrac, double &kFrac)
{
    if ((fracSolid < FULL_LIQUID) || (fracSolid > FULL_SOLID))
        ErrorMsg ("Fraction solid out of range in CalcSlushPosition");

    if ((dirSlush < 0) || (dirSlush > MAX_DIRVALUE))
        ErrorMsg ("Direction out of range in CalcSlushPosition");
	
    switch ( crystalClass[dirSlush] )
	{
		case ((int)(0x100)):
		case ((int)(0x010)):
		case ((int)(0x001)):
			iFrac = Position100(fracSolid, crystalI[dirSlush]);
			jFrac = Position100(fracSolid, crystalJ[dirSlush]);
			kFrac = Position100(fracSolid, crystalK[dirSlush]);
			break;

		case (0x110):
		case (0x101):
		case (0x011):
			iFrac = Position110(fracSolid, crystalI[dirSlush]);
			jFrac = Position110(fracSolid, crystalJ[dirSlush]);
			kFrac = Position110(fracSolid, crystalK[dirSlush]);
			break;

		case (0x111):
			iFrac = Position111(fracSolid, crystalI[dirSlush]);
			jFrac = Position111(fracSolid, crystalJ[dirSlush]);
			kFrac = Position111(fracSolid, crystalK[dirSlush]);
			break;

		default:		// use spherical center with iFrac as radius
			iFrac = PositionHigh(fracSolid);
			jFrac = kFrac = 0.5;
	}; //endswitch

	return;
}; //endfunc		
	

//______________________________________________________
//
//	TMelt
//		Melting temperature of current solid phase
//______________________________________________________
inline double TMelt(const int nodeA)
{
	return ( *(Phase[A(PhaseCodeSolid)]->transitionTemp) );
}; //endinline


//______________________________________________________
//
//	TSolidify
//		Solidification temperature of a solid phase (phaseNum)
//		accesible from the current liquid phase
//______________________________________________________
inline int TSolidify(const int nodeA, const int phaseNum)
{
	return ( *((Phase[A(PhaseCodeLiquid)]->transitionTemp) + phaseNum) );
}; //endinline
