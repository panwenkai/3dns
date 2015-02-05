//  _________________________________________________________
// |
// |   Thermal.cpp   Heat flow related routines
// |
// |
// |   (C) 1998-99  Columbia University, MSME
// |   (C) 2000-01	Harvard University, DEAS
// |__________________________________________________________

#include "3dns.h"
#include "parsefunc.h"
#include "phase.h"
#include "readdata.h"
#include "energy.h"		       //for QLaser

#define EXT_LEVEL
#include "thermal.h"
#undef EXT_LEVEL

#define FORWARD ((int)0)
#define REVERSE ((int)1)

//_________________________________________
// Private functions
//_________________________________________
void CalcCp (DMATRIX cp, DMATRIX qT);
void CalcKEffective (DMATRIX kEff, DMATRIX kTherm, int direction);
void CalcKNode (DMATRIX kTherm, DMATRIX qT);
void CalcDelta (DMATRIX tDelta, DMATRIX tCurrent, int direction);
void CalcTNew(DMATRIX dTnew,
				DMATRIX delToldI, DMATRIX delToldJ, DMATRIX delToldK,
				DMATRIX kOldI, DMATRIX kOldJ, DMATRIX kOldK, 
				DMATRIX kNewI, DMATRIX kNewJ, DMATRIX kNewK, const int dir);
void CheckDT(DMATRIX dT, DMATRIX dtForward, DMATRIX dtReverse);
inline double KNeighbor(int nodeA, int nodeB, DMATRIX kTherm, DMATRIX del1, DMATRIX del2);
//void MatrixShift(DMATRIX newMatx, DMATRIX oldMatx, int direction);

//_________________________________________
// Private variables
//_________________________________________
DMATRIX Cp;			//heat capacity per node [J/K]
DMATRIX dTNew;		//new temperature changes
DMATRIX K;			//thermal conductivity of each node


//___________________________________________________
// CalcCp
//        Calculate heat capacity at each node [K / Joule]
//        LOOK_UP returns inverse heat capacity [(K cm3 / Joule)]
//___________________________________________________
void CalcCp (DMATRIX cp, DMATRIX qT)
{
	int i, j, k;
	int nodeA;

	for (i = I_FIRST; i < I_LAST; i++)
    {
		for (j = J_FIRST; j < J_LAST; j++)
		{
			for (k = K_FIRST; k < K_LAST; k++)
			{
				nodeA = IJKToIndex(i, j, k);
				
				if (A(State) == SOLID)
					A(cp) = LOOK_UP (Phase[A(PhaseCodeSolid)]->heatCapacity, A(qT)) * A(Volume);
		
				else if (A(State) == LIQUID)
					A(cp) = LOOK_UP (Phase[A(PhaseCodeLiquid)]->heatCapacity, A(qT)) * A(Volume);
	  
				else if (A(State) == SLUSH)
					A(cp) = ((LOOK_UP (Phase[A(PhaseCodeSolid)]->heatCapacity, A(qT)) * A(FractionSolid)) +
							 (LOOK_UP (Phase[A(PhaseCodeLiquid)]->heatCapacity, A(qT)) * (1-A(FractionSolid)) ) ) *
							 A(Volume);
				else 
					ErrorMsg("Unrecognized state value at t=%9.5g", Sim.sClock.curTime);
			} //endloop k
		} //endloop j
	} //endloop i
} //endfunc


//___________________________________________________
// CalcKEffective
//		Calculate effective thermal resistance between
//		neighbors
//___________________________________________________
void CalcKEffective (DMATRIX kEff, DMATRIX kTherm, int direction)
{
	int i, j, k;
	int iPlus, jPlus, kPlus; 
	int nodeA, nodeB;

  	switch (direction)		       //outside of i/j loops for speed
    {
		case EAST:
      		for (i = I_FIRST; i < I_LAST; i++)
			{
				if ((i == I_LAST-1) && ! I_PERIODIC)	//adiabatic edge
				{
					for (j = J_FIRST; j < J_LAST; j++)
						for (k = K_FIRST; k < K_LAST; k++)
							kEff[i][j][k] = 0.0;
					continue;
				}

				else
					iPlus = ((i == I_LAST-1) ? I_FIRST : i+1);

				for (j = J_FIRST; j < J_LAST; j++)
				{
					for (k = K_FIRST; k < K_LAST; k++)
					{
						nodeA = IJKToIndex(i, j, k);
						nodeB = IJKToIndex(iPlus, j, k);

	  					A(kEff) = KNeighbor(nodeA, nodeB, kTherm, DelX, AreaYZ);
					} //endloop k
				} //endloop j
			} //endloop i 

			break;

		case SOUTH:
      		for (j = J_FIRST; j < J_LAST; j++)
			{
				if ((j == J_LAST-1) && ! J_PERIODIC)	//adiabatic edge
				{
					for (i = I_FIRST; i < I_LAST; i++)
						for (k = K_FIRST; k < K_LAST; k++)
							kEff[i][j][k] = 0.0;
					continue;
				}

				else
					jPlus = ((j == J_LAST-1) ? J_FIRST : j+1);

				for (i = I_FIRST; i < I_LAST; i++)
				{
					for (k = K_FIRST; k < K_LAST; k++)
					{
						nodeA = IJKToIndex(i, j, k);
						nodeB = IJKToIndex(i, jPlus, k);

	  					A(kEff) = KNeighbor(nodeA, nodeB, kTherm, DelY, AreaXZ);
					} //endloop k
				} //endloop i
			} //endloop j 

			break;

		case BACK:
      		for (k = K_FIRST; k < K_LAST; k++)
			{
				if ((k == K_LAST-1) && ! K_PERIODIC)	//adiabatic edge
				{
					for (i = I_FIRST; i < I_LAST; i++)
						for (j = J_FIRST; j < J_LAST; j++)
							kEff[i][j][k] = 0.0;
					continue;
				}

				else
					kPlus = ((k == K_LAST-1) ? K_FIRST : k+1);

				for (i = I_FIRST; i < I_LAST; i++)
				{
					for (j = J_FIRST; j < J_LAST; j++)
					{
						nodeA = IJKToIndex(i, j, k);
						nodeB = IJKToIndex(i, j, kPlus);

	  					A(kEff) = KNeighbor(nodeA, nodeB, kTherm, DelZ, AreaXY);
					} //endloop j
				} //endloop i
			} //endloop k 

			break;

	} //endswitch
} //endfunc


//___________________________________________________
// CalcKNode
//        Calculate thermal conductivity of each node
//___________________________________________________
void CalcKNode (DMATRIX kTherm, DMATRIX qT)
{
	int i, j, k;
	int nodeA;

	for (i = I_FIRST; i < I_LAST; i++)
	{
		for (j = J_FIRST; j < J_LAST; j++)
		{
			for (k = K_FIRST; k < K_LAST; k++)
			{
				nodeA = IJKToIndex(i, j, k);

				if (A(State) == SOLID)
					A(kTherm) = LOOK_UP (Phase[A(PhaseCodeSolid)]->thermConductivity, A(qT));

				else if (A(State) == LIQUID)
					A(kTherm) = LOOK_UP (Phase[A(PhaseCodeLiquid)]->thermConductivity, A(qT));

				else			       //assume SLUSH
				{
					double kliq, ksol;
					ksol = LOOK_UP (Phase[A(PhaseCodeSolid)]->thermConductivity, A(qT));
					kliq = LOOK_UP (Phase[A(PhaseCodeLiquid)]->thermConductivity, A(qT));

					A(kTherm) = (kliq * ksol) /
							(A(FractionSolid) * (kliq - ksol) + ksol);
				} //endelse
			} //endloop k
		} //endloop j
    }  //endloop i
} //endfunc


//___________________________________________________
// CalcTDelta
//		Calculate delta temperatures
//___________________________________________________
void CalcDelta (DMATRIX mDelta, DMATRIX mSource, int direction)
{
	int i, j, k;
	int iNeighbor, jNeighbor, kNeighbor;

	unsigned int nodeA;
	
	switch (direction)
	{
		case (EAST):
			for (i = 0; iNeighbor = I_EAST, i < I_LAST; i++)
				for (j = 0; j < J_LAST; j++)
					for (k = 0; k < K_LAST; k++)
					{
						nodeA = IJKToIndex(i, j, k);
						A(mDelta) = mSource[iNeighbor][j][k] - A(mSource);
					}; //endloop k
			break;

		case (SOUTH):
			for (i = 0; i < I_LAST; i++)
				for (j = 0; jNeighbor = J_SOUTH, j < J_LAST; j++)
					for (k = 0; k < K_LAST; k++)
					{
						nodeA = IJKToIndex(i, j, k);
						A(mDelta) = mSource[i][jNeighbor][k] - A(mSource);
					}; //endloop k
			break;

		case (BACK):
			for (i = 0; i < I_LAST; i++)
				for (j = 0; j < J_LAST; j++)
					for (k = 0; kNeighbor = K_BACK, k < K_LAST; k++)
					{
						nodeA = IJKToIndex(i, j, k);
						A(mDelta) = mSource[i][j][kNeighbor] - A(mSource);
					}; //endloop k
			break;

		default:
			ErrorMsg("Bad call to CalcDelta-- illegal direction");
	} //endswitch

	return;
} //endfunc

	
//___________________________________________________
// CalcTNew
//        Calculate new temperatures based on direction
//        returns maximum error
//___________________________________________________
void CalcTNew(DMATRIX dTnew,
				DMATRIX delToldI, DMATRIX delToldJ, DMATRIX delToldK,
				DMATRIX kOldI, DMATRIX kOldJ, DMATRIX kOldK, 
				DMATRIX kNewI, DMATRIX kNewJ, DMATRIX kNewK, const int dir)
{
	int i, j, k;
	int iNew, jNew, kNew;		//coordinates in the new direction
	int iOld, jOld, kOld;		//coordinates in the old direction
	int iPrev, jPrev, kPrev;	//previous node locations
	int nodeA, nodeNew, nodeOld;
	double val1, val2;

	int sign = ((dir == FORWARD) ? 1 : -1);

	MatrixZero(dTnew);	//initially no change from Told
	
								//see RA p.746 for logic tables
	for (   i = ((dir == FORWARD) ? I_FIRST : I_LAST-1);
			i != ((dir == FORWARD) ? I_LAST : I_FIRST-1); 
			i += sign )
	{
		iOld = (dir == FORWARD) ? i : ((i + I_LAST - 1) % I_LAST);
		iNew = (dir == FORWARD) ? ((i + I_LAST - 1) % I_LAST) : i;
		iPrev = (dir == FORWARD) ? iNew : ((i + 1) % I_LAST);
		
		for (   j = ((dir == FORWARD) ? J_FIRST : J_LAST-1); 
				j != ((dir == FORWARD) ? J_LAST : J_FIRST-1); 
				j += sign )
		{
			jOld = (dir == FORWARD) ? j : ((j + J_LAST - 1) % J_LAST);
			jNew = (dir == FORWARD) ? ((j + J_LAST - 1) % J_LAST) : j;
			jPrev = (dir == FORWARD) ? jNew : ((j + 1) % J_LAST);

			for (   k = ((dir == FORWARD) ? K_FIRST : K_LAST-1); 
					k != ((dir == FORWARD) ? K_LAST : K_FIRST-1); 
					k += sign )
			{				
				nodeA = IJKToIndex(i, j, k);
				val1 = val2 = 0;

				if ( Geometry.isThick[DIM_I] )  //if has I dimension
				{
					nodeOld = IJKToIndex(iOld, j, k);
					nodeNew = IJKToIndex(iNew, j, k);

					val1 += NEW(kNewI) * dTnew[iPrev][j][k] + (OLD(kOldI) * OLD(delToldI) - 
							NEW(kNewI) * NEW(delToldI)) * sign;	
					val2 += NEW(kNewI);	//denominator
				}; //endif

				if ( Geometry.isThick[DIM_J] )  //if has J dimension
				{
					nodeOld = IJKToIndex(i, jOld, k);
					nodeNew = IJKToIndex(i, jNew, k);

					val1 += NEW(kNewJ) * dTnew[i][jPrev][k] + (OLD(kOldJ) * OLD(delToldJ) - 
							NEW(kNewJ) * NEW(delToldJ)) * sign;	
					val2 += NEW(kNewJ);	//denominator
				}; 

				if ( Geometry.isThick[DIM_K] )
				{
					kOld = (dir == FORWARD) ? k : ((k + K_LAST - 1) % K_LAST);
					kNew = (dir == FORWARD) ? ((k + K_LAST - 1) % K_LAST) : k;
					kPrev = (dir == FORWARD) ? kNew : ((k + 1) % K_LAST);

					nodeOld = IJKToIndex(i, j, kOld);
					nodeNew = IJKToIndex(i, j, kNew);

					val1 += NEW(kNewK) * dTnew[i][j][kPrev] + (OLD(kOldK) * OLD(delToldK) - 
							NEW(kNewK) * NEW(delToldK)) * sign;	
					val2 += NEW(kNewK);	//denominator
				}; 

				A(dTnew) = (A(QInterface) + A(QLaser) +  2.0 * Sim.sClock.curDTime * val1) / 
							(A(Cp) + 2.0 * Sim.sClock.curDTime * val2);
			
			}; //endloop k
		}; //endloop j
	}; //endloop k

  	return;		    //maximum error between Tnew and Tcompare
} //endfunc


//___________________________________________________
// CheckDT
//    Calculate mean temperature changes and check
//		range
//___________________________________________________
void CheckDT(DMATRIX dT, DMATRIX dtForward, DMATRIX dtReverse)
{
	int nodeA;

	Sim.dTmax = 0;		//maximum node temperature change

  	for (nodeA = 0; nodeA < I_LAST*J_LAST*K_LAST; nodeA++)
	{
		A(dT) = 0.5 * (A(dtForward) + A(dtReverse));

		if (fabs(A(dT)) > Sim.dTmax)
		{
			Sim.dTmax = fabs(A(dT));		//maximum temperature change
			Sim.dTnode = nodeA;
		} //endif

	} //endloop nodeA

	return;
} //endfunc


//___________________________________________________
// FractionCalc
//        
//		Calculate fractional position used in temperature
//		interpolation.  This allows better estimates
//		of temperatures across nodes of dissimilar
//		materials.  See LAB notes p.1067 and Dev_kVacuum.nb
//		(JPL 11-13-1)
//___________________________________________________
double FractionCalc (double frac, DMATRIX del, const int nodeA, const int nodeB)
{
	double d;

	d = (B(del) * A(K) + A(del) * B(K));

	if (d > 0)		// everything but vacuum
	{
		if (frac > 0.5)    // in node A
			return ( ((2 * frac - 1.0) * A(del) * B(K)) / d );

		else				// in node B
			return ( (2 * frac * B(del) * A(K) + A(del) * B(K)) / d );
	}

	else		// both nodes A & B are vacuum
	{
		if (frac <= 0.5)
			return (frac + 0.5);
		else
			return (frac - 0.5);
	}	//endif
} //endfunc			



//___________________________________________________
// HeatFlow
//        Calculate heat flow, update temperatures
//___________________________________________________
int HeatFlow ()
{
	DMATRIX delTEast, delTSouth, delTBack;
	DMATRIX kEast, kSouth, kBack;
	DMATRIX dTf, dTr;		//new temperatures calculated in forward & reverse
	
	Sim.calcHeatFlow = true;
	Sim.dTmax = 0;

  	CalcCp (MatrixNew(&Cp), T);	       //heat capacity in each node [J/K]
  	CalcKNode (K, T);			   //thermal conductivity of each node

	if (Geometry.isThick[DIM_I])
	{	
  		CalcKEffective(MatrixNew(&kEast), K, EAST);		//effective thermal resistances
		CalcDelta(MatrixNew(&delTEast), T, EAST);		//temperature gradients
	} //endif

	if (Geometry.isThick[DIM_J])
	{
  		CalcKEffective(MatrixNew(&kSouth), K, SOUTH);	//effective thermal resistances
		CalcDelta(MatrixNew(&delTSouth), T, SOUTH);	//temperature gradients
	} //endif

	if (Geometry.isThick[DIM_K])
	{
		CalcKEffective(MatrixNew (&kBack), K, BACK);
		CalcDelta(MatrixNew(&delTBack), T, BACK);
	} //endif
  	CalcTNew (MatrixNew(&dTf),	 //find new forward temp
			  delTEast, delTSouth, delTBack, 
			  kEast, kSouth, kBack, kEast, kSouth, kBack, FORWARD);  	

  	CalcTNew (MatrixNew (&dTr),	//find new forward temp
			  delTEast, delTSouth, delTBack, 
			  kEast, kSouth, kBack, kEast, kSouth, kBack, REVERSE);  
			  
	MatrixRecycle (delTEast);	//dump unneccesary matrices
	MatrixRecycle (delTSouth);
	MatrixRecycle (delTBack);
	MatrixRecycle (Cp);		       
  	MatrixRecycle (kEast);
	MatrixRecycle (kSouth);
	MatrixRecycle (kBack);

	CheckDT(MatrixNew(&dTNew), dTf, dTr);	   //calculate mean temperature change
	
  	MatrixRecycle (dTf);
  	MatrixRecycle (dTr);

	if (Sim.dTmax > Sim.dTemperatureHi[(THIS_ERA)])
		throw(CHANGE_HI);		

	return(CHANGE_OK);		//stability result			
} //endfunc


//___________________________________________________
// HeatFlowDump
//    Dump heat flow changes
//___________________________________________________
void HeatFlowDump()
{
	if (Sim.calcHeatFlow)
	{
		MatrixRecycle(dTNew);	//dump this calculation
		Sim.calcHeatFlow = false;
	}; //endif

	return;
} //endfunc


//___________________________________________________
// HeatFlowFinalize
//    Finalize heat flow changes
//___________________________________________________
void HeatFlowFinalize()
{
	const double TOO_HOT = 0.95 * MAX_DEGREES;   //maximum allowable temperature
	const double TOO_COLD = 1.05 * MIN_DEGREES; 	//just above absolute zero
	//const double tooSteep = MAX_GRADIENT;		//maximum thermal gradient

	int i, j, k;
	int nodeA;

	for (nodeA=0; nodeA < I_LAST * J_LAST * K_LAST; nodeA++)
	{
		A(T) += A(dTNew);		//add temperature increment

		if (A(T) > TOO_HOT)
		{
			IndexToIJK(nodeA, i, j, k);		//get ijk coordinates

			if (j == Laser.beginAbsorption[i][k])
				ErrorMsg ("Node [%d][%d][%d] is too hot at t=%9.5g, check laser energy",
							i, j, k, Sim.sClock.curTime);
			else
                ErrorMsg ("Node [%d][%d][%d] is too hot at t=%9.5g, likely an instability",
                   			i, j, k, Sim.sClock.curTime);
		} //endif

		if (A(T) < TOO_COLD)
		{
			IndexToIJK(nodeA, i, j, k);		//get ijk coordinates
            		
			ErrorMsg("Node [%d][%d][%d] is too cold at t=%9.5g, likely an instability",
                		i, j, k, Sim.sClock.curTime);
		} //endif

	} //endloop

	MatrixRecycle(dTNew);
	return;
} //endfunc


//______________________________________________________
//
//      Interpolate3D
//              Calculates interpolated temperatures anywhere in a node
//______________________________________________________
double Interpolate3D( double tESD, double tESU, double tEND, double tENU,
					  double tWSD, double tWSU, double tWND, double tWNU,
					  double fracI, double fracJ, double fracK)
{
    return (fracI * fracJ * fracK * tESD +
		    fracI * fracJ * (1-fracK) * tESU +
			fracI * (1-fracJ) * fracK * tEND +
			fracI * (1-fracJ) * (1-fracK) * tENU +
			(1-fracI) * fracJ * fracK * tWSD +
		    (1-fracI) * fracJ * (1-fracK) * tWSU +
			(1-fracI) * (1-fracJ) * fracK * tWND +
			(1-fracI) * (1-fracJ) * (1-fracK) * tWNU);
} //endfunc


//______________________________________________________
//
//	InterpolateT
//		Calculates interpolated temperatures anywhere in a node
//		See LAB notes p.1069 for truth tables
//		Expensive calculation-- used only in slush movement and
//		interface checks!
//______________________________________________________
double InterpolateT (int i, int j, int k, double iFrac, double jFrac, double kFrac)
{
	int iE, iW, jN, jS, kD, kU;	

    if (!Geometry.isThick[DIM_I])
		iE = iW = I_FIRST;

	else			// has multiple i-nodes
	{
		iW = i - ((iFrac <= 0.5) && (i > I_FIRST));
		iE = i + ((iFrac > 0.5) && (i < I_LAST-1));

		iFrac = FractionCalc(iFrac, DelX, IJKToIndex(iW, j, k), IJKToIndex(iE, j, k));
	}

    if (!Geometry.isThick[DIM_J])
		jS = jN = J_FIRST;

	else			// has multiple j-nodes
	{
		jN = j - ((jFrac <= 0.5) && (j > J_FIRST));
		jS = j + ((jFrac > 0.5) && (j < J_LAST-1));

		if ( (jS == Laser.beginAbsorption[i][k]) && 
				(jS < J_LAST-1) && (jFrac <= 0.5) )	//need to extrapolate ?
		{
			if (jS == J_FIRST)  //topmost node
			{
				jS = jN + 1;
				jFrac -= 1;
			}

			else if (K[i][jN][k] < K_ISVACUUM)		//vacuum node above
			{
				jN = jS;		//shift to lower nodes
				jS = jN + 1;	//shift to lower nodes
				jFrac -= 1;
			}	//endif

		} //endif

		jFrac = FractionCalc(jFrac, DelY, IJKToIndex(i, jN, k), IJKToIndex(i, jS, k));
	}

	if (!Geometry.isThick[DIM_K])
		kD = kU = K_FIRST;

	else			// has multiple k-nodes
	{
		kU = k - ((kFrac <= 0.5) && (k > K_FIRST));
		kD = k + ((kFrac > 0.5) && (k < K_LAST-1));

		kFrac = FractionCalc(kFrac, DelZ, IJKToIndex(i, j, kU), IJKToIndex(i, j, kD));
	}
	
	return(Interpolate3D(T[iE][jS][kD], T[iE][jS][kU], T[iE][jN][kD], T[iE][jN][kU],
						T[iW][jS][kD], T[iW][jS][kU], T[iW][jN][kD], T[iW][jN][kU],
    					iFrac, jFrac, kFrac ) );
} //endfunc


//______________________________________________________
//
//	InterpolateT
//		Interpolates to edges of node.  This is also expensive
//		and should be used only for node-node catalysis
//		it cannot handle higher-order (200 etc.) directions
//______________________________________________________
double InterpolateT(const int nodeA, const int dirEdge)
{
	int i, j, k;
	double x;

	IndexToIJK(nodeA, i, j, k);		//rebuild ijk indices !!

	x = InterpolateT(i, j, k, 
		Interface.iCrystal[dirEdge] * 0.5 + 0.5,
		Interface.jCrystal[dirEdge] * 0.5 + 0.5,
		Interface.kCrystal[dirEdge] * 0.5 + 0.5);
	return(x);
}

//___________________________________________________
// KNeighbor
//        Calculate thermal resistance wrt neighbors
//___________________________________________________
inline double KNeighbor(int nodeA, int nodeB, DMATRIX kTherm, DMATRIX del1, DMATRIX del2)
{
	if ((A(kTherm) == 0.0) || (B(kTherm) == 0.0)) 
		return(0.0);
	
	else 
		return(A(del2) * B(kTherm) * A(kTherm) /
	      		(A(del1) * B(kTherm) + B(del1) * A(kTherm)) ); 
} //endinline


//___________________________________________________
// ThermalInit
//        Set up thermal related matrices etc.
//___________________________________________________
void ThermalInit ()
{
  	int p;
	int i, j, k;

	double tSub;

   	if (!InRange( Sim.dTemperatureLo, 0, Sim.eraMax, 0, MAX_TCHANGE))
    	ErrorMsg("TEMPERATURE_CHANGE_LO is out of range, or not specified");	

   	if (!InRange( Sim.dTemperatureHi, 0, Sim.eraMax, 0, MAX_TCHANGE))
    	ErrorMsg("TEMPERATURE_CHANGE_HI is out of range, or not specified");	

	for (p=P_FIRST; p<Sim.numberPhases; p++)
	{

		if (!InRange (Phase[p]->thermConductivity, MIN_DEGREES, MAX_DEGREES, K_MIN, K_MAX))
			ErrorMsg ("Thermal conductivity out of range {%g,%g} for %s %s",
        		K_MIN, K_MAX, Phase[p]->matlClassName, Phase[p]->phaseName);

		if (!InRange (Phase[p]->heatCapacity, MIN_DEGREES, MAX_DEGREES, CP_MIN, CP_MAX))
    		ErrorMsg ("Heat capacity out of range {%g,%g} for %s %s",
        		CP_MIN, CP_MAX, Phase[p]->matlClassName, Phase[p]->phaseName);
	} //endloop p

	
	//________ CHECK THERMAL PROPERTIES _______________
	tSub = Sim.tSubstrate.Evaluate(0);	//will fix later

	if (!InRange (tSub, (double)MIN_DEGREES, (double)MAX_DEGREES))
		ErrorMsg ("TEMP_SUBSTRATE out of range {%d, %d}", MIN_DEGREES, MAX_DEGREES);

  	MatrixNew (&T);			//set up temperature matrix

	for (i = I_FIRST; i < I_LAST; i++)
		for (j = J_FIRST; j < J_LAST; j++)
			for (k = K_FIRST; k < K_LAST; k++)
      			T[i][j][k] = tSub;

	Sim.tEnvironment = tSub;

  	CalcKNode (MatrixNew(&K), T);	//set up thermal conductivity matrix

	return;
} //endfunc
