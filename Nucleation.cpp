//  _________________________________________________________
// |
// |   Nucleate.cpp    Nucleation probability calculation
// |
// |   (C) 1998-99  Columbia University, MSME
// |   (C) 2000-01	Harvard University, DEAS
// |__________________________________________________________

#include "3dns.h"
#include "phase.h"
#include "parsefunc.h"		       //for LOOK_UP
#include "report.h"

#define EXT_LEVEL
#include "nucleation.h"
#undef EXT_LEVEL

//___________________________________
// Private functions
//___________________________________
bool RandomCheck(double meanNum);

//___________________________________
// Private variables
//___________________________________
DMATRIX thresholdHom;	  //nucleation threshold for heterogeneous nucleation to solid
DMATRIX thresholdHet;     //nucleation threshold for homogeneous nucleation to solid
DMATRIX thresholdHomLiquid;	  //nucleation threshold for heterogeneous nucleation to liquid
DMATRIX thresholdHetLiquid;     //nucleation threshold for homogeneous nucleation to liquid
DMATRIX thresholdHetLiquidSurface;     //nucleation threshold for homogeneous nucleation to liquid
double *randT3, *randT2, *randT1, *randT0;  //pre-calculated random values


//______________________________________________________
//
//	NucleateHeterogeneous
//		Check if heterogeneous nucleation will occur in a
//		liquid node with a neighbor from an ALIEN material class
//______________________________________________________
void NucleateHeterogeneous(CELL &cellNew, CELL &cellOld, unsigned int nodeLiquid, unsigned int nodeOther, const int dirToLiq, const double tInterface)
{
	bool willNucleate = false;
	int pIndex, pSolid;
	double densityNuc;		//nucleation probability density

	const double UNIT_CONV = CM_TO_M * CM_TO_M;

	if (!LIQ(CanHetNucleate)	 ||			//non-nucleating
		!OTHER(CatalyzeFreezing) ||			//non-catalytic interface
		(BZONE(dirToLiq) > 0) )				//non-(100) neighbor

		return;
									//loop thru all possible phases
	for (pIndex = 0; (!willNucleate &&
		(pIndex < Phase[LIQ(PhaseCodeLiquid)]->numTransitions)); pIndex++)
	{		
		if (tInterface > Phase[LIQ(PhaseCodeLiquid)]->transitionTemp[pIndex])  //no possibility of nucleation
			continue;

		pSolid = Phase[LIQ(PhaseCodeLiquid)]->transitionTo[pIndex];

		densityNuc = Sim.sClock.curDTime * LOOK_UP(Phase[pSolid]->hetNucleation, tInterface);

		//_______________ STOCHASTIC HETEROGENEOUS ____________________
		if (Sim.modeStochastic)
		{				
			switch(dirToLiq) //solid to liquid direction
			{
				case EAST:
				case WEST:
					densityNuc *= LIQ(AreaYZ) * UNIT_CONV;	//area converted to m^2 
					break;
				case SOUTH:
				case NORTH:
					densityNuc *= LIQ(AreaXZ) * UNIT_CONV;	//area converted to m^2 
					break;
				case BACK:
				case FRONT:
					densityNuc *= LIQ(AreaXY) * UNIT_CONV;	//area converted to m^2 
					break;

				default:		//no other directions are considered
					ErrorMsg("Illegal direction call in NucleateHeterogeneous");
			}; //endswitch
	
			LIQ(DensityHetNuc) = densityNuc;		//reset nucleation probability density

			willNucleate = RandomCheck(LIQ(DensityHetNuc) * -1.0);
		} //endif

		//_______________ DETERMINISTIC HETEROGENEOUS ____________________
		else
		{
			LIQ(DensityHetNuc) += densityNuc;		//increment nucleation density

			willNucleate = 	(LIQ(DensityHetNuc) > LIQ(thresholdHet));	// #/m2 density 
		}; //endif

		//________________ NUCLEATE IF NECESSARY _________________
		if (willNucleate)
		{
			NodeSolidify(cellNew, cellOld, nodeLiquid, dirToLiq | HETEROGENEOUS, 
				tInterface, pIndex, 0.0, ++Sim.numberNucleated);
			
			LIQ(Velocity) = 0.0;	//force zero velocity initially
			
			ReportNucleationEvent(nodeLiquid, dirToLiq | HETEROGENEOUS);
		}; //endif
	}; //endloop 
	
	return;
}; //endfunc


//______________________________________________________
//
//	NucleateHeterogeneousLiquid
//		Check if heterogeneous nucleation will occur in a
//		liquid node with a neighbor from an ALIEN material class
//______________________________________________________
void NucleateHeterogeneousLiquid(CELL &cellNew, CELL &cellOld, unsigned int nodeSolid, unsigned int nodeOther, const int dirToSolid, const double tInterface)
{
	bool willNucleate = false;
	int pIndex, pSolid;
	double densityNuc;		//nucleation probability density

	const double UNIT_CONV = CM_TO_M * CM_TO_M;

	if (!SOL(CanHetNucleateLiquid)	 ||			//non-nucleating
		!OTHER(CatalyzeMelting) ||			//non-catalytic interface
		(BZONE(dirToSolid) > 0) )				//non-(100) neighbor

		return;
									//loop thru all possible phases
	
	if (SOL(T) < *(Phase[SOL(cellOld.cPhaseSolid)]->transitionTemp))
		return;

	densityNuc = Sim.sClock.curDTime * LOOK_UP(Phase[SOL(cellOld.cPhaseSolid)]->hetNucleationLiquid, tInterface);

	//_______________ STOCHASTIC HETEROGENEOUS ____________________
	if (Sim.modeStochastic)
	{				
		switch(dirToSolid) //solid to liquid direction
		{
			case EAST:
			case WEST:
				densityNuc *= SOL(AreaYZ) * UNIT_CONV;	//area converted to m^2 
				break;
			case SOUTH:
			case NORTH:
				densityNuc *= SOL(AreaXZ) * UNIT_CONV;	//area converted to m^2 
				break;
			case BACK:
			case FRONT:
				densityNuc *= SOL(AreaXY) * UNIT_CONV;	//area converted to m^2 
				break;

			default:		//no other directions are considered
				ErrorMsg("Illegal direction call in NucleateHeterogeneous");
		}; //endswitch
	
		SOL(DensityHetNucLiquid) = densityNuc;		//reset nucleation probability density

		willNucleate = RandomCheck(SOL(DensityHetNucLiquid) * -1.0);
	} //endif

	//_______________ DETERMINISTIC HETEROGENEOUS ____________________
	else
	{
		SOL(DensityHetNucLiquid) += densityNuc;		//increment nucleation density

		willNucleate = 	(SOL(DensityHetNucLiquid) > SOL(thresholdHetLiquid));	// #/m2 density 
	}; //endif

	//________________ NUCLEATE IF NECESSARY _________________
	if (willNucleate)
	{
		//To do: incorporate NodeMelt like from free surface
		//NodeSolidify(cellNew, cellOld, nodeSolid, dirToSolid | HETEROGENEOUS, 
		//	tInterface, pIndex, 0.0, ++Sim.numberNucleated);
		NodeMelt(cellNew, cellOld, nodeSolid, Interface.dirInverse[dirToSolid] | HETEROGENEOUS_LIQUID, tInterface, 0.0);	
		SOL(Velocity) = 0.0;	//force zero velocity initially
			
		ReportNucleationEvent(nodeSolid, dirToSolid | HETEROGENEOUS);
	}; //endif
	
	return;
}; //endfunc

//______________________________________________________
//
//	NucleateHeterogeneousLiquid
//		Check if heterogeneous nucleation will occur in a
//		liquid node with a neighbor from an ALIEN material class
//______________________________________________________
void NucleateHeterogeneousLiquidSurface(CELL &cellNew, CELL &cellOld, unsigned int nodeSolid, const int dirToSolid, const double tInterface)
{
	bool willNucleate = false;
	int pIndex, pSolid;
	double densityNuc;		//nucleation probability density

	const double UNIT_CONV = CM_TO_M * CM_TO_M;

	if (!SOL(CanHetNucleateLiquidSurface)	 ||			//non-nucleating
		(BZONE(dirToSolid) > 0) )				//non-(100) neighbor

		return;
									//loop thru all possible phases
	
	if (SOL(T) < *(Phase[SOL(cellOld.cPhaseSolid)]->transitionTemp))
		return;

	densityNuc = Sim.sClock.curDTime * LOOK_UP(Phase[SOL(cellOld.cPhaseSolid)]->hetNucleationLiquidSurface, tInterface);
	double* nucleationrate = Phase[SOL(cellOld.cPhaseSolid)]->hetNucleationLiquidSurface;

	//_______________ STOCHASTIC HETEROGENEOUS ____________________
	if (Sim.modeStochastic)
	{				
		switch(dirToSolid) //solid to liquid direction
		{
			case EAST:
			case WEST:
				densityNuc *= SOL(AreaYZ) * UNIT_CONV;	//area converted to m^2 
				break;
			case SOUTH:
			case NORTH:
				densityNuc *= SOL(AreaXZ) * UNIT_CONV;	//area converted to m^2 
				break;
			case BACK:
			case FRONT:
				densityNuc *= SOL(AreaXY) * UNIT_CONV;	//area converted to m^2 
				break;

			default:		//no other directions are considered
				ErrorMsg("Illegal direction call in NucleateHeterogeneous");
		}; //endswitch
	
		SOL(DensityHetNucLiquidSurface) = densityNuc;		//reset nucleation probability density

		willNucleate = RandomCheck(SOL(DensityHetNucLiquidSurface) * -1.0);
	} //endif

	//_______________ DETERMINISTIC HETEROGENEOUS ____________________
	else
	{
		SOL(DensityHetNuc) += densityNuc;		//increment nucleation density

		willNucleate = 	(SOL(DensityHetNucLiquidSurface) > SOL(thresholdHetLiquidSurface));	// #/m2 density 
	}; //endif

	//________________ NUCLEATE IF NECESSARY _________________
	if (willNucleate)
	{
		//To do: incorporate NodeMelt like from free surface
		//NodeSolidify(cellNew, cellOld, nodeSolid, dirToSolid | HETEROGENEOUS, 
		//	tInterface, pIndex, 0.0, ++Sim.numberNucleated);
		NodeMelt(cellNew, cellOld, nodeSolid, Interface.dirInverse[dirToSolid] | HETEROGENEOUS_LIQUID, tInterface, 0.0);	
		SOL(Velocity) = 0.0;	//force zero velocity initially
			
		ReportNucleationEvent(nodeSolid, dirToSolid | HETEROGENEOUS);
	}; //endif
	
	return;
}; //endfunc

//______________________________________________________
//
//	NucleateHomogeneous
//		Check if homogeneous nucleation will occur
//		and if so, initiate solidification
//______________________________________________________
void NucleateHomogeneous(CELL &cellNew, CELL &cellOld)
{
	bool willNucleate = false;
    int i, j, k;
	int pIndex, pSolid;
	int nodeA;
	double densityNuc;		//nucleation probability density

	const double CM_TO_M_3 = CM_TO_M * CM_TO_M * CM_TO_M;
    
  	if(Sim.numberLiquid == 0)
    	return;

	for (i = I_FIRST; i < I_LAST; i++)	
	{
		if (!(*(Geometry.canChangeI + i)) )
			continue;		//this page cant change

		for (j = J_FIRST; j < J_LAST; j++)
		{
			if (!(*(Geometry.canChangeJ + j)) )
				continue;		//this row cant change
			
			for (k = K_FIRST; k < K_LAST; k++)
			{	
				if (!(*(Geometry.canChangeK + k)) )
					continue;		//this column cant change

				nodeA = IJKToIndex(i, j, k);
				
				if ((A(cellOld.cState) != LIQUID) || !A(CanHomNucleate))		//only liquid nodes can nucleate
					continue;

				if (Sim.modeStochastic)		//reset probability for stochastic mode
					A(DensityHomNuc) = (double)0.0;

				for (pIndex = 0; (!willNucleate && 
					(pIndex < Phase[A(cellOld.cPhaseLiquid)]->numTransitions)); pIndex++)
				{
					if (A(T) > Phase[A(cellOld.cPhaseLiquid)]->transitionTemp[pIndex])
						continue;
					
					pSolid = Phase[A(cellOld.cPhaseLiquid)]->transitionTo[pIndex];

					densityNuc = Sim.sClock.curDTime * LOOK_UP(Phase[pSolid]->homNucleation, A(T));

					//_______________ STOCHASTIC HOMOGENEOUS ____________________
					if (Sim.modeStochastic)
					{	
						A(DensityHomNuc) = densityNuc * A(Volume) * CM_TO_M_3;
						
						willNucleate = RandomCheck(A(DensityHomNuc) * -1.0);
					} //endif

					//_______________ DETERMINISTIC HOMOGENEOUS ____________________
					else
					{
						A(DensityHomNuc) += densityNuc;		//increment nucleation density

						willNucleate = 	( A(DensityHomNuc) > A(thresholdHom) );	//heterogeneous
					}; //endif
					
					//_______________ NUCLEATE IF NECESSARY _____________
					if (willNucleate)
					{
						NodeSolidify(cellNew, cellOld, nodeA, HOMOGENEOUS, A(T), 0, 0.0, ++Sim.numberNucleated);
					
						ReportNucleationEvent(nodeA, HOMOGENEOUS);
					}; //endif
				}; //endloop p
	
			}; //endloop k
		}; //endloop j
	}; //endloop i
	
	return;
}; //endfunc


//______________________________________________________
//
//	NucleateHomogeneousLiquid
//		Check if homogeneous nucleation will occur
//		and if so, initiate solidification
//______________________________________________________
void NucleateHomogeneousLiquid(CELL &cellNew, CELL &cellOld)
{
	bool willNucleate = false;
    int i, j, k;
	int pIndex, pSolid;
	int nodeA;
	double densityNuc;		//nucleation probability density

	const double CM_TO_M_3 = CM_TO_M * CM_TO_M * CM_TO_M;
    
  	if(Sim.numberSolid == 0)
    	return;

	for (i = I_FIRST; i < I_LAST; i++)	
	{
		if (!(*(Geometry.canChangeI + i)) )
			continue;		//this page cant change

		for (j = J_FIRST; j < J_LAST; j++)
		{
			if (!(*(Geometry.canChangeJ + j)) )
				continue;		//this row cant change
			
			for (k = K_FIRST; k < K_LAST; k++)
			{	
				if (!(*(Geometry.canChangeK + k)) )
					continue;		//this column cant change

				nodeA = IJKToIndex(i, j, k);
				
				if ((A(cellOld.cState) != SOLID) || !A(CanHomNucleateLiquid))		//only solid nodes can nucleate
					continue;

				if (Sim.modeStochastic)		//reset probability for stochastic mode
					A(DensityHomNucLiquid) = (double)0.0;

				if (A(T) < *(Phase[A(cellOld.cPhaseSolid)]->transitionTemp))
					continue;

				densityNuc = Sim.sClock.curDTime * LOOK_UP(Phase[A(cellOld.cPhaseSolid)]->homNucleationLiquid, A(T));

				//_______________ STOCHASTIC HOMOGENEOUS ____________________
				if (Sim.modeStochastic)
				{	
					A(DensityHomNucLiquid) = densityNuc * A(Volume) * CM_TO_M_3;
						
					willNucleate = RandomCheck(A(DensityHomNucLiquid) * -1.0);
				} //endif

				//_______________ DETERMINISTIC HOMOGENEOUS ____________________
				else
				{
					A(DensityHomNuc) += densityNuc;		//increment nucleation density

					willNucleate = 	( A(DensityHomNucLiquid) > A(thresholdHomLiquid) );	//heterogeneous
				}; //endif
					
				//_______________ NUCLEATE IF NECESSARY _____________
				if (willNucleate)
				{
					//To incoporate NodeMelting
					//NodeSolidify(cellNew, cellOld, nodeA, HOMOGENEOUS, A(T), 0, 0.0, ++Sim.numberNucleated);
					NodeMelt(cellNew, cellOld, nodeA, HETEROGENEOUS_LIQUID, A(T), 0.0);
					//To incoporate liquid nucleation reporting
					ReportNucleationEvent(nodeA, HOMOGENEOUS);
				}; //endif
	
			}; //endloop k
		}; //endloop j
	}; //endloop i
	
	return;
}; //endfunc

//______________________________________________________
//
//      NucleationCleanup
//       Increments nucleation probabilities in each node
//______________________________________________________
void NucleationCleanup()
{
  	MatrixFree(CatalyzeFreezing);
  	MatrixFree(CatalyzeMelting);
  	MatrixFree(DensityHetNuc);
    MatrixFree(DensityHomNuc);
	MatrixFree(DensityHetNucLiquid);
    MatrixFree(DensityHomNucLiquid);

  	MatrixFree(CanHetNucleate);
  	MatrixFree(CanHomNucleate);
	MatrixFree(CanHetNucleateLiquid);
  	MatrixFree(CanHomNucleateLiquid);

	MatrixFree(thresholdHet);
	MatrixFree(thresholdHom);
	MatrixFree(thresholdHetLiquid);
	MatrixFree(thresholdHomLiquid);
	MatrixFree(thresholdHetLiquidSurface);
}; //endfunc


//______________________________________________________
//
//      NucleationInit
//       Increments nucleation probabilities in each node
//______________________________________________________
void NucleationInit()
{
	int i, j, k, p, r;
	int nodeA;
	int phaseCode;
	int typeVal;

	string typeMsg;

  	MatrixNew(&CanHetNucleate);
  	MatrixNew(&CanHomNucleate);
	MatrixNew(&CanHetNucleateLiquid);
  	MatrixNew(&CanHomNucleateLiquid);
	MatrixNew(&CanHetNucleateLiquidSurface);
  	MatrixNew(&thresholdHet);
    MatrixNew(&thresholdHom);
	MatrixNew(&thresholdHetLiquid);
    MatrixNew(&thresholdHomLiquid);
	MatrixNew(&thresholdHetLiquidSurface);
  	MatrixZero(MatrixNew(&DensityHetNuc));
  	MatrixZero(MatrixNew(&DensityHomNuc));
	MatrixZero(MatrixNew(&DensityHetNucLiquid));
  	MatrixZero(MatrixNew(&DensityHomNucLiquid));
	MatrixZero(MatrixNew(&DensityHetNucLiquidSurface));

	//_________ CHECK REGION INFORMATION _______________
	for (r=0; r<Geometry.numberRegions; r++)
	{
		typeMsg = (r < Geometry.jZones) ? "LAYER" : "OVERLAY";
		typeVal = (r < Geometry.jZones) ? r : r - Geometry.jZones;
		
		if (Region[r]->canChange)
		{
			if (!InRange (Region[r]->thresholdHet, (double)0.0, 1.001 * DOUBLE_INFINITY))
    			ErrorMsg ("HET_THRESHOLD out of range in %s %d", typeMsg, typeVal);

			if (!InRange (Region[r]->thresholdHetLiquid, (double)0.0, 1.001 * DOUBLE_INFINITY))
    			ErrorMsg ("HET_THRESHOLD_LIQUID out of range in %s %d", typeMsg, typeVal);

			if (!InRange (Region[r]->thresholdHetLiquidSurface, (double)0.0, 1.001 * DOUBLE_INFINITY))
    			ErrorMsg ("HET_THRESHOLD_LIQUID_SURFACE out of range in %s %d", typeMsg, typeVal);

			if (!InRange (Region[r]->thresholdHom, (double)0.0, 1.001 * DOUBLE_INFINITY))
    			ErrorMsg ("HOM_THRESHOLD out of range in %s %d", typeMsg, typeVal);

			if (!InRange (Region[r]->thresholdHomLiquid, (double)0.0, 1.001 * DOUBLE_INFINITY))
    			ErrorMsg ("HOM_THRESHOLD_LIQUID out of range in %s %d", typeMsg, typeVal);

								//must have finite thresholds to nucleate
			Region[r]->canHetNucleate = (Region[r]->thresholdHet < DOUBLE_INFINITY);
			Region[r]->canHomNucleate = (Region[r]->thresholdHom < DOUBLE_INFINITY);
			Region[r]->canHetNucleateLiquid = (Region[r]->thresholdHetLiquid < DOUBLE_INFINITY);
			Region[r]->canHomNucleateLiquid = (Region[r]->thresholdHomLiquid < DOUBLE_INFINITY);
			Region[r]->canHetNucleateLiquidSurface = (Region[r]->thresholdHetLiquidSurface < DOUBLE_INFINITY);

			for (p=0; p < Phase[Region[r]->phaseLiquid]->numTransitions; p++)
			{
				phaseCode = *(Phase[Region[r]->phaseLiquid]->transitionTo+p);
        		
				if (Region[r]->canHetNucleate)
				{
					if (Phase[phaseCode]->hetNucleation == NULL)
						ErrorMsg("NUCLEATION_HET not defined for %s %s", 
							Phase[phaseCode]->matlClassName, Phase[phaseCode]->phaseName);

					if (Phase[phaseCode]->hetNucleationLiquid == NULL)
						ErrorMsg("NUCLEATION_HET_LIQUID not defined for %s %s", 
							Phase[phaseCode]->matlClassName, Phase[phaseCode]->phaseName);

					if (!InRange (Phase[phaseCode]->hetNucleation, MIN_DEGREES, 
							*(Phase[phaseCode]->transitionTemp), 0.0, HET_RATE_MAX))
    					ErrorMsg("Heterogeneous soild nucleation rate out of range in %s %s",
							Phase[phaseCode]->matlClassName, Phase[phaseCode]->phaseName);

					if (!InRange (Phase[phaseCode]->hetNucleationLiquid, *(Phase[phaseCode]->transitionTemp), 
							MAX_DEGREES, 0.0, HET_RATE_MAX))
    					ErrorMsg("Heterogeneous liquid nucleation rate out of range in %s %s",
							Phase[phaseCode]->matlClassName, Phase[phaseCode]->phaseName);

					if (!InRange (Phase[phaseCode]->hetNucleation, *(Phase[phaseCode]->transitionTemp), 
							MAX_DEGREES, 0.0, NUCRATE_MIN))
    					ErrorMsg("Heterogeneous nucleation rate must be zero above transition temperature in %s %s",
							Phase[phaseCode]->matlClassName, Phase[phaseCode]->phaseName);

					if (!InRange (Phase[phaseCode]->hetNucleationLiquid, MIN_DEGREES, 
							*(Phase[phaseCode]->transitionTemp), 0.0, NUCRATE_MIN))
    					ErrorMsg("Heterogeneous nucleation rate must be zero above transition temperature in %s %s",
							Phase[phaseCode]->matlClassName, Phase[phaseCode]->phaseName);
				}; //endif

				if (Region[r]->canHomNucleate)
				{
        			if (Phase[phaseCode]->homNucleation == NULL)
						ErrorMsg("NUCLEATION_HOM not defined for %s %s", 
							Phase[phaseCode]->matlClassName, Phase[phaseCode]->phaseName);

					if (Phase[phaseCode]->homNucleationLiquid == NULL)
						ErrorMsg("NUCLEATION_HOM_LIQUID not defined for %s %s", 
							Phase[phaseCode]->matlClassName, Phase[phaseCode]->phaseName);

					if (!InRange (Phase[phaseCode]->homNucleation, MIN_DEGREES, 
							*(Phase[phaseCode]->transitionTemp), 0.0, HOM_RATE_MAX))
    					ErrorMsg("Homogeneous solid nucleation rate out of range in %s %s",
							Phase[phaseCode]->matlClassName, Phase[phaseCode]->phaseName);

					if (!InRange (Phase[phaseCode]->homNucleationLiquid, *(Phase[phaseCode]->transitionTemp), 
							MAX_DEGREES, 0.0, HOM_RATE_MAX))
    					ErrorMsg("Homogeneous liquid nucleation rate out of range in %s %s",
							Phase[phaseCode]->matlClassName, Phase[phaseCode]->phaseName);

					if (!InRange (Phase[phaseCode]->homNucleation, *(Phase[phaseCode]->transitionTemp), 
							MAX_DEGREES, 0.0, NUCRATE_MIN))
    					ErrorMsg("Homogeneous solid nucleation rate must be zero above transition temperature in %s %s",
							Phase[phaseCode]->matlClassName, Phase[phaseCode]->phaseName);

					if (!InRange (Phase[phaseCode]->homNucleationLiquid, MIN_DEGREES, 
							*(Phase[phaseCode]->transitionTemp), 0.0, NUCRATE_MIN))
    					ErrorMsg("Homogeneous liquid nucleation rate must be zero above transition temperature in %s %s",
							Phase[phaseCode]->matlClassName, Phase[phaseCode]->phaseName);
				}; //endif
			}; //endloop p;
		}; //endif

		for (i = 0; i < Region[r]->numberI; i++)
		{
			for (j = 0; j < Region[r]->numberJ; j++)
			{
				for (k = 0; k < Region[r]->numberK; k++)
				{
					nodeA = IJKToIndex(Region[r]->iLocations[i], 
							Region[r]->jLocations[j], Region[r]->kLocations[k]);
					
					A(CanHetNucleate) = Region[r]->canHetNucleate;
					A(CanHomNucleate) = Region[r]->canHomNucleate;
					A(thresholdHet)	= Region[r]->thresholdHet;
					A(thresholdHom) = Region[r]->thresholdHom;

					A(CanHetNucleateLiquid) = Region[r]->canHetNucleateLiquid;
					A(CanHomNucleateLiquid) = Region[r]->canHomNucleateLiquid;
					A(CanHetNucleateLiquidSurface) = Region[r]->canHetNucleateLiquidSurface;
					A(thresholdHetLiquid) = Region[r]->thresholdHetLiquid;
					A(thresholdHomLiquid) = Region[r]->thresholdHomLiquid;
					A(thresholdHetLiquidSurface) = Region[r]->thresholdHetLiquidSurface;
				}; //endloop k
			}; //endloop j
		}; //endloop i

	}; //endloop r

	randT3 = new double[RAND_MAX+2];
	randT2 = new double[RAND_MAX+2];
	randT1 = new double[RAND_MAX+2];
	randT0 = new double[RAND_MAX+2];

	if ((randT3 == NULL) || (randT2 == NULL) ||
			(randT1 == NULL) || (randT0 == NULL) )
		ErrorMsg("Out of memory");

	const double R3 = 1.0 / (double)(RAND_MAX+1);
	const double R2 = pow(R3,2);
	const double R1 = pow(R3,3);
	const double R0 = pow(R3,4);
    
	for (i = 0; i <= RAND_MAX+1; i++)
	{
		if (i > 0)
			*(randT3+i) = log((double)i * R3);

		*(randT2+i) = (double)i * R2;
		*(randT1+i) = (double)i * R1;
		*(randT0+i) = (double)i * R0;

	}; //endloop i

	*(randT3) = -DOUBLE_INFINITY;

	if (Sim.modeStochastic)
	{
		if (Sim.seedStochastic == INTFLAG_ALL)
			Sim.seedStochastic = int(time(NULL));

		else if (Sim.seedStochastic <= 0)
			ErrorMsg("Only positive integer values (or 'All') are allowed for SEED_STOCHASTIC");

		srand((unsigned)Sim.seedStochastic);	//seed the generator
		InfoMsg("\nStochastic seed is [%d] ", Sim.seedStochastic);
	}; //endif

	return;
};  //endfunc


inline double chkRange(double &randVal, double &testVal)
{
	if (*((double*)(&randVal) + 1) <= testVal)
		return (0.0);		//no chance

	return (testVal - randVal);
};


//______________________________________________________
//
//   RandomCheck
//      Checks nucleation probabilities against 
//		random number function to determine nucleation
//______________________________________________________
bool RandomCheck(double meanNum)
{
	int r3, r2, r1, r0;		//random numbers
	static double pDelta=0;	
	const double R3 = 1.0 / (double)(RAND_MAX + 1);
	const double MIN_ARG = -1.0e-14;		//minimum argument calculable in exp

	r3 = rand();						//level 3 check
	pDelta = chkRange(randT3[r3], meanNum);
	if (pDelta == 0.0) return(false);	//certainly not nucleation
	if (pDelta < 0.0) return (true);		//certainly nucleation

	pDelta = (meanNum > MIN_ARG) ? 
		(1.0 - R3*r3) + meanNum :		//true probability for nucleation
		exp(meanNum) - R3*r3;	  
	
	r2 = rand();					//level 2 check
	pDelta = chkRange(randT2[r2], pDelta);
	if (pDelta == 0.0) return(false);
	if (pDelta < 0.0) return (true);

	r1 = rand();					//level 1 check
	pDelta = chkRange(randT1[r1], pDelta);
	if (pDelta == 0.0) return(false);
	if (pDelta < 0.0) return (true);

	r0 = rand();					//level 0 check
	pDelta = chkRange(randT0[r0], pDelta);
	if (pDelta >= 0.0) return(false);

	return(true);		//it will nucleate
}; //endfunc	
