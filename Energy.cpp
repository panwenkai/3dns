//  _________________________________________________________
// |
// |   Energy.cpp	Energy deposition routines
// |
// |
// |   (C) 1998-99  Columbia University, MSME
// |   (C) 2000-01	Harvard University, DEAS
// |__________________________________________________________

#include "3dns.h"
#include "phase.h"
#include "report.h"

#define EXT_LEVEL
#include "energy.h"
#undef EXT_LEVEL

#define SQ(x) ((x)*(x))
#define EMISSIVITY_MIN ((double) 0.0)
#define EMISSIVITY_MAX ((double) 1.0)
#define INDEXN_MIN ((double) 0.0)
#define INDEXN_MAX ((double) 10.0)
#define INDEXK_MIN ((double) 0.0)
#define INDEXK_MAX ((double) 10.0)
#define LASER_THRESHOLD ((double) 1E-5)		//threshold for nonzero energy

//_________________________________________
// Private functions
//_________________________________________
int Compare (double *x, int size, double val);
void CreateEgyTemporal(ENTRY &entry);
double LaserEgy (int i);
double PropertyAverage(unsigned int nodeA, double *propSolid, double *propLiquid);
void Radiate();
double Transmit(int i, int k);

//_________________________________________
// Private variables
//_________________________________________
bool *isAbsorbing;	//if phase is absorbing at some temperature

MATRIX2D colEgyCurrent;
MATRIX2D colEgyTotal;


//______________________________________________________
//
// 	Absorb
//		Calculates energy absorbed in a single node with
//		energy incident at the top
//______________________________________________________
double Absorb(unsigned int nodeA, double egyIncident)
{
	double alpha;

	alpha = PropertyAverage(nodeA, Phase[A(PhaseCodeSolid)]->optAbsorption, 
				Phase[A(PhaseCodeLiquid)]->optAbsorption);

	return (egyIncident * (1 - exp(-alpha * A(DelY)) ));
}; //endfunc


//______________________________________________________
//
//	LaserInit
//		Initialize laser parameters
//______________________________________________________
void LaserInit()
{
	int i, k;
	int numOffsets=0;	//number of offsets beyond simulation time
	int phaseCode;

	double curIntensity;
	double timeOffset;
	double thisEgy, totEgy;
	double egyX, egyZ;		//temporary energy values

	string typeMsg;
	
	isAbsorbing = new bool[Sim.numberPhases]; //to find non-absorbing phases
	memset(isAbsorbing, 0, sizeof(bool) * Sim.numberPhases);

    if (!InRange(Laser.waveLength, M_TO_CM*MIN_LAMBDA, M_TO_CM*MAX_LAMBDA))
    	ErrorMsg("LASER_WAVELENGTH is out of range, or not specified");

    if (!InRange(Laser.energyDensity, MIN_ENERGY, MAX_ENERGY))
    	ErrorMsg("ENERGY_DENSITY is out of range");

	for (phaseCode = P_FIRST; phaseCode < P_LAST; phaseCode++)
	{
		if (!InRange (Phase[phaseCode]->emissivity, MIN_DEGREES, MAX_DEGREES,
    			EMISSIVITY_MIN, EMISSIVITY_MAX))
    		ErrorMsg ("EMISSIVITY out of range {%g,%g} for %s %s",
        		EMISSIVITY_MIN, EMISSIVITY_MAX, Phase[phaseCode]->matlClassName,
				Phase[phaseCode]->phaseName);

		if (!InRange (Phase[phaseCode]->indexN, MIN_DEGREES, MAX_DEGREES,
    			INDEXN_MIN, INDEXN_MAX))
    		ErrorMsg ("INDEX_N out of range {%g,%g} for %s %s",
        		INDEXN_MIN, INDEXN_MAX, Phase[phaseCode]->matlClassName,
				Phase[phaseCode]->phaseName);

		if (!InRange (Phase[phaseCode]->indexK, MIN_DEGREES, MAX_DEGREES,
    			INDEXK_MIN, INDEXK_MAX))
    		ErrorMsg ("INDEX_K out of range {%g,%g} for %s %s",
        		INDEXK_MIN, INDEXK_MAX,  Phase[phaseCode]->matlClassName,
				Phase[phaseCode]->phaseName);
		
		Phase[phaseCode]->optAbsorption = new double[MAX_DEGREES];

		for (i = MIN_DEGREES; i < MAX_DEGREES; i++)  //build absorption table
		{
			Phase[phaseCode]->optAbsorption[i] = Phase[phaseCode]->indexK[i] * 
					4.0 * PI / Laser.waveLength;

			isAbsorbing[phaseCode] |= (Phase[phaseCode]->indexK[i] > 0.0);
		}; //endloop i

	}; //endloop phaseCode

	MatrixZero(MatrixNew(&QLaser));				//laser energy deposited in each node

	colEgyCurrent.Reset(I_LAST, K_LAST, 0);
	colEgyTotal.Reset(I_LAST, K_LAST, 0);

	Laser.egySpatial.Reset(I_LAST, K_LAST, 0);  //initialize laser spatial array
	Laser.beginAbsorption = (int **)Create2DArray(I_LAST, K_LAST, sizeof(int));
	
	for (i = I_FIRST; i < I_LAST; i++)			//construct 2D spatial density function
	{
		for (k = K_FIRST; k < K_LAST; k++)
		{
			
			egyX = (Laser.dataSpatialX.dataType == T_UNKNOWN) ? 1.0 :
				Laser.dataSpatialX.Evaluate((Geometry.xMap[i] + 0.5 * DelX[i][J_FIRST][K_FIRST]) * CM_TO_M,
											(Geometry.zMap[k] + 0.5 * DelZ[I_FIRST][J_FIRST][k]) * CM_TO_M );

			egyZ = (Laser.dataSpatialZ.dataType == T_UNKNOWN) ? 1.0 :
				Laser.dataSpatialZ.Evaluate((Geometry.zMap[k] + 0.5 * DelZ[I_FIRST][J_FIRST][k]) * CM_TO_M);

			Laser.egySpatial[i][k] = egyX * egyZ;
		}; //endloop k
	}; //endloop i

	if (!Normalize(Laser.egySpatial, 0, 1))		//peak value must be 1
	{
		WarningMsg("Cant normalize laser SPATIAL profile, so using 0 everywhere... ");
		Laser.egySpatial.Reset(I_LAST, K_LAST, 0);
	}

	Floor(Laser.egySpatial, LASER_THRESHOLD, 0);	//drop values below LASER_THRESHOLD

	for (i = I_FIRST; i < I_LAST; i++)			//pre-calculate energy per node
	{
		for (k = K_FIRST; k < K_LAST; k++)
		{
			Laser.egySpatial[i][k] *= AreaXZ[i][J_FIRST][k];

			Laser.beginAbsorption[i][k] = 0;
			
			while (!isAbsorbing[PhaseCodeSolid[i][Laser.beginAbsorption[i][k]][k]] && 
					!CanChange[i][Laser.beginAbsorption[i][k]][k])
			{
				Laser.beginAbsorption[i][k]++;
				
				if (Laser.beginAbsorption[i][k] == J_LAST)
						break;
			}; //endwhile
		
		}; //endloop k
	}; //endloop i 
	
	//____________ TEMPORAL ___________________________________
	if (Laser.dataTemporal.dataType == T_UNKNOWN)
		ErrorMsg("LASER_TEMPORAL was not specified");

	totEgy = 0.0;		//total integrated energy in pulse
	Laser.peakTemporal = 0.0;	
	timeOffset = 0.0;	//offset time into waveform

	do		// full integration of temporal function
	{
		thisEgy = 0.0;

		do	// integrate over simulation time
		{			
			curIntensity = Laser.dataTemporal.Evaluate(Sim.sClock.curTime + timeOffset);
		
			if (!InRange(curIntensity, 0.0, 1.0/LASER_THRESHOLD))
				ErrorMsg("LASER_TEMPORAL out of range {%g,%g} at t=%11.5G sec", 
					0.0, 1.0/LASER_THRESHOLD, Sim.sClock.curTime + timeOffset );
		
			Laser.peakTemporal = max(Laser.peakTemporal, curIntensity);  //find peak value
			thisEgy += curIntensity * Sim.sClock.curDTime;		//newtons method
		
		} while ( Sim.sClock.Increment() ); 

		totEgy += thisEgy;				//add this section to total energy
		timeOffset += Sim.sClock.curTime;		//move offset
		numOffsets++;

		Sim.sClock.Reset();		//reset clock for more trials

		if (numOffsets == 3)
			WarningMsg("Temporal function is nonzero beyond simulation range ");

	} while ( (thisEgy > (totEgy * 1000.0 * LASER_THRESHOLD)) && //still appreciable energy
			(timeOffset < MAX_TIME) );	//still reasonable times

	if (timeOffset > 2 * MAX_TIME)
		ErrorMsg("Laser temporal profile is not normalizable, check that it is 0 at long times");

	if (totEgy < LASER_THRESHOLD * LASER_THRESHOLD * MAX_TIME)
	{
		WarningMsg("Laser temporal profile is considered to be 0.0");
		Laser.normTemporal = 0;
	}

	else 
		Laser.normTemporal = 1 / totEgy;		//normalization factor for temporal

	return;
}; //endfunc


//______________________________________________________
//
//      LaserCleanup
//       Frees all arrays associated with the laser
//______________________________________________________
void LaserCleanup ()
{
	MatrixFree (QLaser);

	delete [] isAbsorbing;

	return;
}; //endfunc


//______________________________________________________
//
//      LaserInput
//      Distributes laser energy into the sample via
//              the MATRIX QLaser
//______________________________________________________
void LaserInput()
{
	int i, j, k;

  	double egyTemporal;
	double egyRemaining;
	double egyThisNode;

	unsigned int nodeA;

  	MatrixZero(QLaser);
  	
	egyTemporal = Laser.dataTemporal.Evaluate(Sim.sClock.curTime);  //evaluate temporal function

	if (egyTemporal < (Laser.peakTemporal * LASER_THRESHOLD))	//laser is effectively off
		return;
	
	egyTemporal *= Sim.sClock.curDTime * Laser.normTemporal * Laser.energyDensity;    // [Joules/cm2] this clock

	for (i = I_FIRST; i < I_LAST; i++)
	{
		for (k = K_FIRST; k < K_LAST; k++)
		{
			if ((Laser.egySpatial[i][k] == 0.0) || (Laser.beginAbsorption[i][k] == J_LAST))
				continue;							//no laser energy here

											//calculate amount transmitted into sample
			colEgyCurrent[i][k] = egyTemporal * Laser.egySpatial[i][k] * Transmit(i,k);
			egyRemaining = colEgyCurrent[i][k];

			j = Laser.beginAbsorption[i][k];

			while ((j < J_LAST) && 
				(egyRemaining > (colEgyCurrent[i][k] * LASER_THRESHOLD)))
			{
				nodeA = IJKToIndex(i, j, k);

				egyThisNode = Absorb(nodeA, egyRemaining);
				egyRemaining -= egyThisNode;

				A(QLaser) += egyThisNode;
				j++;
			};
				
		}; //endloop k
	}; //endloop i

	if (Sim.modeRadiation) 
		Radiate();	//negative flux as radiation from the surface

}; //endfunc
 

//______________________________________________________
//
// 	PropertyAverage
// 		Computes average properties in mixed-phase nodes
//______________________________________________________
double PropertyAverage(unsigned int nodeA, double *propSolid, double *propLiquid)
{
	switch (A(State))
	{
		case SOLID:
			return (LOOK_UP (propSolid, A(T)) );
			break;
		case LIQUID:
			return (LOOK_UP (propLiquid, A(T)) );
			break;
		default:	//slush
			return (LOOK_UP (propSolid, A(T)) * A(FractionSolid) +
					LOOK_UP (propLiquid, A(T)) * (1 - A(FractionSolid)) );
	}; //endswitch

	return(0.0);
}; //endfunc


//______________________________________________________
//
// 	Transmit
// 	Calculates the fraction of energy transmitted into the first absorbing
//	node.  Thin film effects for a single non-absorbing dielectric are treated
//	nodeA is first dielectric node, nodeB is first absorbing node
//______________________________________________________
double Transmit(int i, int k)
{
  	double t, rad;
	double n2, n3, k3;
	double sqN2, sqN3, sqK3;
	
	unsigned int nodeA, nodeB;

	nodeB = IJKToIndex(i, Laser.beginAbsorption[i][k], k);

	n3 = PropertyAverage(nodeB, Phase[B(PhaseCodeSolid)]->indexN, 
				Phase[B(PhaseCodeLiquid)]->indexN);		//index in absorbing medium
	k3 = PropertyAverage(nodeB, Phase[B(PhaseCodeSolid)]->indexK,
				Phase[B(PhaseCodeLiquid)]->indexK);

	if (Laser.beginAbsorption[i][k] > J_FIRST)		//a dielectric cap exists
	{
		nodeA = IJKToIndex(i, Laser.beginAbsorption[i][k] - 1, k);
	
		n2 = LOOK_UP(Phase[A(PhaseCodeSolid)]->indexN, A(T));				//index in dielectric medium

		if (n2 != 1.0)		// a dielectric
		{
			rad = 4.0 * Geometry.yMap[Laser.beginAbsorption[i][k]] * n2 * PI / 
				Laser.waveLength;		//optical length (assumes n2 constant throughout dielectric)

			sqN2 = SQ (n2);
			sqN3 = SQ (n3);
			sqK3 = SQ (k3);

			t = 8 * sqN2 * n3;		//numerator
			t /= sqK3 + sqN2 + sqK3 * sqN2 + SQ (sqN2) + 4.0 * sqN2 * n3 + sqN3 + sqN2 * sqN3 +
				(sqN2 - 1.0) * ((sqK3 - sqN2 + sqN3) * cos (rad) - 2.0 * k3 * n2 * sin (rad));
	
  			return(t);
		}; //endif
	}; //endif

	return((4.0 * n3) / (SQ(k3) + SQ (n3 + 1.0)) );		//simple planar transmission

}; //endfunc


//______________________________________________________
//
// 	Radiate
//	Radiate energy using the Stefan-Boltzmann formula from
//	the first absorbing node in each column
//______________________________________________________
void Radiate()
{
	int i, j, k;
	unsigned int nodeA;

	double emissivity;

  	for (i = I_FIRST; i < I_LAST; i++)
	{
		for (k = K_FIRST; k < K_LAST; k++)
		{
      		j = Laser.beginAbsorption[i][k];
      		nodeA = IJKToIndex(i, j, k);

			emissivity = PropertyAverage(nodeA, Phase[A(PhaseCodeSolid)]->emissivity, 
						Phase[A(PhaseCodeLiquid)]->emissivity);

      		QLaser[i][j][k] -= 5.6697e-12 * emissivity * pow (A(T) - Sim.tEnvironment, 4.0) *
        					Sim.sClock.curDTime * A(AreaXZ);
		}; //endloop k
    }; //endloop i
}; //endfunc

