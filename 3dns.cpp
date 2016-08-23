//  _____________________________________________________
// |
// |    3dns.cpp          3D numerical simulation
// |
// |    Simulates the heat transfer and phase
// |    transformation behaviour of a thin silicon film
// |    under pulsed laser irradiation.
// |
// |      v1.0    Vikas Gupta                    1995
// |      v2.0    J.Leonard  (WIN95/BC++5.0)     2-10-98
// |      v2.1    J.Leonard / S.Dubler           7-16-98
// |      v2.2    J.Leonard / S.Dubler           7-24-98
// |      v2.3    J.Leonard / S.Dubler           8-01-98
// |      v2.4    J.Leonard / S.Dubler           8-07-98
// |	  v2.5	  J.Leonard						11-20-98
// | 	  v2.6	  J.Leonard 					12-15-98
// | 	  v2.7	  J.Leonard						01-02-99
// |	  v2.8	  J.Leonard						02-15-99
// |	  v3.0	  J.Leonard / Visual C++6, 3D	07-28-99
// |	  v3.7    J.Leonard						11-14-01
// |
// |    (C) 1995-99      Columbia University, MSME
// |    (C) 2000-01	Harvard University, DEAS
// |_____________________________________________________

#define EXT_LEVEL
#include "3dns.h"
#undef EXT_LEVEL

#include "energy.h"
#include "nucleation.h"
#include "parsefunc.h"		       //for GenTable
#include "phase.h"
#include "report.h"
#include "readdata.h"
#include "thermal.h"
#include "files.h"   //temporary
#include <iostream>
#include <io.h>
//_________________________________________________
// private functions
//_________________________________________________
void CleanUp (void);
void GeometryInit (void);
void TimeStepAdjust();

//_________________________________________________
// main
//_________________________________________________
int main (int argc, char **argv)
{
	string datFileName;
    time_t t;

	ZeroStructures ();		      //set all globals structs to known state (0)

	std::cout<< "argc = " << argc << std::endl;
	std::cout<< "argv[0] = " << argv[0] << std::endl;
	std::cout<< "argv[1] = " << argv[1] << std::endl;

	if ((argc == 1) || (argc > 3))
  	{
    	ErrorMsg ("Usage %s [-b] data_file\n", argv[0]);
  	}

  	else if ((argc == 3) && (!strcmp(argv[1], "-b")))
  	{
    	datFileName = argv[2];
    	Sim.modeError = ERR_BATCH;
  	}

  	else if (argc == 2)
    {
    	datFileName = argv[1];
    	Sim.modeError = ERR_CTRLC;
  	}

	else if (argc == 3)
  	{
    	ErrorMsg ("Usage %s [-b] data_file\n", argv[0]);
  	}

	SetErrMode(Sim.modeError);
	SetInfoMode(INFO_ON);		//info messages are reported

	Report.dirRoot = ParseRoot(argv[0]);			//directory of executable
	Report.dirOutput = ParseOutput(datFileName);	//directory of input (.dat) file
	Report.baseName = ParseBase(datFileName);		//base name

	LoadEnvironmentFile();			//load the initialization file
  	LogInit(Report.baseName);	//initialize log file

	InfoMsg("   __________________________________________________ \n");
	InfoMsg("  /                                                 /|\n"); 
	InfoMsg(" /_________________________________________________/ |\n");
	InfoMsg(" |  3DNS  v%s.%s Rapid Solidification Simulator | |\n", VERSION, BUILD);
	InfoMsg(" |                                                 | |\n");
	InfoMsg(" |  J.P. Leonard, A.B. Limanov, James S. Im        | |\n");
	InfoMsg(" |  Materials Science & Engineering                | |\n");
	InfoMsg(" |                                                 | |\n");
	InfoMsg(" |  (C) 1994-2000, Columbia University             | /\n");
	InfoMsg(" |_________________________________________________|/ \n");

  	LoadParametersFile (datFileName);
  	InfoMsg ("...parameters done\n");

	InfoMsg ("Initializing geometry... ");
  	GeometryInit ();
  	InfoMsg ("geometry done\n");

    ReportProperties(Report.baseName);   // dump lookup tables if necessary

	InfoMsg ("Initializing phase-change cellular automata matrices... ");
  	InterfaceInit ();
	PhaseInit ();
  	InfoMsg ("phase done\n");

  	InfoMsg ("Initializing heat-flow finite differences matrices... ");
  	ThermalInit ();
  	InfoMsg ("thermal done\n");

  	InfoMsg ("Initializing laser spatial and temporal profiles... ");
  	LaserInit ();
  	InfoMsg ("laser done\n");

  	InfoMsg ("Initializing nucleation monte-carlo matrices... ");
  	NucleationInit ();
  	InfoMsg ("nucleation done\n");

 	InfoMsg ("Opening report files... ");
  	ReportInit ();
  	InfoMsg ("report done\n");

    time(&t);			       //get current time
	InfoMsg("Begin simulation at: %s\n", ctime(&t));


	//___________ SIMULATION ___________________
  	try
	{
		do		
		{
			ReportData();				//Output data if necessary
      	
			try
			{
				LaserInput();	//laser absorption -> QLaser
				IntraNode();	//internal interface motion -> QInterface
 				HeatFlow();		//heat flow
			}

			catch (int)			//change was too high
			{
				IntraNodeDump();
				HeatFlowDump();
			
				if (!Sim.sClock.Decrement() )
					ErrorMsg("Could not reverse simulation clock");

				TimeStepAdjust();		

				continue;
			}

			IntraNodeFinalize();	//finalize changes to interface
			InterNode();			//inter-node motion, melting, nucleation
			HeatFlowFinalize();		//finalize temperature changes Tnew->T

			TimeStepAdjust();
			
			Report.prevStep = Sim.sClock.curDTime;

		} while ( Sim.sClock.Increment(Report.sClock) ); 
		ReportData();			   //force final report
	
	} //endtry

	catch (int)		//termination of program
	{
		InfoMsg("\n");
	} //endcatch

    time(&t);			       //get current time
	InfoMsg("End simulation at: %s\n", ctime(&t));

  	CleanUp ();			       //free all arrays
	system("PAUSE");
  	return 0;
}


// ______________________________________________________________
// TimeStepAdjust
//     Adjust time step if necessary
// ______________________________________________________________
void TimeStepAdjust()
{
	static int numMax=0;    //number of times exceeded a max
	static int numCalm=0;   //number of times calm condition exceeded
	
	int i, j, k;

	double facIntLo, facIntHi;		//factors relative to interface motion limits
	double facTempLo, facTempHi;	//factors relative to temperature change limits
	double facTimeLo;			    //factors relative to time step
	double factor;

	bool autoStep = (Sim.dTimeLo[THIS_ERA] < Sim.dTimeHi[THIS_ERA]);

	//______ CHECK FOR UPPER BOUND VIOLATIONS __________
	if (Sim.dImax > Sim.dInterfaceHi[THIS_ERA])	//above interface motion limit
	{
		numMax++;
		factor = Sim.dImax / Sim.dInterfaceHi[THIS_ERA];
		
		if(!Sim.sClock.StepAdjust(0.8 / factor) || (numMax >= MAX_RETRIES))
		{
			IndexToIJK(Sim.dInode, i, j, k);
			ErrorMsg("Could not repair interface instability at dI[%d][%d][%d] at time %g", 
					i, j, k, Sim.sClock.curTime);
		} //endif
		
		return;
	}

	else if (Sim.dTmax > Sim.dTemperatureHi[THIS_ERA]) //above temperature change limit
	{
		numMax++;
		factor = Sim.dTmax / Sim.dTemperatureHi[THIS_ERA];

		if(!Sim.sClock.StepAdjust(0.8 / factor) || (numMax >= MAX_RETRIES))
		{
			IndexToIJK(Sim.dTnode, i, j, k);
			ErrorMsg("Could not repair temperature instability dT[%d][%d][%d] at time %g", 
					i, j, k, Sim.sClock.curTime);
		} //endif
		
		return;
	}

	numMax = 0;

	if (Sim.sClock == Report.sClock)  //report clocks are not standard so dont adjust
		return;

	facTimeLo = Sim.sClock.curDTime / Sim.dTimeLo[THIS_ERA];
	facTempLo = Sim.dTmax / Sim.dTemperatureLo[THIS_ERA];
	facIntLo = Sim.dImax / Sim.dInterfaceLo[THIS_ERA];

	facTempHi = Sim.dTmax / Sim.dTemperatureHi[THIS_ERA];
	facIntHi = Sim.dImax / Sim.dInterfaceHi[THIS_ERA];

	//______ RECOVER FROM PREVIOUS SUBCLOCKS ___________
	if ((facTimeLo < 1.0) && (facTimeLo > 0.0) && 
		(facTimeLo > facTempHi) && (facTimeLo > facIntHi) && !autoStep)
	{
//InfoMsg("recovering\n");
		if (!Sim.sClock.StepAdjust(CTRL_PROP / facTimeLo))
			ErrorMsg("Could not adjust clock step");

		return;
	}

	//______ AUTOMATIC ACCELERATION ____________________
	if ((facTimeLo > 0.0) && autoStep)
	{
//InfoMsg("Adjusting dT=%g, dI=%g\n", Sim.dTmax, Sim.dImax);
		if ((facTempLo > 0.0) && (facIntLo > 0.0))  //temperature and interface
			factor = exp(0.333 * log(facTempLo * facIntLo * facTimeLo));

		else if ((facTempLo > 0.0) || (facIntLo > 0.0))	   //one or none
			factor = exp(0.50 * log((facTempLo + facIntLo) * facTimeLo));

		else //no change at all this clock
			factor = .01;

		factor = (factor < 0.01) ? 0.01 : (factor > 100) ? 100 : factor; //chop
		factor = (factor > 1.0) ? CTRL_PROP * factor : factor / CTRL_PROP;

		if ( ((facTempHi / factor) < 1.0) && ((facIntHi / factor) < 1.0) )
			if (!Sim.sClock.StepAdjust(1.0 / factor) )
				ErrorMsg("Could not adjust clock step");
	} //endif

	return;

} //endfunc
	

// ______________________________________________________________
// GeometryInit
//   Initialize simulation globals and geometry
// ______________________________________________________________
void GeometryInit (void)
{
	int i, j, k, r;
	int ir, jr, kr;
	unsigned int nodeA;

	double xPos, yPos, zPos;

	char *typeMsg;
	int typeVal;

	if (!InRange(Geometry.nodesI, 0, Geometry.iZones, 1, MAX_I))
    	ErrorMsg("NODES_I out of range, or not specified");

    if (!InRange(Geometry.nodesJ, 0, Geometry.jZones, 1, MAX_J))
    	ErrorMsg("NODES_J out of range, or not specified");
	
	if (!InRange(Geometry.nodesK, 0, Geometry.kZones, 1, MAX_K))
    	ErrorMsg("NODES_K out of range, or not specified");

	if (!InRange(I_LAST, 1, MAX_I))
		ErrorMsg("Must have at least 1 and less than %d I-nodes", MAX_I);

	if (!InRange(J_LAST, 1, MAX_J))
		ErrorMsg("Must have at least 1 and less than %d J-nodes", MAX_J);

	if (!InRange(K_LAST, 1, MAX_K))
		ErrorMsg("Must have at least 1 and less than %d K-nodes", MAX_K);
	
	Geometry.isThick = new bool[3];
	Geometry.isThick[DIM_I] = (I_LAST >= 1);	//boolean flags
	Geometry.isThick[DIM_J] = (J_LAST >= 1);
	Geometry.isThick[DIM_K] = (K_LAST >= 1);
												//packed integer
	Geometry.dimThick = (int)Geometry.isThick[DIM_I] * THICK_I +
						(int)Geometry.isThick[DIM_J] * THICK_J +
						(int)Geometry.isThick[DIM_K] * THICK_K;

	if ( (I_LAST * J_LAST * K_LAST) > MAX_NODES )
		ErrorMsg("Total number of nodes in simulation must be less than %d", MAX_NODES);

	if (!(Geometry.isThick[DIM_I] || Geometry.isThick[DIM_J] || Geometry.isThick[DIM_K]))
		ErrorMsg("Cant do 0-dimensional simulations, check NODES_I, NODES_J, NODES_K");

    if (!InRange(Geometry.sizeI, 0, Geometry.iZones, MIN_NODESIZE * CM_TO_M, MAX_NODESIZE * CM_TO_M))
    	ErrorMsg("SIZE_I out of range, or not specified");

    if (!InRange(Geometry.sizeJ, 0, Geometry.jZones, MIN_NODESIZE * CM_TO_M, MAX_NODESIZE * CM_TO_M))
    	ErrorMsg("SIZE_J out of range, or not specified");

    if (!InRange(Geometry.sizeK, 0, Geometry.kZones, MIN_NODESIZE * CM_TO_M, MAX_NODESIZE * CM_TO_M))
    	ErrorMsg("SIZE_K out of range, or not specified");

	// compute delX and delY from the data read in
  	Geometry.delX = new double[Geometry.iZones];
	Geometry.delY = new double[Geometry.jZones];
	Geometry.delZ = new double[Geometry.kZones];
  	
	for( ir = 0; ir < Geometry.iZones; ir++)	//compute width of nodes
    {
		Geometry.delX[ir] = Geometry.sizeI[ir] * 
			M_TO_CM /( double) Geometry.nodesI[ir];
        
		if (Geometry.delX[ir] < 2*MIN_NODESIZE)
        	InfoMsg("\nWarning: Node width less than %9.5Gm in I_REGION %d",
            	2 * MIN_NODESIZE * CM_TO_M, ir);
	} //endloop ir
		  	
	for( jr = 0; jr < Geometry.jZones; jr++)	//compute height of nodes
	{
    	Geometry.delY[jr] = Geometry.sizeJ[jr] * 
			M_TO_CM /( double) Geometry.nodesJ[jr];
			
		if (Geometry.delY[jr] < 2*MIN_NODESIZE)
        	InfoMsg("\nWarning: Node height less than %9.5Gm in J_REGION %d",
				2 * MIN_NODESIZE *CM_TO_M, jr);

		Region[jr]->iLocations = new int [I_LAST];	
		Region[jr]->jLocations = new int [Geometry.nodesJ[jr]];
		Region[jr]->kLocations = new int [K_LAST];

		Region[jr]->numberI = I_LAST;
		Region[jr]->numberJ = Geometry.jBound[jr+1] - Geometry.jBound[jr];
		Region[jr]->numberK = K_LAST;

		for (i = I_FIRST; i < I_LAST; i++) //explicit list of i-nodes
			Region[jr]->iLocations[i] = i;

		for (j = Geometry.jBound[jr]; j < Geometry.jBound[jr+1]; j++)
			Region[jr]->jLocations[j - Geometry.jBound[jr]] = j;	//explicit list of j-nodes this layer

		for (k = K_FIRST; k < K_LAST; k++)
			Region[jr]->kLocations[k] = k;

	} //endloop jr

	for( kr = 0; kr < Geometry.kZones; kr++)	//comput depth of nodes
	{
    	Geometry.delZ[kr] = Geometry.sizeK[kr] * 
			M_TO_CM /( double) Geometry.nodesK[kr];
			
		if (Geometry.delZ[kr] < 2*MIN_NODESIZE)
        	InfoMsg("\nWarning: Node height less than %9.5Gm in K_REGION %d",
				2 * MIN_NODESIZE * CM_TO_M, kr);
	} //endloop kr

	MemoryInit ();		//allocate IMATRIX[] and DMATRIX[]

	MatrixNew (&DelX);	//width of each node
  	MatrixNew (&DelY);	//height of each node
	MatrixNew (&DelZ);	//depth of each node
  	MatrixNew (&DelXY);	//diagonal distance XY
	MatrixNew (&DelXZ);	//diagonal distance XZ
	MatrixNew (&DelYZ);	//diagonal distance YZ
	MatrixNew (&DelXYZ); //diagonal distance XYZ
	
	MatrixNew (&AreaXY);	//area in xy plane
	MatrixNew (&AreaXZ);	//area in xz plane
	MatrixNew (&AreaYZ);	//area in yz plane
	MatrixNew (&Volume);	//node volume

	Geometry.xMap = new double[I_LAST];
	Geometry.yMap = new double[J_LAST];
	Geometry.zMap = new double[K_LAST];

	for (i = I_FIRST, ir = 0; i < I_LAST; i++)
	{
		ir += (i == Geometry.iBound[ir+1]);

		for (j = J_FIRST, jr = 0; j < J_LAST; j++)
		{
			jr += (j == Geometry.jBound[jr+1]);

			for (k = K_FIRST, kr = 0; k < K_LAST; k++)
			{
				kr += (k == Geometry.kBound[kr+1]);

	   	 		nodeA = IJKToIndex(i,j,k);

				A(DelX) = Geometry.delX[ir];
	    		A(DelY) = Geometry.delY[jr];
 				A(DelZ) = Geometry.delZ[kr];

	    		A(DelXY) = sqrt(A(DelX)*A(DelX) + A(DelY)*A(DelY));
				A(DelXZ) = sqrt(A(DelX)*A(DelX) + A(DelZ)*A(DelZ));
				A(DelYZ) = sqrt(A(DelY)*A(DelY) + A(DelZ)*A(DelZ));
				A(DelXYZ) = sqrt(A(DelX)*A(DelX) + A(DelY)*A(DelY) + A(DelZ)*A(DelZ));

				A(AreaXY) = A(DelX)*A(DelY);
				A(AreaXZ) = A(DelX)*A(DelZ);
				A(AreaYZ) = A(DelY)*A(DelZ);
				A(Volume) = A(DelX)*A(DelY)*A(DelZ);
	  		} //endloop k
		} //endloop j
	} //endloop i

	for (i = I_FIRST, xPos=0; i < I_LAST; i++)	//x positions of nodes
	{
		Geometry.xMap[i] = xPos;
		xPos += DelX[i][0][0];
	} //endloop i

	for (j = J_FIRST, yPos=0; j < J_LAST; j++)	//y positions of nodes
	{
		Geometry.yMap[j] = yPos;
		yPos += DelY[0][j][0];
	} //endloop j

	for (k = K_FIRST, zPos=0; k < K_LAST; k++)	//z positions of nodes
	{
		Geometry.zMap[k] = zPos;
		zPos += DelZ[0][0][k];
	} //endloop k
	
	for (r = 0; r < Geometry.numberRegions; r++)
	{
		typeMsg = (r < Geometry.jZones) ? ":LAYER" : "OVERLAY";
		typeVal = (r < Geometry.jZones) ? r : r - Geometry.jZones;

		if (!InRange(Region[r]->numberI, 1, I_LAST+1))
			ErrorMsg("%s %d must have between 1 and %d i-nodes", typeMsg, typeVal, I_LAST);

		if (!InRange(Region[r]->numberJ, 1, J_LAST+1))
			ErrorMsg("%s %d must have between 1 and %d j-nodes", typeMsg, typeVal, J_LAST);

		if (!InRange(Region[r]->numberK, 1, K_LAST+1))
			ErrorMsg("%s %d must have between 1 and %d k-nodes", typeMsg, typeVal, K_LAST);

		if (!InRange(Region[r]->iLocations, 0, Region[r]->numberI, I_FIRST, I_LAST))
			ErrorMsg("Bad i-node numbers in %s %d", typeMsg, typeVal);

		if (!InRange(Region[r]->jLocations, 0, Region[r]->numberJ, J_FIRST, J_LAST))
			ErrorMsg("Bad j-node numbers in %s %d", typeMsg, typeVal);

		if (!InRange(Region[r]->kLocations, 0, Region[r]->numberK, K_FIRST, K_LAST))
			ErrorMsg("Bad k-node numbers in %s, %d", typeMsg, typeVal);
	} //endloop r
			
	if (!InRange( Sim.eraMax, 1, MAX_ERAS))
		ErrorMsg("Must have between 1 and %d simulation eras", MAX_ERAS);
	
	if (!InRange( Sim.eraDuration, 0, Sim.eraMax, MIN_CLOCK, MAX_TIME) )
		ErrorMsg("TIME_ERA entries must be between %g and %g seconds", MIN_CLOCK, MAX_TIME);

   	if (!InRange( Sim.dTimeLo, 0, Sim.eraMax, MIN_CLOCK, MAX_CLOCK))
    	ErrorMsg("TIME_STEP_LO is out of range, or not specified");	

   	if (!InRange( Sim.dTimeHi, 0, Sim.eraMax, MIN_CLOCK, MAX_CLOCK))
    	ErrorMsg("TIME_STEP_HI is out of range, or not specified");	

	if (!Sim.sClock.Reset(Sim.eraMax, Sim.eraDuration, 
			Sim.dTimeLo, Sim.dTimeHi))
		ErrorMsg("Cant initialize clock, check time parameters");

}	//endfunc


// ______________________________________________________________
// CleanUp
//              Clean up arrays and files upon termination of program
// ______________________________________________________________
void CleanUp (void)
{
  delete[]Geometry.iBound;
  delete[]Geometry.jBound;
  delete[]Geometry.delX;
  delete[]Geometry.delY;

  MatrixFree (T);

  ReportCleanup ();
  exit (0);
}


//___________________________________________________
// ZeroStructures
//    Set all global struct elements to 0 
//___________________________________________________
void ZeroStructures ()
{
	memset (&Environment, 0, sizeof (Environment));
	memset (&Geometry, 0, sizeof (Geometry));
	memset (&Laser, 0, sizeof (Laser));
	memset (&Report, 0, sizeof (Report));
	memset (&Sim, 0, sizeof (Sim));
} //endfunc

