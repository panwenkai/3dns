//  _________________________________________________________
// |
// |   ReadData.cpp    Input routines
// |
// |
// |   (C) 1998-99  Columbia University, MSME
// |   (C) 2000-01	Harvard University, DEAS
// |__________________________________________________________

#include "3dns.h"
#include "parsefunc.h"
#include "files.h"
#include "report.h"
#include "phase.h"
#include <string.h>

#define EXT_LEVEL
#include "readdata.h"
#undef EXT_LEVEL

#define DELIM2 "\t ,;}" 		// delimeters for a list of values
#define NO_PATH ((string *)0)

//_________________________________________
// Private functions
//_________________________________________
void ResizeOverlay( );
int  LoadMaterialFile( char *fn, int phase);
void LoadOverlayBlock(FILEINFO &fInfo, ENTRY &entry);
REGION *RegionNew(int i);

//_________________________________________
// Private variables
//_________________________________________

//______________________________________________________
//
//	LoadEnvironmentFile
//		Read environment file(s) in output, then root directory
//______________________________________________________
void LoadEnvironmentFile()
{
	ENTRY entry;
    FILEINFO fInfo;
	string *envPath=0;
	int i, curLine;
	bool foundIni;

	Report.pathSearch = StringListNew(2);
	
	Report.pathSearch[0] = Report.dirOutput;
	Report.pathSearch[1] = Report.dirRoot;

	for (i = 0; i < 2; i++)		//search thru both directories
	{
		fInfo.fileName = ParseFileName(Report.pathSearch[i], ENVIRONMENT_NAME);

		SetErrMode(ERR_NONE);		//quiet search for files
		SetInfoMode(INFO_OFF);

		foundIni = FileOpenRead(fInfo);

		SetErrMode(Sim.modeError);	//restore mode	
		SetInfoMode(INFO_ON);

		if (foundIni)
		{
			curLine = 0;

			while( ReadLine(fInfo, entry) )		//it was opened successfully now get data
			{
				curLine++;

				if( entry.keyWord == "PATH" )
				{
					ParseLine(entry, T_STRING, NUM_ANY, NULL, NO_PATH);	//parse to words
	  				entry.value.Load(&envPath);			//list of directories
 
					Report.pathSearch = StringListCat(Report.pathSearch, envPath);

					StringListDelete(&envPath);		//clean deletion of dynamic array
				} //endif

       			else
	       			ErrorMsg("Unrecognized entry in %s, line %d", fInfo.fileName, curLine);
				
				entry.Zero();
			} //endwhile

			fclose( fInfo.fPtr);
		} //endif
	} //endloop i

	return;
} //endfunc
	

//______________________________________________________
//
// LoadMaterialFile
//		Reads in a series of phases in a file
//		stores the index list in Region[regionNum]
//		returns the total number 
//______________________________________________________
int LoadMaterialFile(string &fileName, int regionNum)
{
	PHASE *thisPhase;
	REGION *thisRegion = Region[regionNum];	

	static int *fileIndex = new int [MAX_MATERIAL_CLASS];
	static string *matlClassList = StringListNew(MAX_MATERIAL_CLASS);
	static string *fileNameList = StringListNew(MAX_MATERIAL_CLASS);
	
	static int curMatlClass=0;	//total number of material classes
	static int curPhase=1;		//total number of phases
	
	string str;
	string matlClassName;
	string *transitionList=0;

	int i, liqPhase;
	int matlClassNum;
	int phasesFound=0;			//number of phases found in this file
	int numTransitions=0;		//number of transitions specified in the liquid phase

	bool inBlock = false;	//inside phase block
	bool inClass = false;	//inside class definition
	bool inLiquid = false;	//inside liquid phase

	ENTRY entry;
  	FILEINFO fInfo;

	//_______________ FILE RESOLUTION ___________________
	i = StringListIndex(&fileNameList, fileName);

	if (i < 0)		//this is a new file
	{
		i = StringListCount(fileNameList);		
		StringListAdd(fileNameList, fileName, MAX_MATERIAL_CLASS);		//add new name to list
		fileIndex[i] = regionNum;

		fInfo.fileName = fileName;
		InfoMsg("Reading MATERIAL file");
		FileOpenRead(fInfo, Report.pathSearch);		//open material file for reading
	}

	else		//already read the file, so use info from that region
	{
		thisRegion->phaseIndex =				//use same list of phases
			Region[ fileIndex[i] ]->phaseIndex;

		thisRegion->phaseLiquid = Region[fileIndex[i]]->phaseLiquid;
												//copy liquid phase # if any
		return(curPhase);
	} //endif

	//______________ READ NEW FILE ___________________
	while( ReadLine( fInfo, entry))		//loop thru every line in file
   	{
      	if( entry.keyWord == "MATERIAL_CLASS_NAME" )
		{
			if (inClass)
				ErrorMsg("MATERIAL_CLASS_NAME cant be re-entered on line %d", fInfo.lineNumber);

			inClass = true;
					
	        ParseLine(entry, T_STRING, 1, NULL, NO_PATH);
		  	entry.value.Load(&matlClassName);
			
			matlClassNum = StringListIndex(&(matlClassList), matlClassName);

			if (matlClassNum < 0)		//no match so new class
			{
				matlClassNum = StringListCount(matlClassList);
				if(!StringListAdd(matlClassList, matlClassName, MAX_MATERIAL_CLASS))
					ErrorMsg("Couldn't initialize new material class");
				curMatlClass++;
			} //endif
			
			if (curMatlClass == MAX_MATERIAL_CLASS)
				ErrorMsg("Too many unique material classes in line %d", fInfo.lineNumber);

			continue;
		} //endif
			
		if (!inClass)
			ErrorMsg("Must define MATERIAL_CLASS_NAME before PHASE in line %d", fInfo.lineNumber);

		if( (entry.keyWord == "PHASE") )
        {
			switch (entry.typeFound)
			{
				case T_BEGIN:
					if (inBlock)
						ErrorMsg("PHASE END is needed in line %d", fInfo.lineNumber);

					else
					{
						inBlock = true;
						
						thisPhase = PhaseNew(curPhase + phasesFound);
						thisPhase->matlClass = matlClassNum + 1;	//class 0 is null class
						thisPhase->matlClassName = matlClassName;
					} //endif

					break;
				
				case T_END:
					if (!inBlock)
						ErrorMsg("PHASE END should not be in line %d", fInfo.lineNumber);
					else
					{
						inBlock = false;
						phasesFound++;					
					} //endif
					
					break;
						
				default:  

					ErrorMsg("Need BEGIN or END in line %d", fInfo.lineNumber);
			} //endswitch
		
			continue;		//get next line
		} //endif

		if (!inBlock)		//need to have initiated a phase
			ErrorMsg("Need PHASE BEGIN in line %d", fInfo.lineNumber);

		if( entry.keyWord == "PHASE_NAME" )
		{
			ParseLine(entry, T_STRING, 1, NULL, NO_PATH);
			entry.value.Load(&(thisPhase->phaseName));

			if (thisPhase->phaseName == "LIQUID")	//a liquid phase
			{
				inLiquid = true;
				thisPhase->phaseType = LIQUID;
			}

			else
			{
				inLiquid = false;
				thisPhase->phaseType = SOLID;
			}
				
			continue;
		} //endif

		if (!(bool)(thisPhase->phaseName))
			ErrorMsg("Must define PHASE NAME first in line %d", fInfo.lineNumber);

		if( entry.keyWord == "EMISSIVITY" )
        {
            ParseLine(entry, T_CONST | T_FUNCTION | T_RPN | T_MATRIX2D, NULL, VAL_ANY, Report.pathSearch);
	  		GenTable(entry.value, &(thisPhase->emissivity), MIN_DEGREES, MAX_DEGREES, 1.0);
        }

		else if( entry.keyWord == "ENTHALPY" )
        {
			ParseLine(entry, T_CONST | T_FUNCTION | T_RPN | T_MATRIX2D, NULL, VAL_ANY, Report.pathSearch);
	  		GenTable(entry.value, &(thisPhase->enthalpyH), MIN_DEGREES, MAX_DEGREES, 1.0);
        }

      	else if( entry.keyWord == "HEAT_CAPACITY" )
        {
            ParseLine(entry, T_CONST | T_FUNCTION | T_RPN | T_MATRIX2D, NULL, VAL_ANY, Report.pathSearch);
	  		GenTable(entry.value, &(thisPhase->heatCapacity), MIN_DEGREES, MAX_DEGREES, 1.0);
        }

      	else if( entry.keyWord == "INTERFACE_RESPONSE" )
	  	{
            ParseLine(entry, T_CONST | T_FUNCTION | T_RPN | T_MATRIX2D, NULL, VAL_ANY, Report.pathSearch);
	  		GenTable(entry.value, &(thisPhase->interfaceResponse), MIN_DEGREES, MAX_DEGREES,
				M_TO_CM);
        }

      	else if( entry.keyWord == "INDEX_K" )
	  	{
            ParseLine(entry, T_CONST | T_FUNCTION | T_RPN | T_MATRIX2D, NULL, VAL_ANY, Report.pathSearch);
	  		GenTable(entry.value, &(thisPhase->indexK), MIN_DEGREES, MAX_DEGREES, 1.0);
        }

      	else if( entry.keyWord == "INDEX_N" )
	  	{
            ParseLine(entry, T_CONST | T_FUNCTION | T_RPN | T_MATRIX2D, NULL, VAL_ANY, Report.pathSearch);
	  		GenTable(entry.value, &(thisPhase->indexN), MIN_DEGREES, MAX_DEGREES, 1.0);
        }

		else if( entry.keyWord == "NUCLEATION_HET" )
	  	{
			ParseLine(entry, T_CONST | T_FUNCTION | T_RPN | T_MATRIX2D, NULL, VAL_ANY, Report.pathSearch);
	  		GenTable(entry.value, &(thisPhase->hetNucleation), MIN_DEGREES, MAX_DEGREES, 1.0);
       	}

      	else if( entry.keyWord == "NUCLEATION_HOM" )
	 	{
            ParseLine(entry, T_CONST | T_FUNCTION | T_RPN | T_MATRIX2D, NULL, VAL_ANY, Report.pathSearch);
	  		GenTable(entry.value, &(thisPhase->homNucleation), MIN_DEGREES, MAX_DEGREES, 1.0);
        }

      	else if( entry.keyWord == "THERMAL_CONDUCTIVITY" )
	  	{
            ParseLine(entry, T_CONST | T_FUNCTION | T_RPN | T_MATRIX2D, NULL, VAL_ANY, Report.pathSearch);
	  		GenTable(entry.value, &(thisPhase->thermConductivity), MIN_DEGREES, MAX_DEGREES, 1.0);
        }

      	else if( entry.keyWord == "TRANSITION_TEMP" )
		{
			if (inLiquid)
				ParseLine(entry, T_INT, &numTransitions, VAL_ANY, NO_PATH);   //can have multiple transitions
			else
				ParseLine(entry, T_INT, 1, VAL_ANY, NO_PATH);	// only single value allowed for solid phases

			entry.value.Load(&(thisPhase->transitionTemp));
        }

		else if( entry.keyWord == "TRANSITION_TO" )
		{
			if (inLiquid)
			{
				ParseLine(entry, T_STRING, &numTransitions, NULL, NO_PATH);   //can have multiple transitions
				entry.value.Load(&(transitionList));		//load in transition list for later resolution
				thisPhase->numTransitions = entry.numElements;
			}

			else
			{
				ParseLine(entry, T_STRING, 1, NULL, NO_PATH);	// only single value allowed for solid phases
				entry.value.Load(&str);

				if (!(str == "LIQUID"))
					ErrorMsg("Solid phase must transition to liquid in line %d", fInfo.lineNumber);
			} //endif
        }

      	else
			ErrorMsg( "Unrecognized token <%s> in line %d", entry.keyWord, fInfo.lineNumber);
    } //endwhile

	if (inBlock)
		ErrorMsg("Missing PHASE END in material file <%s>", fInfo.fileName);

	FileClose(fInfo);

	if (!(phasesFound))
		ErrorMsg("No phases found in material file <%s>", fInfo.fileName);
	
	thisRegion->phaseIndex = new int [phasesFound+1];
	
	for (i = 0; i < phasesFound; i++)				//copy phase index
		thisRegion->phaseIndex[i] = curPhase + i;

	thisRegion->phaseIndex[i] = -1;			//negative terminated list

	//______ RESOLVE TRANSITION LINKS __________
	liqPhase = PhaseFind(Region[regionNum]->phaseIndex, "LIQUID");

	if (liqPhase)		//do liquid phase
	{
		if ((phasesFound - 1) != Phase[liqPhase]->numTransitions)
			ErrorMsg("Bad TRANSITION_TO in line %d, check number of solid phases", fInfo.fileName);
		
		thisRegion->phaseLiquid = liqPhase;					
		Phase[liqPhase]->transitionTo = new int[numTransitions];

		for (i=0; i<numTransitions; i++)  //resolve each transition to
		{
			Phase[liqPhase]->transitionTo[i] = 
				PhaseFind(Region[regionNum]->phaseIndex, transitionList[i]);

			if (Phase[liqPhase]->transitionTo[i] == 0)
				ErrorMsg("Couldnt resolve TRANSITION_TO in liquid phase");
		} //endloop i
	} //endif

	StringListDelete(&transitionList);	

	for (i = 0; i < phasesFound; i++)
	{
		if (Phase[curPhase + i]->phaseType == SOLID)	//a solid phase
		{
			Phase[curPhase + i]->transitionTo = new int [1];
			Phase[curPhase + i]->transitionTo[0] = liqPhase;
		}
	} //endloop i

	curPhase += phasesFound;		//update current phase number
	
	InfoMsg("done\n");
	
	return(curPhase);		//total number of phases loaded
} //endfunc


//______________________________________________________
//
//      LoadOverlayBlock
//
//______________________________________________________
void LoadOverlayBlock(FILEINFO &fInfo, ENTRY &entry, int regionNum)
{
	REGION *thisRegion;

	string startName;

	if (regionNum == (MAX_REGIONS-1) )
		ErrorMsg("Too many overlays or layers");

	if (entry.typeFound != T_BEGIN)
		ErrorMsg("Bad overlay block in line %d, must use 'OVERLAY BEGIN'", fInfo.lineNumber);

	if (LAYER_LAST == 0)
		ErrorMsg("OVERLAY block appears before LAYER definition in line %d", Sim.curLine);

	thisRegion = RegionNew(regionNum);	//allocate new region

  	while( ReadLine(fInfo, entry) )    //read until reach end or start of new block
    {
    	if (entry.keyWord == "OVERLAY")
		{
			if (entry.typeFound == T_END)
				break;
			else
				ErrorMsg("Couldnt find end of overlay block in line %d", fInfo.lineNumber);
		}

      	else if( entry.keyWord == "OVERLAY_LOCATIONS_I" )
        {
            ParseLine(entry, T_INT, NUM_ANY, VAL_NON_NEGATIVE, NO_PATH);
			entry.value.Load(&(thisRegion->iLocations) );
			thisRegion->numberI = entry.numElements;
        }

      	else if( entry.keyWord == "OVERLAY_LOCATIONS_J" )
        {
            ParseLine(entry, T_INT, NUM_ANY, VAL_NON_NEGATIVE, NO_PATH);
			entry.value.Load(&(thisRegion->jLocations) );
			thisRegion->numberJ = entry.numElements;
        }

		else if( entry.keyWord == "OVERLAY_LOCATIONS_K" )
        {
            ParseLine(entry, T_INT, NUM_ANY, VAL_NON_NEGATIVE, NO_PATH);
			entry.value.Load(&(thisRegion->kLocations) );
			thisRegion->numberK = entry.numElements;
        }

      	else if( entry.keyWord == "OVERLAY_CAN_CHANGE" )
        {
            ParseLine(entry, T_BOOL, 1, NULL, NO_PATH);
			entry.value.Load(&(thisRegion->canChange) );
        }

      	else if( entry.keyWord == "OVERLAY_CATALYZE_FREEZING" )
        {
            ParseLine(entry, T_BOOL, 1, NULL, NO_PATH);
			entry.value.Load(&(thisRegion->catalyzeFreezing) );
        }

      	else if( entry.keyWord == "OVERLAY_CATALYZE_MELTING" )
		{
            ParseLine(entry, T_BOOL, 1, NULL, NO_PATH);
			entry.value.Load(&(thisRegion->catalyzeMelting) );
        }

      	else if( entry.keyWord == "OVERLAY_HET_THRESHOLD" )
		{
            ParseLine(entry, T_FLOAT, 1, VAL_NON_NEGATIVE | VAL_FLAG_INFINITY, NO_PATH);
			entry.value.Load(&(thisRegion->thresholdHet) );
        }

        else if( entry.keyWord == "OVERLAY_HOM_THRESHOLD" )
        {
            ParseLine(entry, T_FLOAT, 1, VAL_NON_NEGATIVE | VAL_FLAG_INFINITY, NO_PATH);
			entry.value.Load(&(thisRegion->thresholdHom) );
        }

      	else if( entry.keyWord == "OVERLAY_MATERIAL" )
		{
			string temp;

			ParseLine(entry, T_STRING, 1, NULL, NO_PATH);
	  		entry.value.Load(&temp);

			Sim.numberPhases = LoadMaterialFile(temp, regionNum);
		}
		
		else if( entry.keyWord == "OVERLAY_PHASE_START" )
		{
			ParseLine(entry, T_STRING, 1, NULL, NO_PATH);
	  		entry.value.Load(&(startName));
		}

      	else
			ErrorMsg( "Unrecognized overlay block token <%s> line %d",
		  		entry.keyWord, fInfo.lineNumber);
    } //endwhile

	thisRegion->phaseStart = PhaseFind(thisRegion->phaseIndex, startName);

	if (thisRegion->phaseStart == 0)
		ErrorMsg("Cant resolve starting phase for OVERLAY block line %d", Sim.curLine);

	return;
} //endfunc


//______________________________________________________
//
//   LoadParametersFile
//       Reads main parameter file
//______________________________________________________
void LoadParametersFile(string &datFileName)
{
    int i, j, k, r;

	int curRegion = 0;		//current region being read

	FILEINFO fInfo;
  	ENTRY entry;
	REGION *thisRegion;

	InfoMsg ("\nReading PARAMETER file", datFileName);
	fInfo.fileName = datFileName;
	FileOpenRead(fInfo);		//open file (without need for search path)
	InfoMsg("\n");
  
	Phase = (PHASE **)new void*[MAX_PHASES];
	memset(Phase, 0, sizeof(void *) * MAX_PHASES);

	Region = (REGION **)new void*[MAX_REGIONS];
	memset(Region, 0, sizeof(void *) * MAX_REGIONS);

  	while( ReadLine( fInfo, entry) )		//loop thru every line in file
    {
      	if( entry.keyWord == "ENERGY_DENSITY" )
        {
            ParseLine(entry, T_FLOAT, 1, VAL_NON_NEGATIVE, NO_PATH);
			entry.value.Load(&Laser.energyDensity);
        }

		else if( entry.keyWord == "INTERFACE_MOTION_LO" )
        {
            ParseLine(entry, T_FLOAT, &(Sim.eraMax), VAL_POSITIVE, NO_PATH);
			entry.value.Load(&(Sim.dInterfaceLo));
		}

		else if( entry.keyWord == "INTERFACE_MOTION_HI" )
        {
            ParseLine(entry, T_FLOAT, &(Sim.eraMax), VAL_POSITIVE, NO_PATH);
			entry.value.Load(&(Sim.dInterfaceHi));
		}

      	else if( entry.keyWord == "LASER_SPATIAL_X" )
        {
            ParseLine(entry, T_CONST | T_FUNCTION | T_RPN | T_MATRIX1D | T_MATRIX2D, 1, VAL_ANY, Report.pathSearch);
			Laser.dataSpatialX = entry.value;	//generic-- finish in LaserInit
        }

		else if( entry.keyWord == "LASER_SPATIAL_Z" )
        {
            ParseLine(entry, T_CONST | T_FUNCTION | T_RPN | T_MATRIX1D, 1, VAL_ANY, Report.pathSearch);
			Laser.dataSpatialZ = entry.value;	//generic-- finish in LaserInit
        }


        else if( entry.keyWord == "LASER_TEMPORAL" )
        {
            ParseLine(entry, T_CONST | T_FUNCTION | T_RPN | T_MATRIX1D, 1, VAL_ANY, Report.pathSearch);
			Laser.dataTemporal = entry.value;	//generic-- finish in LaserInit
        }

      	else if( entry.keyWord == "LASER_WAVELENGTH" )
		{
            ParseLine(entry, T_FLOAT, 1, VAL_POSITIVE, NO_PATH);
        	entry.value.Load(&(Laser.waveLength));
			Laser.waveLength *= M_TO_CM;
        }

      	else if( entry.keyWord == "LAYER_CAN_CHANGE" )
        {
			bool *temp;
			ParseLine(entry, T_BOOL, &(Geometry.jZones), NULL, NO_PATH);
			entry.value.Load(&temp);
			
			for (r = LAYER_FIRST; r < LAYER_LAST; r++)	//distribute info to regions
			{
				thisRegion = RegionNew(r);
				thisRegion->canChange = temp[r]; 
			} //endloop 
			
			delete [] temp;
        }

      	else if( entry.keyWord == "LAYER_CATALYZE_FREEZING" )
        {
            bool *temp;
			ParseLine(entry, T_BOOL, &(LAYER_LAST), NULL, NO_PATH);
			entry.value.Load(&temp);

			for (r = LAYER_FIRST; r < LAYER_LAST; r++)	//distribute info to regions
			{
				thisRegion = RegionNew(r);
				thisRegion->catalyzeFreezing = temp[r]; 
			} //endloop 

			delete [] temp;
		}

      	else if( entry.keyWord == "LAYER_CATALYZE_MELTING" )
        {
            bool *temp;
            ParseLine(entry, T_BOOL, &(LAYER_LAST), NULL, NO_PATH);
			entry.value.Load(&temp);

			for (r = LAYER_FIRST; r < LAYER_LAST; r++)	//distribute info to regions
			{
				thisRegion = RegionNew(r);
				thisRegion->catalyzeMelting = temp[r]; 
			} //endloop 
			
			delete [] temp;
        }

        else if( entry.keyWord == "LAYER_HET_THRESHOLD" )
        {
            double *temp;
			
			ParseLine(entry, T_FLOAT, &(Geometry.jZones), VAL_NON_NEGATIVE | VAL_FLAG_INFINITY, NO_PATH);
			entry.value.Load(&temp);

			for (i=LAYER_FIRST; i<LAYER_LAST; i++)	//distribute info to regions
			{
				thisRegion = RegionNew(i);
				thisRegion->thresholdHet = temp[i]; 
			} //endloop 

			delete [] temp;
        }

      	else if( entry.keyWord == "LAYER_HOM_THRESHOLD" )
        {
            double *temp;
			
			ParseLine(entry, T_FLOAT, &(Geometry.jZones), VAL_NON_NEGATIVE | VAL_FLAG_INFINITY, NO_PATH);
			entry.value.Load(&temp);

			for (i=LAYER_FIRST; i<LAYER_LAST; i++)	//distribute info to regions
			{
				thisRegion = RegionNew(i);
				thisRegion->thresholdHom = temp[i]; 
			} //endloop 

			delete [] temp;
        }

      	else if( entry.keyWord == "LAYER_MATERIALS" )
		{
			string *temp=0;

			ParseLine(entry, T_STRING, &(LAYER_LAST), NULL, NO_PATH);
	  		entry.value.Load(&temp);

	  		for(r = LAYER_FIRST; r < LAYER_LAST; r++)
			{
				thisRegion = RegionNew(r);
				Sim.numberPhases = LoadMaterialFile(temp[r], r);
			} //endloop 
			
			StringListDelete(&temp);
		}

		else if( entry.keyWord == "LAYER_PHASE_START" )
		{
			string *temp=0;

			ParseLine(entry, T_STRING, &(Geometry.jZones), NULL, NO_PATH);
	  		entry.value.Load(&(temp));

	  		for(i=LAYER_FIRST; i<LAYER_LAST; i++)
			{
				thisRegion = RegionNew(i);		//new regions if necessary
				thisRegion->phaseStart = PhaseFind(thisRegion->phaseIndex, temp[i]);
				
				if (thisRegion->phaseStart == 0)
					ErrorMsg("Cant resolve starting phases in line %d, check that material files are loaded", Sim.curLine);
			} //endloop 

			StringListDelete(&temp);
		}

		else if( entry.keyWord == "MODE_RADIATION" )
        {
        	ParseLine(entry, T_BOOL, 1, NULL, NO_PATH);
			entry.value.Load(&(Sim.modeRadiation));			
        }

		else if( entry.keyWord == "MODE_RESOLIDIFY_STOP" )
        {
        	ParseLine(entry, T_BOOL, 1, NULL, NO_PATH);
			entry.value.Load(&(Sim.modeResolidify));			
        }

        else if( entry.keyWord == "MODE_STOCHASTIC" )
        {
        	ParseLine(entry, T_BOOL, 1, NULL, NO_PATH);
			entry.value.Load(&(Sim.modeStochastic));			
        }


		else if( entry.keyWord == "NODES_I" )
		{
            ParseLine(entry, T_INT, &(Geometry.iZones), VAL_POSITIVE, NO_PATH);
			entry.value.Load(&(Geometry.nodesI));

	  		Geometry.iBound = new int[Geometry.iZones + 1];
	  		Geometry.iBound[0] = 0;

	  		for( i = 1; i <= Geometry.iZones; i++)
				Geometry.iBound[i] = Geometry.iBound[i - 1] + Geometry.nodesI[i - 1];
		}

      	else if( entry.keyWord == "NODES_J" )
		{
            ParseLine(entry, T_INT, &(Geometry.jZones), VAL_POSITIVE, NO_PATH);
			entry.value.Load(&(Geometry.nodesJ));

	  		Geometry.jBound = new int[Geometry.jZones + 1];
	  		Geometry.jBound[0] = 0;

	  		for( j = 1; j <= Geometry.jZones; j++)
	    		Geometry.jBound[j] = Geometry.jBound[j - 1] + Geometry.nodesJ[j - 1];
		}

      	else if( entry.keyWord == "NODES_K" )
		{
            ParseLine(entry, T_INT, &(Geometry.kZones), VAL_POSITIVE, NO_PATH);
			entry.value.Load(&(Geometry.nodesK));

	  		Geometry.kBound = new int[Geometry.kZones + 1];
	  		Geometry.kBound[0] = 0;

	  		for( k = 1; k <= Geometry.kZones; k++)
	    		Geometry.kBound[k] = Geometry.kBound[k - 1] + Geometry.nodesK[k - 1];
		}

      	else if( entry.keyWord == "OVERLAY" )
        {
			if (Geometry.jZones == 0)
				ErrorMsg("Must define layers before overlays in line %d", Sim.curLine);
			
			LoadOverlayBlock(fInfo, entry, curRegion++);
		}

/*      else if( entry.keyWord == "PERIODIC_IJK" )
        {
            ParseLine(entry, T_BOOL, 3, NULL, NO_PATH);
			Geometry.modePeriodic = new bool[3];
			entry.value.Load(&(Geometry.modePeriodic));
        }
*/

		else if( entry.keyWord == "REPORT_HISTORY" )
        {
            ParseLine(entry, T_BOOL, 1, NULL, NO_PATH);
			entry.value.Load(&(Report.watchHistory));
        }

      	else if( entry.keyWord == "REPORT_INTERFACE" )
        {
            ParseLine(entry, T_BOOL, 1, NULL, NO_PATH);
			entry.value.Load(&(Report.watchInterface));
        }

      	else if( entry.keyWord == "REPORT_INTERVAL" )
		{
        	ParseLine(entry, T_FLOAT, &(Sim.eraMax), VAL_NON_NEGATIVE, NO_PATH);
			entry.value.Load(&(Report.tInterval));
        }

      	else if( entry.keyWord == "REPORT_LASER" )
		{
            ParseLine(entry, T_BOOL, 1, NULL, NO_PATH);
			entry.value.Load(&(Report.watchLaser));
        }

        else if( entry.keyWord == "REPORT_NUCLEATION" )
        {        
            ParseLine(entry, T_BOOL, 1, NULL, NO_PATH);
			entry.value.Load(&(Report.watchNucleation));
        }

		else if( entry.keyWord == "REPORT_PROBABILITY" )
        {
            ParseLine(entry, T_BOOL, 1, NULL, NO_PATH);
			entry.value.Load(&(Report.watchProbability));
        }

		else if( entry.keyWord == "REPORT_PROPERTIES" )
        {
            ParseLine(entry, T_BOOL, 1, NULL, NO_PATH);
			entry.value.Load(&(Report.watchProperties));
        }

      	else if( entry.keyWord == "REPORT_TEMPERATURE" )
        {
            ParseLine(entry, T_BOOL, 1, NULL, NO_PATH);
			entry.value.Load(&(Report.watchTemperature));
		}

		else if( entry.keyWord == "SEED_STOCHASTIC" )
		{	
			ParseLine(entry, T_INT, 1, VAL_POSITIVE | VAL_FLAG_ALL, NO_PATH);
			entry.value.Load(&(Sim.seedStochastic));
		}

      	else if( entry.keyWord == "SIZE_I" )
		{
            ParseLine(entry, T_FLOAT, &(Geometry.iZones), VAL_POSITIVE, NO_PATH);
			entry.value.Load(&(Geometry.sizeI));
        }

        else if( entry.keyWord == "SIZE_J" )
    	{
            ParseLine(entry, T_FLOAT, &(Geometry.jZones), VAL_POSITIVE, NO_PATH);
			entry.value.Load(&(Geometry.sizeJ));			
        }

        else if( entry.keyWord == "SIZE_K" )
    	{
            ParseLine(entry, T_FLOAT, &(Geometry.kZones), VAL_POSITIVE, NO_PATH);
			entry.value.Load(&(Geometry.sizeK));			
        }

		else if( entry.keyWord == "TEMPERATURE_CHANGE_LO" )
        {
            ParseLine(entry, T_FLOAT, &(Sim.eraMax), VAL_POSITIVE, NO_PATH);
			entry.value.Load(&(Sim.dTemperatureLo));
		}

		else if( entry.keyWord == "TEMPERATURE_CHANGE_HI" )
        {
            ParseLine(entry, T_FLOAT, &(Sim.eraMax), VAL_POSITIVE, NO_PATH);
			entry.value.Load(&(Sim.dTemperatureHi));
		}

        else if( entry.keyWord == "TEMP_SUBSTRATE" )
        {
        	ParseLine(entry, (T_INT | T_FLOAT | T_MATRIX3D), 1, VAL_POSITIVE, NO_PATH);
			Sim.tSubstrate = entry.value;
        }

		else if( entry.keyWord == "TIME_STEP_LO" )
        {
            ParseLine(entry, T_FLOAT, &(Sim.eraMax), VAL_POSITIVE, NO_PATH);
			entry.value.Load(&(Sim.dTimeLo));
		}

		else if( entry.keyWord == "TIME_STEP_HI" )
        {
            ParseLine(entry, T_FLOAT, &(Sim.eraMax), VAL_POSITIVE, NO_PATH);
			entry.value.Load(&(Sim.dTimeHi));
		}

      	else if( entry.keyWord == "TIME_ERA" )
		{
        	ParseLine(entry, T_FLOAT, &(Sim.eraMax), VAL_POSITIVE, NO_PATH);
			entry.value.Load(&(Sim.eraDuration));
        }

		else if( entry.keyWord == "WATCH_INTERFACE_I" )
        {
        	ParseLine(entry, T_INT, NUM_ANY, VAL_NON_NEGATIVE | VAL_FLAG_ALL, NO_PATH);
			entry.value.Load(&(Report.interfaceWatch[DIM_I]));
        }

		else if( entry.keyWord == "WATCH_INTERFACE_J" )
        {
        	ParseLine(entry, T_INT, NUM_ANY, VAL_NON_NEGATIVE | VAL_FLAG_ALL, NO_PATH);
			entry.value.Load(&(Report.interfaceWatch[DIM_J]));
        }

		else if( entry.keyWord == "WATCH_INTERFACE_K" )
        {
        	ParseLine(entry, T_INT, NUM_ANY, VAL_NON_NEGATIVE | VAL_FLAG_ALL, NO_PATH);
			entry.value.Load(&(Report.interfaceWatch[DIM_K]));
        }

		else if( entry.keyWord == "WATCH_PHASE_I" )
        {
        	ParseLine(entry, T_INT, NUM_ANY, VAL_NON_NEGATIVE | VAL_FLAG_ALL, NO_PATH);
			entry.value.Load(&(Report.phaseWatch[DIM_I]));
        }

      	else if( entry.keyWord == "WATCH_PHASE_J" )
        {
        	ParseLine(entry, T_INT, NUM_ANY, VAL_NON_NEGATIVE | VAL_FLAG_ALL, NO_PATH);
			entry.value.Load(&(Report.phaseWatch[DIM_J]));
       	}

		else if( entry.keyWord == "WATCH_PHASE_K" )
        {
        	ParseLine(entry, T_INT, NUM_ANY, VAL_NON_NEGATIVE | VAL_FLAG_ALL, NO_PATH);
			entry.value.Load(&(Report.phaseWatch[DIM_K]));
       	}

      	else if( entry.keyWord == "WATCH_TEMPERATURE_I" )
		{
        	ParseLine(entry, T_INT, NUM_ANY, VAL_NON_NEGATIVE | VAL_FLAG_ALL, NO_PATH);
			entry.value.Load(&(Report.tempWatch[DIM_I]));
        }

      	else if( entry.keyWord == "WATCH_TEMPERATURE_J" )
		{
        	ParseLine(entry, T_INT, NUM_ANY, VAL_NON_NEGATIVE | VAL_FLAG_ALL, NO_PATH);
			entry.value.Load(&(Report.tempWatch[DIM_J]));
		}

      	else if( entry.keyWord == "WATCH_TEMPERATURE_K" )
		{
        	ParseLine(entry, T_INT, NUM_ANY, VAL_NON_NEGATIVE | VAL_FLAG_ALL, NO_PATH);
			entry.value.Load(&(Report.tempWatch[DIM_K]));
		}
		/*
      	else
			ErrorMsg( "Unrecognized parameter <%s> line %d in %s",
		  		entry.keyWord, fInfo.lineNumber, fInfo.fileName);
		*/
		//we comment this out to enable parameters used only in the TR simulation but not
		//3dns heat folow simulation - 2014.11.07 - JW
		curRegion = (curRegion == 0) ? Geometry.jZones : curRegion;
	} //endwhile

  	FileClose(fInfo);

	Geometry.numberRegions = curRegion;		//total number of layers plus overlays
} //endfunc


//______________________________________________________
//
//	RegionNew
//      Allocates a new region or connects to existing
//______________________________________________________
REGION *RegionNew(int i)
{
	if (Region[i] == NULL)
	{
		Region[i] = new REGION;
		memset(Region[i], 0, sizeof(REGION));
	} //endif

	return(Region[i]);
} //endfunc


