//  _________________________________________________________
// |
// |   Files.cpp    File handling functions
// |
// | (C) 1998-99  Columbia University, MSME
// | (C) 2000-01	Harvard University, DEAS
// |__________________________________________________________

#include "3dns.h"
#include "parsefunc.h"

#define EXT_LEVEL
#include "files.h"
#undef EXT_LEVEL

//_________________________________________
// Private functions
//_________________________________________

//______________________________________________________
//
//  FileClose
//    Clean close of a file
//______________________________________________________
void FileClose(FILEINFO &fInfo)
{
  	fclose(fInfo.fPtr);
	return;
}; //endfunc

void FileClose(FILE *fPtr)
{
  	fclose(fPtr);
	return;
}; //endfunc


//______________________________________________________
//
//  FileOpenRead
//    Searches for and opens a file 
//______________________________________________________
bool FileOpenRead(FILEINFO &fInfo, const string *searchPath)
{
	string tempFileName;

	tempFileName = fInfo.fileName;		
	memset(&fInfo, 0, sizeof(FILEINFO));	//clear out junk

  	if (searchPath != NULL)		//try users path
		fInfo.fileName = SearchThruPath(searchPath, tempFileName);
	else		//no path !
		fInfo.fileName = tempFileName;

   	if (!(bool)fInfo.fileName)	//couldnt find it, quietly return
		ErrorMsg("File %s not found-- check search path", tempFileName);			

	fInfo.fPtr = fopen( fInfo.fileName, "r");

	if (fInfo.fPtr == NULL)  //no file found
	{
		ErrorMsg( "Could not open %s.", (char *)fInfo.fileName);
		return(false);
	}; //endif

	InfoMsg(" %s...", fInfo.fileName);

  	fInfo.lineNumber = 0;

	return(true);	//successfully opened it
}; //endfunc


// ______________________________________________________________
// FileOpenWrite
//		Routine to open a file given mode and path
// ______________________________________________________________
bool FileOpenWrite(FILE **fptr, char *path)
{
  	if((*(fptr) = fopen(path, "w")) != NULL)
		return(true);

	else
		ErrorMsg("Can't open file %s.\n", path);

	return(false);
}; //endfunc

	
//______________________________________________________
//
//  ReadLine
//    Retrieve the next well structured line in the parameter file and
//    break up the line and store it in string.
//______________________________________________________
bool ReadLine(FILEINFO &fInfo, ENTRY &entry)
{
	int i=0,j=0;

	static char newLine[NUMCHAR_LINE];		//raw string data
	static string newList[MAX_COLUMNS];	//parsed into words

	char *rtnVal=0;
	char *linePtr=0;

	if (fInfo.fPtr == NULL)
		return(false);

	do
    {
      	rtnVal = fgets( newLine, NUMCHAR_LINE, fInfo.fPtr);

        if(rtnVal == NULL)		//end of file or error
			return(false);

      	fInfo.lineNumber++;

    } while( IsBlank(newLine) || IsComment(newLine));  //enddo

	StringListZero(newList);	//zero out return string list
//    Sim.curLine = fInfo.lineNumber;

	EraseComment(newLine);					//drop off ending comments
	linePtr = newLine;
	Environment.lineNumber = fInfo.lineNumber;	//use for error reporting

	while (ParseWord(&linePtr, newList[i++]))	//parse into separate strings
	{
		if (i == MAX_COLUMNS)
			ErrorMsg("Too many columns at line number %d in %s", fInfo.lineNumber, fInfo.fileName);
	}; //endwhile

	if (i - 1)		//if some data
	{
		ParseEntry(entry, newList);		//copy into entry structure
		return(true);
	}		

	else
		return(false);
}; //endfunc


//______________________________________________________
//
//  ReadMatrix2D
//    Reads in a matrix of 2D data, but stops at the end
//    of file or next keyword
//______________________________________________________
bool ReadMatrix2D(FILEINFO &fInfo, MATRIX2D &newMatx, ENTRY &entry)
{
	FILE *tempPtr = fInfo.fPtr;	//remember previous location
	MATRIX2D tempMatx;
	double *sizeH=0;
	
	static double tempData[MAX_2D];

	int i;
	int rowNum=0;
	int maxColumns, maxRows;
	int numColumns=0;
	
	entry.Zero();		//clear out existing data

	while (ReadLine(fInfo, entry))
	{
		if (entry.keyWord == "PAGE")	   //another page
			break;					
		
		ParseLine(entry, T_RAW, &numColumns, VAL_ANY, (string *)0);  //if non-float data

		if (!rowNum)   //first row, so initialize
		{
			maxColumns = entry.numElements;
			maxRows = MAX_2D / maxColumns;  
			
			if (maxRows > MAX_ROWS)
				maxRows = MAX_ROWS;

//			tempMatx.Reset(maxRows, maxColumns, 0.0);
		}; //endif

		if (rowNum >= maxRows-1)
		{
			InfoMsg("Too many rows of data at line %d", fInfo.lineNumber);
			break;		//just keep whats been read
		}; //endif
		
		entry.value.Load(&sizeH);
		
		if (sizeH)			//some data was obtained
		{
			for (i=0; i<entry.numElements; i++)
				*(tempData + (maxRows * i) + rowNum) = *(sizeH + i);	//transpose into array

			delete [] sizeH;
		};
		
		entry.Zero();	//dump everything
		
		rowNum++;
	}; //endwhile

	if (rowNum)
	{
		newMatx.Reset(maxColumns, rowNum);
		
		for (i=0; i<maxColumns; i++)
			newMatx.Overlay(tempData + (i * maxRows), rowNum, i, 0);
	
		return(true);
	}; //endif

	return(false);

}; //endfunc

//__________ reads a single file _______________
bool ReadMatrix2D(MATRIX2D &newMatx, const char *fileName, string *pathSearch)
{
	bool b;
	FILEINFO fInfo;
	ENTRY entry;

	fInfo.fileName = fileName;
	InfoMsg("  Reading DATA file");
	
	if (!FileOpenRead(fInfo, pathSearch) )		//open data file
		return(false);		

	b = ReadMatrix2D(fInfo, newMatx, entry);
	
	FileClose(fInfo);
	InfoMsg("done\n");

	return(b);

}; //endfunc


//______________________________________________________
//
//  ReadMatrix3D
//    Reads in a matrix of 3D data
//______________________________________________________
bool ReadMatrix3D(MATRIX3D &newMatx, const char *fileName, string *pathSearch)
{
	int curMap = 0;
	int curPage = 0;

	double *val;

	FILEINFO fInfo;
	ENTRY entry;

	fInfo.fileName = fileName;
	InfoMsg("  Reading DATA file");
	
	if (!FileOpenRead(fInfo, pathSearch) )		//open data file
		return(false);		

	newMatx.Reset(MAX_PAGES);

	while (ReadLine(fInfo, entry))
	{		
		while (entry.keyWord == "PAGE")
		{
			if (curPage == MAX_PAGES)
				ErrorMsg("Too many pages in line %d", fInfo.lineNumber);

			ParseLine(entry, T_FLOAT | T_INT, NUM_ANY, VAL_ANY, (string *)0);  //get page ID (if any)
			entry.value.Load(&val);	

			if ( !(newMatx.mapArray[0].Valid()) )
			{
				newMatx.mapArray[0].Reset(1, MAX_PAGES);
				newMatx.mapArray[0].Name() = "PAGEMAP";
			}; //endif

			(newMatx.mapArray[0])[0][curPage] = 
				(entry.numElements > 0) ? *val : (double) curPage;

			ReadMatrix2D(fInfo, newMatx.Page(curPage++), entry);
		}; //endwhile

		if ((bool)entry.keyWord && strcmp(entry.keyWord,"PAGE") != 0)		//a map entry
		{
			if (++curMap == MAX_MAPS)
				ErrorMsg("To many mapping arrays at line %d", fInfo.lineNumber);

			ParseLine(entry, T_FLOAT | T_INT, NUM_ANY, VAL_ANY, (string *)0);  //get floats

			if (entry.numElements > 0)	//some data 
			{
				entry.value.Load(&val);		//load into temporary array

				newMatx.mapArray[curMap].Reset(1, entry.numElements);
				newMatx.mapArray[curMap].Overlay(val, entry.numElements, M_DOMAIN, 0);  //store in matrix
				newMatx.mapArray[curMap].Name() = entry.keyWord;

				delete [] val;		//dump temporary array
			}; //endif
		}; //endif
	}; //endwhile
	
	FileClose(fInfo);
	
	if (curPage)
	{
		newMatx = newMatx.SubMatrix(0, curPage);	//resize 3d matrix
		newMatx.mapArray[0] = newMatx.mapArray[0].SubMatrix(0, 0, 1, curPage);
	}; //endif


	InfoMsg("done\n");

	return(true);

}; //endfunc


// ______________________________________________________________
// SearchPath
//
//    	Search through path list for the location of a file
//		will return full path if found
// ______________________________________________________________
string SearchThruPath(const string *searchPath, const string fileName)
{
	int i=0;
    string rtnParse;

    while(*(searchPath + i) )
    {
    	rtnParse = ParseFileName(*(searchPath + i++), fileName); 

        if (FileExists(rtnParse))
        	return(rtnParse);
    }; //endwhile

    rtnParse.Zero();	//set to null string
	return(rtnParse);
}; //endfunc
