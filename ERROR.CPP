//  _________________________________________________________
// |
// |   Error.cpp   Error reporting routines
// |
// |
// |   (C) 1998-99  Columbia University, MSME
// |   (C) 2000-01	Harvard University, DEAS
// |__________________________________________________________

#include <stdarg.h>
#include <iostream>
#include <stdio.h>

#include "3dns.h"  //also includes error.h
#include "files.h"

#define EXT_LEVEL
#include "error.h"
#undef EXT_LEVEL

int modeErr;		//error reporting mode
int modeInfo;		//info reporting mode

//_________________________________________
// Private functions
//_________________________________________
void Msg (char *fmt, va_list ap);


// ______________________________________________________________
// ErrorMsg
//
// ______________________________________________________________
void ErrorMsg (char *fmt, ...)
{
	static char strMessage[255];
	
  	va_list ap;
	va_start (ap, fmt);
	
	memset(strMessage, 0, 255);	//ignore old string

	switch (modeErr)
	{
		case (ERR_CTRLC):	//report error then wait for keypress
  			strcpy (strMessage, "\nError: ");
  			strcat (strMessage, fmt);
  			strcat (strMessage, "\n");
  	
			Msg (strMessage, ap);	//print message to stdio
  			va_end (ap);

			Msg("Press CTRL-C or close window to exit\n", 0);

			EndProgram();
			break;

		case (ERR_MATHLINK):
			vsprintf(strMessage, fmt, ap);  //print to a string
  			va_end (ap);
			throw (strMessage);		//throw back to mathlink handlers
			
			break;

		case (ERR_BATCH):
			strcpy (strMessage, "\nError: ");
  			strcat (strMessage, fmt);
  			strcat (strMessage, "\n");
  	
			Msg (strMessage, ap);	//print message to stdio
  			va_end (ap);

	  		exit(EXIT_FAILURE);
			
			break;
		
		case (ERR_NONE):	//just ignore it
			return;
			break;

		default:
			exit (EXIT_FAILURE);
	}; //endswitch
};	//endfunc


// ______________________________________________________________
// InfoMsg
//
// ______________________________________________________________
void InfoMsg (char *fmt,...)
{
	if (modeInfo != INFO_ON)
		return;

	va_list ap;
  	va_start (ap, fmt);

  	Msg (fmt, ap);
  	va_end (ap);
}; //endfunc


// ______________________________________________________________
// InRange
//              Checks range of a list of strings
// ______________________________________________________________
bool InRange (double param, double vMin, double vMax)
{
	return ((param >= vMin) && (param < vMax));
}; //endfunc

bool InRange (int param, int iMin, int iMax)
{
	return ((param >= iMin) && (param < iMax));
}; //endfunc

bool InRange (double *array, int iStart, int iStop, double vMin, double vMax)
{
	int i;
	bool inRange = true;

	if (array == NULL)
    	return (false);

  	for (i = iStart; ((i < iStop) && inRange); i++)
    	inRange = ((array[i] >= vMin) && (array[i] < vMax));

  	return (inRange);
}; //endfunc

bool InRange (int *array, int iStart, int iStop, int vMin, int vMax)
{
	int i;
	bool inRange = true;

  	if (array == NULL)
    	return (false);

    i = iStart;
    if (iStop == -1)       //negative terminated array
    	while (array[i] != -1 && inRange)
        {
        	inRange = ((array[i] >= vMin) && (array[i] < vMax));
            i++;
        }

    else
		for (i = iStart; ((i < iStop) && inRange); i++)
    		inRange = ((array[i] >= vMin) && (array[i] < vMax));

  	return (inRange);
}; //endfunc


// ______________________________________________________________
// ModeSet
//	Sets reporting modes for info and error messages
// ______________________________________________________________
void SetErrMode(int mErr)
{
	modeErr	= mErr;

	return;
}; //endfunc

void SetInfoMode(int mInfo)
{
	modeInfo = mInfo;

	return;
}; //endfunc


// ______________________________________________________________
// Msg
//
// ______________________________________________________________
void Msg (char *fmt, va_list ap)
{
	vprintf (fmt, ap);		  
	
	if (Environment.fpLog != NULL)
	{
		vfprintf(Environment.fpLog, fmt, ap);
		fflush(Environment.fpLog);
	}; //endif

	return;
}; //endfunc


// ______________________________________________________________
// WarningMsg
//
// ______________________________________________________________
void WarningMsg (char *fmt,...)
{
	char *strMessage = new char[255];
	va_list ap;

	va_start (ap, fmt);

	strcpy (strMessage, "\nWarning: ");
	strcat (strMessage, fmt);
	Msg (strMessage, ap);

	va_end (ap);
}; //endfunc
