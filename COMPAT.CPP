//  _________________________________________________________
// |
// |   Compat.cpp    Routines subject to compatibility issues
// |   				J.Leonard, S.Dubler, M.Crowder
// |
// |  (C) 1998-99  Columbia University, MSME
// |  (C) 2000-01	Harvard University, DEAS
// |__________________________________________________________

#include "compat.h"
#include "error.h"		

#ifdef __unix__             		//GCC / unix
#include <unistd.h>
#define DIR_SEPARATOR "/"
#define CHAR_DIR_SEPARATOR '/'

#else 		       	//Borland C++ v5.0, Visual C++ 6.0
#include <windows.h>
#include <io.h>
#define DIR_SEPARATOR "\\"
#define CHAR_DIR_SEPARATOR '\\'
#endif

#ifdef __EMX__
#define FIXSLASHES
#endif

#ifdef __CYGWIN32__
#define FIXSLASHES
void _fnslashify(char *);
#endif


// ______________________________________________________________
// EndProgram
//
//    Use an elegant method to exit from the program
// ______________________________________________________________
void EndProgram(void)
{
#ifdef __unix__
	while(true);

#else
	WaitMessage();	//WIN32
#endif

	return;
}; //endfunc


// ______________________________________________________________
// FileExists
//
//    Checks for existence of a file
// ______________________________________________________________
bool FileExists(const char *fileName)
{
	return (_access(fileName, 0) == 0);
}; //endfunc


// ______________________________________________________________
// ParseBase
//
//    Get basename
// ______________________________________________________________
char *ParseBase(char *arg)
{
  	char *baseName = new char[strlen(arg) + 1];

  	char *strPtr1=NULL;
  	char *strPtr2=NULL;
  	char *strPtr3=NULL;

#ifdef FIXSLASHES
  	_fnslashify(arg);
#endif

  	strPtr1 = arg;		       //set to first value

  	while((strPtr2 = strchr(strPtr1, CHAR_DIR_SEPARATOR)) != NULL)
    	strPtr1 = strPtr2 + 1;	       //search to final token name

  	strcpy(baseName, strPtr1);

  	if((strPtr2 = strchr(baseName, '.')) != NULL)
    	*(strPtr2) = '\0';		       //drop off extension if necessary

  	strPtr3 = new char[strlen(baseName) + 1];	//make right length

  	strcpy(strPtr3, baseName);

  	delete [] baseName;

  	return(strPtr3);
};


//_________________________________________________
// ParseFileName
//
//  Constructs filename from a list of directories
//      and names.  Will insert separators as necessary
//_________________________________________________
char *ParseFileName(const char *str1, const char *str2)
{
  	char *strPtr1=0;

  	if((str1 == NULL) ||(str2 == NULL))
    	ErrorMsg("Cant parse NULL strings in ParseFileName");

  	if((*(str1 + strlen(str1) - 1) == CHAR_DIR_SEPARATOR) ||
     	(*(str1 + strlen(str1) - 1) == '.') ||
     	(*(str2) == CHAR_DIR_SEPARATOR) ||
     	(*(str2) == '.'))
    {
      	strPtr1 = new char[strlen(str1) + strlen(str2) + 1];
      	strcpy(strPtr1, str1);
      	strcat(strPtr1, str2);
    }

  	else
    {
      	strPtr1 = new char[strlen(str1) + strlen(str2) + 2];
      	strcpy(strPtr1, str1);
      	strcat(strPtr1, DIR_SEPARATOR);
      	strcat(strPtr1, str2);
    };

  	return(strPtr1);
}; //endfunc

char *ParseFileName(const char *str1, const char *str2, const char *str3)    //handle fewer arguments
{
  	return(ParseFileName(str1, ParseFileName(str2, str3)));
}; //endfunc

char *ParseFileName(const char *str1, const char *str2, const char *str3, const char *str4)
{
  	return(ParseFileName(str1, ParseFileName(str2, str3, str4)));
}; //endfunc


// ______________________________________________________________
// ParseOutput
//
//       Get output directory, returns ptr to a new string that
//   contains the output directory with terminating '\\' or '/'
// ______________________________________________________________
char *ParseOutput(char *arg)
{

  	char *strPtr1=0;
  	char *strPtr2=0;
  	char *strPtr3=0;

#ifdef FIXSLASHES
  	_fnslashify(arg);
#endif

  	strPtr1 = arg;
  	strPtr2 = arg;		       //set to first value

  	while((strPtr3 = strchr(strPtr2, CHAR_DIR_SEPARATOR)) != NULL)
    	strPtr2 = strPtr3 + 1;	       //search to final token name

  	strPtr3 = new char[strPtr2 - strPtr1 + 1];
  	strncpy(strPtr3, strPtr1,(strPtr2 - strPtr1));

  	strPtr3[strPtr2 - strPtr1] = '\0';

#ifdef __unix__
  	if(strPtr1 == strPtr2)
    	strcpy(strPtr3, "./");
#endif

  	return(strPtr3);
}; //endfunc


// ______________________________________________________________
// ParseRoot
//
//   Get root directory, returns ptr to a new string that
//   contains the root directory with DIR_SEPARATOR at the
//       end.
// ______________________________________________________________
char *ParseRoot(char *arg)
{
//____________ GCC / unix ___________________
#ifdef __unix__
  	arg = 0;
  	char *rootName = getenv("TWODPC_ROOT");

#ifdef FIXSLASHES
  	_fnslashify(rootName);
#endif

  	if((rootName == NULL) || strlen(rootName) >= NUMCHAR_DIR - 1)
    	ErrorMsg("Can't get root directory specification");

  	char *strPtr2 = new char[strlen(rootName) + 2];

  	strcpy(strPtr2, rootName);
  	strcat(strPtr2, DIR_SEPARATOR);

  	return(strPtr2);

	
//_____ Borland C++ 5.0 / Visual C++ 6.0 ___________________
#else

	static char *strPtr1=0;
	static char *strPtr2=0;

	char *rootName = new char[strlen(arg) + 1];
  	strcpy(rootName, _strupr(arg));     //arg is the full path & exe name

	strPtr1 = strstr(rootName, PROGRAM_NAME);	//find beginning of program name

	if(strPtr1 == (char *)NULL)
    	ErrorMsg("Cant get root directory specification");

  	*(strPtr1) = '\0';		       //chop off program name

	strPtr2 = new char[strlen(rootName) + 1];
  	strcpy(strPtr2, rootName);

	delete rootName;

	return(strPtr2);
#endif



}; //endfunc




//_________________________________________
// _strupr
//_________________________________________
#ifdef __unix__
char *_strupr(char *s)
{
  	char *c = s;
  	while(*s != 0)
    {
      	if(*s >= 97)
			*s = *s - 32;
      	s++;
    }; //endwhile

  	return c;
}; //endfunc
#endif


//_________________________________________
// _stricmp
//_________________________________________
#ifdef __unix__
int _stricmp(char *s1, char *s2)
{
  	while((*s1 != '\0') &&(*s2 != '\0'))
    {
      	if(toupper(*s1) != toupper(*s2))
			return 1;
      	s1++;
      	s2++;
    }; //endwhile

  	if((*s1 != '\0') ||(*s2 != '\0'))
    	return 1;

  	return 0;
}; //endfunc
#endif


//_________________________________________
// _fnslashify
//_________________________________________

#ifdef __CYGWIN32__
void _fnslashify(char *arg)
{
  	for(; *arg != '\0'; arg++)
    	if(*arg == '\\')
      		*arg = '/';
}; //endfunc
#endif


