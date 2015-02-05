//  _________________________________________________________
// |
// |   Compat.h
// |
// |
// |   (C) 1998-99  Columbia University, MSME
// |   (C) 2000-01	Harvard University, DEAS
// |__________________________________________________________
#ifndef _COMPATH
#define _COMPATH

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef __unix__             		//GCC / unix
#define DIR_INCLUDE "library"
#else 		       //for Borland C++ v5.0, Visual C++ 6.0
#define PROGRAM_NAME "3DNS.EXE"
#define ENVIRONMENT_NAME "3DNS.INI"
#define WIN32_LEAN_AND_MEAN
#endif

#ifndef _WINDEF_
#define max(x,y) ((x>y) ? x : y)
#endif

#define DOUBLE_INFINITY 1E+300		//a large finite number < MAXDOUBLE
#define PI 3.14159265325

#ifdef __EMX__
#include <float.h>
#endif

#define NUMCHAR_KEY  ((int)32)	       //keyword name length
#define NUMCHAR_LINE ((int)8192)		//maximum line length for entry
#define NUMCHAR_DIR  ((int)255)		   //maximum directory line length
#define NUMDIR_PATH  ((int)16)		   //maximum number of directories in path

//_________________________________________
// typedefs
//_________________________________________

//_________________________________________
// Public functions
//_________________________________________

#ifdef __linux__
char *_strupr (char *s);
int _stricmp (char *s1, char *s2);
#endif

void EndProgram(void);
bool FileExists(const char *fileName);
char *ParseBase (char *arg);
char *ParseFileName (const char *path, const char *base, const char *ext);
char *ParseFileName (const char *str1, const char *str2);
char *ParseFileName (const char *str1, const char *str2, const char *str3);
char *ParseFileName (const char *str1, const char *str2, const char *str3, const char *str4);
char *ParseOutput (char *arg);
char *ParseRoot (char *arg);
char *ParseDirectory (char *dirDest, char *dirSrc);

#endif