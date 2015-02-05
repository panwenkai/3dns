//  _________________________________________________________
// |
// |   Files.h    File handling header
// |
// |   (C) 1998-99  Columbia University, MSME
// |   (C) 2000-01	Harvard University, DEAS
// |__________________________________________________________
#ifndef _FILESH
#define _FILESH

//_________________________________________
// Typedefs
//_________________________________________
typedef struct
{
    int lineNumber;
    string fileName;
    FILE *fPtr;
} FILEINFO;

typedef struct        // RPN operator
{
	int errMode;	//error reporting mode
	int infoMode;	//info reporting mode
	int lineNumber;

	string *searchPath;		//directory search path
	FILE *fpLog;
} ENVIRONMENT;

typedef struct
{
	int numMembers;
	char *member[MAX_MEMBERS];
} STRINGLIST;

//_________________________________________
// Public functions
//_________________________________________
bool FileOpenRead(FILEINFO &fInfo, const string *searchPath=0);
bool FileOpenWrite(FILE ** fptr, char *path);
void FileClose(FILEINFO &fInfo);
void FileClose(FILE *fPtr);
bool ReadLine(FILEINFO &fInfo, ENTRY &entry);
bool ReadMatrix2D(FILEINFO &fInfo, MATRIX2D &newMatx, ENTRY &entry);
bool ReadMatrix2D(MATRIX2D &newMatx, const char *fileName, string *pathSearch=0);
bool ReadMatrix3D(MATRIX3D &newMatx, const char *fileName, string *pathSearch=0);
string SearchThruPath(const string *searchPath, const string fileName);

//_________________________________________________
// Public variables
//_________________________________________________

#ifdef EXT_LEVEL
#define EXTERN
#else
#define EXTERN extern
#endif

#undef EXT_LEVEL
EXTERN ENVIRONMENT Environment;		//environment variables
#undef EXTERN

#endif				  

