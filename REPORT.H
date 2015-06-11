//  _________________________________________________________
// |
// |   Report.h
// |
// |
// |   (C) 1998-99  Columbia University, MSME
// |   (C) 2000-01	Harvard University, DEAS
// |__________________________________________________________
#ifndef _REPORTH
#define _REPORTH

//_________________________________________
// Public functions
//_________________________________________
void LogInit(char *fileName);
void ReportCleanup (void);
void ReportData(void);
void ReportInit (void);
void ReportNucleationEvent(int nodeA, int dirToLiq);
void ReportProperties(char *baseName);

char *const SUFFIX_HETEROGENEOUS = ".het";
char *const SUFFIX_HISTORY = ".hist";
char *const SUFFIX_HOMOGENEOUS = ".hom";
char *const SUFFIX_INTERFACE = ".int";
char *const SUFFIX_LASER = ".las";
char *const SUFFIX_LOG = ".log";
char *const SUFFIX_NUCLEATION = ".nuc";
char *const SUFFIX_PROPERTIES = ".prop";
char *const SUFFIX_TEMPERATURE = ".temp";

#endif
