//  _________________________________________________________
// |
// |   ReadData.h
// |
// |
// |   (C) 1998-99  Columbia University, MSME
// |   (C) 2000-01	Harvard University, DEAS
// |__________________________________________________________
#ifndef _READDATAH
#define _READDATAH
//_________________________________________
// Public functions
//_________________________________________
void LoadEnvironmentFile();
void LoadParametersFile(string &fileName);
void  RestoreTemperature( char *baseName, double resumeTime);

#endif

