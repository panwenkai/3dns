//  _____________________________________________________
// |
// |    Server3DNS.cpp          Remote server for simulations
// |
// |    Watches the content of a directory to see if 
// |	simulations are available to be run.  Processes
// |	the files as necessary.
// |
// |      v1.0    Scott Dubler                  09-04-98
// |	  v3.0	  J.Leonard / Visual C++6		07-28-99
// |	  v3.5	   added tray icon				08-01-00
// |
// |	Tray Icon courtesy of: 
// |		Tray42 - System tray example
// |		Michael T. Smith / R.A.M. Technology.Copyright 1997 
// |
// |    (C) 1998-2000      Columbia University, MSME
// |_____________________________________________________

#include "Server3DNS.h"
#include <windows.h>  // main Windows header
#include <windowsx.h> // message cracker macros, etc.
#include <shellapi.h> // Shell32 api
#include <process.h>
#include <stdio.h>

//_________________________________________________
// Private Variables
//_________________________________________________
LRESULT CALLBACK wndProc(HWND hwnd, UINT message, WPARAM wParam, LPARAM lParam);
HINSTANCE ghInst;
int gi = 0; // index into ToolTip table -- not thread-safe, but who cares
#define UWM_SYSTRAY (WM_USER + 1) // Sent to us by the systray
#define NUM_TOOLTIPS 10           // Number of ToolTips in our table
char *tttable[NUM_TOOLTIPS] =
{ "Would you like a small drink?", "No, thank you.", "I wasn't talking to you.",
  "Then who were you talking to?", "Your mom.", "Oh. How is she?",
  "She's good.", "I'm glad to hear.", "I'm sure you already knew.", "OOOHHHH!"};

#define ICON_ON 1
#define ICON_OFF 0

//_________________________________________________
// Private functions
//_________________________________________________
BOOL LaunchProcess (HINSTANCE servInstance, LPCTSTR sExe, LPTSTR sCmdLine, 
					LPTSTR curJobName, BOOL bWait /*=FALSE*/, UINT nTimeout /*=INFINITE*/);

void TrayIcon(HINSTANCE servInstance, UINT command, LPSTR curJobName=0);

//_________________________________________________
// WinMain
//_________________________________________________
int APIENTRY WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance,
						LPSTR lpCmdLine, int nCmdShow)
{
	const int FATAL = 0x001;
	const int MINUTE = 60000;	//milliseconds

	const void* HANDLE_NULL = (void *)-1;

	bool fatalErr = false;
	int cOffset;
	char i;

	LPWIN32_FIND_DATA fileFound = new _WIN32_FIND_DATAA;
	HANDLE rtnHandle;	//return handle from function

//	LPTSTR nameSimFile;		//simulation file name

	LPTSTR cmdLine = new char[MAX_PATH];
	LPTSTR dirOutput = new char[MAX_PATH];
	LPTSTR dirSpool = new char[MAX_PATH];
	LPTSTR execFile = new char[MAX_PATH];
	LPTSTR simSpool = new char[MAX_PATH];
	LPTSTR simOutput = new char[MAX_PATH];
    LPTSTR wildCardFile = new char[MAX_PATH];

	LPTSTR dirRoot = GetCommandLine();
	*(strrchr(dirRoot, '\\')) = '\0';		//remove execname
	*(strrchr(dirRoot, '\\')) = '\0';		//remove .. directory

	lstrcpy(dirOutput, dirRoot+1);		//output directory
	lstrcat(dirOutput, "\\Server\\Output\\");

	lstrcpy(dirSpool, dirRoot+1);			//spool directory
	lstrcat(dirSpool, "\\Server\\Spool\\");

	lstrcpy(execFile, dirRoot+1);			//spool directory
	lstrcat(execFile, "\\3dns.exe");

	lstrcpy(wildCardFile, dirSpool);	//wild card search file
	lstrcat(wildCardFile, "*.sim");		
 
	try
	{
		while (true)	//loop through all files
		{
			rtnHandle = FindFirstFile(wildCardFile, fileFound);

			if (rtnHandle != HANDLE_NULL)
			{
				lstrcpy(simSpool, dirSpool);			//source filename
				lstrcat(simSpool, fileFound->cFileName);

				lstrcpy(simOutput, dirOutput);
				lstrcat(simOutput, fileFound->cFileName);

				if (!CopyFile(simSpool, simOutput, true))	//fail if file exists
				{
					cOffset = lstrlen(simOutput)+1;		//end of default string name

					i = 0;
					lstrcat(simOutput, ".0");			//create new extension

					while (i <= 9)
					{
						
						if (!CopyFile(simSpool, simOutput, true))
						{
							i++;		//try another extension
							*(simOutput + cOffset)+= 1;
						}

						else	//it copied okay
						{
							DeleteFile(simSpool);	//get it out of way in spool
							break;				
						}
					};

					if (i > 9)			//some other problem !
						Sleep(5 * MINUTE);	//wait a while to see what happens

					continue;			//try while loop again
				};

				lstrcpy(cmdLine, execFile);
				lstrcat(cmdLine, " -b ");		//runs in batch mode
				lstrcat(cmdLine, simOutput);	//build command line

				if (!LaunchProcess(hInstance, execFile, cmdLine, 
					fileFound->cFileName, true, INFINITE))
					throw(FATAL);

				DeleteFile(simSpool);	//done, so dump file from spool
			} //endif

			else	//nothing in spool so wait
			{
				Sleep(1 * MINUTE);
			};
		
		} //end infinte while loop
	} //endtry

	catch (int errCode)
	{
		return (errCode);		//died
	}

	return(0);

}; //endfunc

//_________________________________________________
// LaunchProcess
//_________________________________________________
BOOL LaunchProcess (HINSTANCE servInstance, LPCTSTR sExe, LPTSTR sCmdLine, 
					LPSTR curJobName, BOOL bWait /*=FALSE*/, UINT nTimeout /*=INFINITE*/)
{
    PROCESS_INFORMATION processInfo;
    STARTUPINFO         startupInfo;

    ZeroMemory (&startupInfo, sizeof (STARTUPINFO));
    startupInfo.cb = sizeof (STARTUPINFO);

	TrayIcon(servInstance, ICON_ON, curJobName);	//put icon in tray

    BOOL bSuccess = CreateProcess (sExe, sCmdLine, NULL, NULL, TRUE, 
		IDLE_PRIORITY_CLASS | CREATE_NO_WINDOW, NULL,
						NULL, &startupInfo, &processInfo);
    if (bWait)
        WaitForSingleObject (processInfo.hProcess, nTimeout);
 
	TrayIcon(servInstance, ICON_OFF);	//remove icon from tray

    CloseHandle (processInfo.hThread);
    CloseHandle (processInfo.hProcess);

    return bSuccess;
}; //endfunc


//_________________________________________________
// TrayIcon
//_________________________________________________
void TrayIcon(HINSTANCE servInstance, UINT command, LPSTR curJobName)
{
	static WNDCLASSEX wc;		//dummy window description class			
	static HWND hWnd=0;			//handle to dummy window
	static NOTIFYICONDATA nid;	//class to handle icon actions
	
	char *classname = "Noise42.NOTIFYICONDATA.hWnd";

	switch (command)
	{
	case ICON_ON:
		wc.cbSize = sizeof(WNDCLASSEX);		//set up dummy window with icon
		wc.style = 0;
		wc.lpfnWndProc = wndProc;		//pointer to call back function
		wc.cbClsExtra = wc.cbWndExtra = 0;
		wc.hInstance = servInstance;
		wc.hIcon = LoadIcon(servInstance, MAKEINTRESOURCE(IDI_42));
		wc.hCursor = LoadCursor(NULL, IDC_ARROW);
		wc.hbrBackground = (HBRUSH)(COLOR_WINDOW + 1);
		wc.lpszMenuName = NULL;
		wc.lpszClassName = classname;
		wc.hIconSm = (HICON)LoadImage(servInstance, MAKEINTRESOURCE(IDI_42), IMAGE_ICON,
                        GetSystemMetrics(SM_CXSMICON),
                        GetSystemMetrics(SM_CYSMICON), 0);
	
		RegisterClassEx(&wc);
		hWnd = CreateWindowEx(0, classname, classname, WS_POPUP, CW_USEDEFAULT, 0,
			CW_USEDEFAULT, 0, NULL, NULL, servInstance, NULL);

		nid.cbSize = sizeof(NOTIFYICONDATA); // size
		nid.hWnd = hWnd; // window to receive notifications
		nid.uID = 1;     // application-defined ID for icon (can be any UINT value)
		nid.uFlags = NIF_MESSAGE |  // nid.uCallbackMessage is valid, use it
                NIF_ICON |     // nid.hIcon is valid, use it
                NIF_TIP;       // nid.szTip is valid, use it
		nid.uCallbackMessage = UWM_SYSTRAY; // message sent to nid.hWnd
		nid.hIcon = (HICON)LoadImage(servInstance, MAKEINTRESOURCE(IDI_42), IMAGE_ICON,
                        GetSystemMetrics(SM_CXSMICON),
                        GetSystemMetrics(SM_CYSMICON), 0); // 16x16 icon
		
		strcpy(nid.szTip, curJobName);

		// NIM_ADD: Add icon; NIM_DELETE: Remove icon; NIM_MODIFY: modify icon
		Shell_NotifyIcon(NIM_ADD, &nid); // This adds the icon
		break;

	case ICON_OFF:
		nid.cbSize = sizeof(NOTIFYICONDATA);
		nid.hWnd = hWnd;
		nid.uID = 1;
		nid.uFlags = NIF_TIP; // not really sure what to put here, but it works
		Shell_NotifyIcon(NIM_DELETE, &nid); // This removes the icon
		
		DestroyWindow(hWnd);		//destroy the dummy window
		break;
	} //endswitch

	return;
}; //endfunc


//this really isnt activated yet ??
LRESULT CALLBACK wndProc(HWND hwnd, UINT message, WPARAM wParam, LPARAM lParam)
{
  POINT pt;
  HMENU hmenu, hpopup;
  NOTIFYICONDATA nid;

  switch (message) {
    case WM_CREATE:
      SetTimer(hwnd, 0x29a, 2000, NULL); // 2-second timer
      return TRUE;

    case WM_TIMER:
      // Every 2 seconds, we will change something -- in this case, the ToolTip
      if (gi >= NUM_TOOLTIPS - 1) gi = 0;
      else gi++;
      strcpy(nid.szTip, tttable[gi]);
      nid.cbSize = sizeof(NOTIFYICONDATA);
      nid.hWnd = hwnd;
      nid.uID = 1;
      nid.uFlags = NIF_TIP; // Only nid.szTip and nid.uID are valid, change tip
      Shell_NotifyIcon(NIM_MODIFY, &nid); // Modify tooltip
      return TRUE;

    case WM_DESTROY:
      nid.cbSize = sizeof(NOTIFYICONDATA);
      nid.hWnd = hwnd;
      nid.uID = 1;
      nid.uFlags = NIF_TIP; // not really sure what to put here, but it works
      Shell_NotifyIcon(NIM_DELETE, &nid); // This removes the icon
      PostQuitMessage(0);
      KillTimer(hwnd, 0x29a);
      return TRUE;

    case UWM_SYSTRAY: // We are being notified of mouse activity over the icon
      switch (lParam) {
        case WM_RBUTTONUP: // Let's track a popup menu
          GetCursorPos(&pt);
          hmenu = LoadMenu(ghInst, MAKEINTRESOURCE(IDM_CONTEXTMAIN));
          hpopup = GetSubMenu(hmenu, 0);

          /* SetForegroundWindow and the ensuing null PostMessage is a
             workaround for a Windows 95 bug (see MSKB article Q135788,
             http://www.microsoft.com/kb/articles/q135/7/88.htm, I think).
             In typical Microsoft style this bug is listed as "by design".
             SetForegroundWindow also causes our MessageBox to pop up in front
             of any other application's windows. */
          SetForegroundWindow(hwnd);
          /* We specifiy TPM_RETURNCMD, so TrackPopupMenu returns the menu
             selection instead of returning immediately and our getting a
             WM_COMMAND with the selection. You don't have to do it this way.
          */
          switch (TrackPopupMenu(hpopup,            // Popup menu to track
                                 TPM_RETURNCMD |    // Return menu code
                                   TPM_RIGHTBUTTON, // Track right mouse button?
                                 pt.x, pt.y,        // screen coordinates
                                 0,                 // reserved
                                 hwnd,              // owner
                                 NULL))             // LPRECT user can click in
                                                    // without dismissing menu
          {
            case IDM_EXIT: DestroyWindow(hwnd); break;
            case IDM_MESSAGEBOX:
              MessageBox(hwnd, "Heh heh heh...", "Death", MB_TASKMODAL);
              break;
          }
          PostMessage(hwnd, 0, 0, 0); // see above
          DestroyMenu(hmenu); // Delete loaded menu and reclaim its resources
          break;

        case WM_LBUTTONDBLCLK:
          SetForegroundWindow(hwnd); // Our MessageBox pops up in front
          MessageBox(hwnd, "You're cool.", "Hey", MB_TASKMODAL);
          break;
        // Other mouse messages: WM_MOUSEMOVE, WM_MBUTTONDOWN, etc.
      }
      return TRUE; // I don't think that it matters what you return.
  }
  return DefWindowProc(hwnd, message, wParam, lParam);
}