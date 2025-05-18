//
// Created on 2023/01/06
//

#ifndef DEBUGGINGCONSOLE_H
#define DEBUGGINGCONSOLE_H

// DEBUGGING
#if _WIN32

#define WIN32_LEAN_AND_MEAN
#include <windows.h>

inline void attachDebugConsole()
{
#if _DEBUG
    if ( ::GetStdHandle( STD_OUTPUT_HANDLE ) == NULL ) {
        // these next few lines create and attach a console
        // to this process.  note that each process is only allowed one console.
        AllocConsole();
        AttachConsole( GetCurrentProcessId() );
        (void)freopen( "CON", "w", stdout );
        (void)freopen( "CON", "w", stderr );
        (void)freopen( "CON", "r+", stdin );

        // Prevent killing the whole app by closing the console 
        SetConsoleCtrlHandler( NULL, true );
        HWND hwnd = ::GetConsoleWindow();
        if ( hwnd != NULL )
        {
            HMENU hMenu = ::GetSystemMenu( hwnd, FALSE );
            if ( hMenu != NULL ) DeleteMenu( hMenu, SC_CLOSE, MF_BYCOMMAND );
        }

        SetConsoleCP( CP_UTF8 );
        SetConsoleOutputCP( CP_UTF8 );

        printf( "DEBUGGING console output. Actual CP-%d, Console CP-%d\n", GetACP(), GetConsoleOutputCP() );
    }
#endif
}

#endif

#include "CONSOLE.h"

#endif // DEBUGGINGCONSOLE_H
