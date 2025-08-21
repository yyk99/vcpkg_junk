#pragma once

#ifndef NDEBUG
#   define CONSOLE(x) do { std::cout << __func__ << ":" << x << '\n';  } while(0)
#   define CONSOLE_THR(x) do { std::cout << __func__ << ":" << std::this_thread::get_id() << ":" << x << '\n';  } while(0)
#else
#   define CONSOLE(x)
#   define CONSOLE_THR( x )
#endif

#define CONSOLE_EVAL(x) CONSOLE(#x << " : " << (x))

#define STATUS(x) do { \
    std::ostringstream ss; \
    ss << x; \
    CONSOLE ( ss.str() ); \
    GM_LOG_MSG_DIRECT(GM_Log_Status, ss.str().c_str()); \
} while(0)


