//
//
//

#ifndef TEST_UNSORTED
#define TEST_UNSORTED


#ifndef CONSOLE
#ifndef NDEBUG
#define CONSOLE(x)                                                             \
    do {                                                                       \
        std::cout << __func__ << ":" << x << '\n';                             \
    } while (0)
#else
#define CONSOLE(x)
#endif

#define CONSOLE_EVAL(x) CONSOLE(#x << " : " << (x))
#endif

#endif
