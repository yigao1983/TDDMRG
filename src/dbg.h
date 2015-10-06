#ifndef dbg_h
#define dbg_h

#include <stdio.h>
#include <string.h>
#include <errno.h>

#ifdef NDBG
#define debug(M, ...)
#else
#define debug(M, ...) fprintf(stderr, "[DEBUG] (%s:%d)" M "\n", __FILE__, __LINE__, ##__VAR_ARGS__)
#endif

#define freeup(A) { if(A) free(A); A = NULL; }

#define error_msg() (errno == 0 ? "None" : strerror(errno))

#define print_error(M, ...) fprintf(stderr, "[ERROR] (%s:%d: errno: %s) " M "\n", __FILE__, __LINE__, error_msg(), ##__VA_ARGS__)

#define check(A, M, ...) if(!(A)) { print_error(M, ##__VA_ARGS__); errno = 0; goto error; }

#define check_mem(A, M) if(!(A)) { check(A, "Memory error: " M); errno = 0; goto error; }

#define sentinel(M, ...) { print_error(M, ##__VA_ARGS__); errno = 0; goto error; }

#endif
