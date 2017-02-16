#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(xdgges)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(xdggev)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(xdggsvd)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(xdsygv)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(xzgges)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(xzggev)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(xzggsvd)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(xzhegv)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
    {"xdgges",  (DL_FUNC) &F77_NAME(xdgges),  16},
    {"xdggev",  (DL_FUNC) &F77_NAME(xdggev),  15},
    {"xdggsvd", (DL_FUNC) &F77_NAME(xdggsvd), 24},
    {"xdsygv",  (DL_FUNC) &F77_NAME(xdsygv),  10},
    {"xzgges",  (DL_FUNC) &F77_NAME(xzgges),  16},
    {"xzggev",  (DL_FUNC) &F77_NAME(xzggev),  15},
    {"xzggsvd", (DL_FUNC) &F77_NAME(xzggsvd), 25},
    {"xzhegv",  (DL_FUNC) &F77_NAME(xzhegv),  11},
    {NULL, NULL, 0}
};

void R_init_geigen(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
