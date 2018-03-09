#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Fortran calls */
extern void F77_NAME(filmom)(double *T, double *U, int *KMAX, int *K);

extern void F77_NAME(ix3)(double *P3INDX, int *KMAX, int *K, int *N, double *KSTAT, double *S, int *SSIZE, double *A, double *B, double *C, double *T, double *U, double *D, double *E, double *F, int *TEXT, int *error);

extern void F77_NAME(ix3dvs)(double *P3INDX, double *DP3DX, int *KMAX, int *K, int *N, double *KSTAT, double *DKDX, double *S, double *DSA, double *DSB, double *DSC, int *SSIZE, double *A, double *B, double *C, double *T, double *U, double *BCT, double *ACT, double *ABT, double *AAT, double *BBT, double *CCT, double *ABCU, double *AACU, double *AABU, double *BBCU, double *ABBU, double *BCCU, double *ACCU, double *AAAU, double *BBBU, double *CCCU, double *D, double *E, double *F, int *TEXT, int *error);

extern void F77_NAME(moment)(double *X, int *KMAX, int *NMAX, int *K, int *N, double *T, double *U);

extern void F77_NAME(sphcor)(double *DATMAT, int *MAXROW, int *MAXCOL, int *K, int *N, double *XCEN, double *VAR, double *COR, double *EVALUE, double *EV, double *EET, double *TRNMAT, double *INVMAT, int *error, int *DOMXRW, int *ISUPPZ, int *LWORK, double *WORK, int *LIWORK, int *IWORK, int *INFO, double *LAVAL, double *LAVEC);

extern void F77_NAME(trimsu)(double *DATMAT, int *KMAXD, int *NMAXD, int *K, int *N, double *LIMIT, int *ACTION, double *DATCP, int *COPYN, int *error);

static const R_FortranMethodDef FortranEntries[] = {
    {"filmom", (DL_FUNC) &F77_NAME(filmom),  4},
    {"ix3",    (DL_FUNC) &F77_NAME(ix3),    17},
    {"ix3dvs", (DL_FUNC) &F77_NAME(ix3dvs), 38},
    {"moment", (DL_FUNC) &F77_NAME(moment),  7},
    {"sphcor", (DL_FUNC) &F77_NAME(sphcor), 23},
    {"trimsu", (DL_FUNC) &F77_NAME(trimsu), 10},
    {NULL, NULL, 0}
};

void R_init_PP3(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}
