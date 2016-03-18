#ifndef o8para_h
#define o8para_h

/* **************************************************************************** */
/*                                  o8para.h                                    */
/* **************************************************************************** */
#include <cstdio>
#include <cstdlib>

#include <cmath>

#include "patType.h"

/* Macros */

// WARNING: conflict with the compiler
//#define max(A,B)    ((A) > (B) ? (A) : (B))
//#define min(A,B)    ((A) < (B) ? (A) : (B))

#define myfabs(A) ((A)<0 ? (-(A)) : (A))
/* Types */

typedef int         INTEGER;
typedef int         LOGICAL;
enum                {FALSE,TRUE};
typedef float       REAL;
typedef patReal      DOUBLE;

/* Dimensions */

#define NX          110
#define NRESM       1500
#define MAXIT       4000
#define NSTEP       40

#endif
