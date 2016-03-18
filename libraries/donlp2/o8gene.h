#ifndef o8gene_h
#define o8gene_h
/* **************************************************************************** */
/*                                  o8gene.h                                    */
/* **************************************************************************** */

#include "o8para.h"

#define X

#include "o8comm.h"
#include "o8fint.h"

#undef X

DOUBLE zero = 0.e0,one = 1.e0,two = 2.e0,three = 3.e0,four = 4.e0,
five = 5.e0,six = 6.e0,seven = 7.e0,
twom2 = .25e0,twop4 = 16.e0,twop11 = 2048.e0,onep3 = 1.3e0,
onep5 = 1.5e0,p2 = .2e0,p4 = .4e0,p5 = .5e0,p7 = .7e0,p8 = .8e0,p9 = .9e0,
c45 = 45.e0,tm12 = 1.e-12,
tm10 = 1.e-10,tm9 = 1.e-9,tm8 = 1.e-8,tm7 = 1.e-7,tm6 = 1.e-6,tm5 = 1.e-5,
tm4 = 1.e-4,tm3 = 1.e-3,tm2 = 1.e-2,tm1 = 1.e-1,tp1 = 1.e1,tp2 = 1.e2,
tp3 = 1.e3,tp4 = 1.e4;
        
static INTEGER  qpterm;
static INTEGER  fcount;

#define NDUALM  NX+NRESM
#define MDUALM  NRESM*2

static INTEGER  ndual,mi,me,iq;
static DOUBLE   np[NDUALM+1],rnorm,rlow,xj[NDUALM+1][NDUALM+1],
                ddual[NDUALM+1],r[NDUALM+1][NDUALM+1],ud[MDUALM+1],
                ud1[MDUALM+1];

static INTEGER  iptr,iqtr,aitr[MDUALM+1];
static DOUBLE   sstr,riitr;

#endif
