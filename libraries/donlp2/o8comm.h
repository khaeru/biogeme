/* **************************************************************************** */
/*                                  o8comm.h                                    */
/* **************************************************************************** */
                                      
#include "o8fuco.h"

X REAL      runtim;
X REAL      optite;
X DOUBLE    accinf[MAXIT+1][33];
X INTEGER   itstep,phase;

X DOUBLE    upsi,upsi0,upsi1,upsist,psi,psi0,
            psi1,psist,psimin,
            phi,phi0,phi1,phimin,fx,fx0,fx1,
            fxst,DONLPfmin,b2n,b2n0,xnorm,x0norm,sig0,dscal,dnorm,d0norm;
X DOUBLE    sig,sigmin,dirder,cosphi,upsim;
X DOUBLE    x[NX+1],x0[NX+1],x1[NX+1],xmin[NX+1],d[NX+1],d0[NX+1],
            dd[NX+1],difx[NX+1],resmin[NRESM+1];

X DOUBLE    gradf[NX+1],gfn,qgf[NX+1],gphi0[NX+1],gphi1[NX+1],
            gres[NX+1][NRESM+1],gresn[NRESM+1];

X INTEGER   perm[NX+1],perm1[NX+1],colno[2*NRESM+1],donlp2rank;
X DOUBLE    qr[NX+1][NRESM+1],betaq[NRESM+1],diag[NRESM+1],cscal[NRESM+1],
            colle[NRESM+1];

/* colno also used o8qpso with double length ! */

X DOUBLE    a[NX+1][NX+1],scalm,scalm2,diag0[NX+1],matsc;

X INTEGER   violis[NSTEP*NRESM+1],DONLP2alist[NRESM+1],bind[NRESM+1],
            bind0[NRESM+1],DONLP2sort[NRESM+1];
                        
 X DOUBLE    u[NRESM+1],u0[NRESM+1],o8w[NRESM+1],w1[NRESM+1],res[NRESM+1], 
             res0[NRESM+1],res1[NRESM+1],resst[NRESM+1],scf,scf0, 
             yu[NRESM+1],slack[NRESM+1],infeas,work[NRESM+1]; 

X INTEGER   iterma,ifill1;
X DOUBLE    del,del0,del01,delmin,tau0,tau,ny;
X DOUBLE    smalld,smallw,rho,rho1,eta,epsx,c1d,
            scfmax,updmy0,tauqp,taufac,taumax;

X DOUBLE    alpha,beta,theta,sigsm,sigla,delta,stptrm;
X DOUBLE    delta1,stmaxl;

X DOUBLE    level;
X INTEGER   clow,lastdw,lastup,lastch;

X FILE      *prou,*meu;

X DOUBLE    ug[NX+1],og[NX+1],delfac[NRESM+1];

X INTEGER   nreset;

X DOUBLE    xst[NX+1];

