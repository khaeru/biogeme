
/* **************************************************************************** */
/*                                  o8fuco.h                                    */
/* **************************************************************************** */

    
X LOGICAL   val[NRESM+1],gconst[NRESM+1],llow[NX+1],lup[NX+1];
X INTEGER   gunit[4][NRESM+1];

X INTEGER   n,nr,nres,nh,ng;

X DOUBLE    epsmac,tolmac,deldif;

X char      name[41];

X DOUBLE    epsdif;

X LOGICAL   intakt,inx,DONLP2std,te0,te1,te2,te3,singul;
X LOGICAL   ident,eqres,silent,analyt,cold;

X INTEGER   icf,icgf,cfincr,cres[NRESM+1],cgres[NRESM+1];

X LOGICAL   ffuerr,cfuerr[NRESM+1];

