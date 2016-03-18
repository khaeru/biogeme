#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "patMath.h"
#include "o8gene.h"
#include "patDisplay.h"
#include "patFileExists.h"
#include "patIterationBackup.h"
#include "donlp2Parameters.h"

#include <cstring>

/* Conditions of use:                                                           */

/* 1. donlp2 is under the exclusive copyright of P. Spellucci                   */
/*    (e-mail:spellucci@mathematik.tu-darmstadt.de)                             */
/*    "donlp2" is a reserved name                                               */
/* 2. donlp2 and its constituent parts come with no warranty, whether ex-       */
/*    pressed or implied, that it is free of errors or suitable for any         */
/*    specific purpose.                                                         */
/*    It must not be used to solve any problem, whose incorrect solution        */
/*    could result in injury to a person , institution or property.             */
/*    It is at the users own risk to use donlp2 or parts of it and the          */
/*    author disclaims all liability for such use.                              */
/* 3. donlp2 is distributed "as is". In particular, no maintenance, support     */
/*    or trouble-shooting or subsequent upgrade is implied.                     */
/* 4. The use of donlp2 must be acknowledged, in any publication which contains */
/*    results obtained with it or parts of it. Citation of the authors name     */
/*    and netlib-source is suitable.                                            */
/* 5. The free use of donlp2 and parts of it is restricted for research purposes*/
/*    commercial uses require permission and licensing from P. Spellucci.       */

/* d o n l p  2                                                                 */

/*    version 29/05/98 (*)                                                      */
/*    tauqp dependent on scalres only                                           */
/*    weights computed in a modified version in the singular case               */
/*    some comparisons relative to zero changed                                 */
/*    error-return for function evaluation as added feature                     */
/*    termination of QP-solver changed (d not set to zero)                      */
/*    new version of BFGS: if nr  =  0 take Powell's update                     */
/*    no suppression of update beforehand (with exception dg = 0)               */
/*    plus some minor corrections                                               */

/*    for consistency reasons variable names cf and cgf changed into            */
/*    icf and icgf                                                              */

/*    added feature numerical differentiation of order 1,2,6                    */
/*    requires new parameters epsfcn, taubnd, difftype                          */

/*    added feature of external blockwise evaluation of functions               */
/*    rather than individual evaluation                                         */
/*    requires new parameter bloc                                               */

/*    new feature of user-supplied scaling of x                                 */

/* **************************************************************************** */
/* (*) direct translation of -, from f77 to ANSI C language                     */
/*     by S. Schoeffert, ASB, Bourges, France                                   */

patString donlp2(donlp2Parameters theParameters, patIterationBackup* theBackup) {

    void    o8st   (void);
    void    o8fin  (void);
    patString    o8opti (donlp2Parameters theParameters, patIterationBackup* theBackup);
    DOUBLE  o8vecn (INTEGER nl,INTEGER nm,DOUBLE x[]);
    void    esgradf(DOUBLE x[],DOUBLE gradf[]);
    void    esgradh(INTEGER i,DOUBLE x[],DOUBLE gradhi[]);
    void    esgradg(INTEGER i,DOUBLE x[],DOUBLE gradgi[]);
    
    void    setup0 (void);
    void    setup  (donlp2Parameters theParameters);
    
    static DOUBLE   yy[NX+1];
    static INTEGER  i,j;
    
    /* default settings of new parameters */
    
    bloc     = FALSE;
    analyt   = TRUE;
    valid    = FALSE;
    epsfcn   = 1.e-16;
    difftype = 3;
    taubnd   = 1.;
    for (i = 1 ; i <= NX ; i++) {
        xsc[i] = one;
        xtr[i] = zero;
    }
    epsdif = tm8;

    for (i = 0 ; i <= NRESM ; i++) {
    
        /* presently, the test-cases in the testenvironment do not differentiate*/
        /* between linear and                                                   */
        /* nonlinear constraints. for later extensions or users functions this: */

        gconst[i] = FALSE;
        val[i]    = FALSE;
        if ( i > 0 ) gresn[i] = one;
    }
    for (i = 1 ; i <= NRESM ; i++) {
        cfuerr[i] = FALSE;
        
        /* the error indicator of the function evaluation */
    }
    ffuerr = FALSE;

    /* some standard initialization which may be overwritten by */
    /* setup0 or setup                                          */

    inx    = FALSE;
    silent = FALSE;
    
    /* the interactive input feature is no longer supported here. for the sake */
    /* of easy revision the variable is supplied here however                  */
    /* if intakt is set TRUE, output to protocol file is copied to stdout in   */
    /* addition                                                                */
    
    intakt = FALSE;
    te0    = FALSE;
    te1    = FALSE;
    te2    = FALSE;
    te3    = FALSE;
    cold   = TRUE;
    
    /* setup0 must initialize analyt, epsdif, del0, tau0 , n , nh , ng , gunit */
    /* gconst , epsfcn , taubnd , analyt , bloc , difftype                     */
    /* and the initial value for x                                             */
    /* may also change settings of all variables initialized above             */

    setup0();

    for (i = 1 ; i <= n ; i++) {
        xst[i] = x[i];
        if ( xsc[i] == zero ) {
            fprintf(stderr,"scaling variable %i is zero\n",i);
            exit(1);
        } 
        x[i] = x[i]/xsc[i];
    }
    nreset = n;

    o8st();
    
    /* setup may change standard settings of parameters  */
    /* and add some computations in the user environment */

    for (i = 1 ; i <= n ; i++) {
    
        /* ug and og have been evaluted for the original variables */
        /* here we use them internally for the scaled ones         */
        
        if ( llow[i] ) ug[i] = ug[i]/xsc[i];
        if ( lup[i]  ) og[i] = og[i]/xsc[i];
    }
    setup(theParameters);
    
    /* preevaluation of gradients of linear functions */
    /* done only once                                 */

    for (i = 0 ; i <= nres ; i++) {
        if ( gunit[1][i] != 1 && gconst[i] ) {
        
            /* evaluate gradient once */
            
            if ( i == 0 ) {
                val[0] = TRUE;
                
                esgradf(x,gradf);
                
            } else {
                val[i] = TRUE;
                if ( i <= nh ) {
                
                    esgradh(i,x,yy);
                    
                } else {
                
                    esgradg(i-nh,x,yy);
                }
                for (j = 1 ; j <= n ; j++) {
                    gres[j][i] = yy[j];
                }
                gresn[i] = patMax(one,o8vecn(1,n,yy));
            }
        }
    }
    runtim = clock();
    
    /* call the optimizer */

    patString status = o8opti(theParameters,theBackup);
    
    runtim = (clock()-runtim)/CLOCKS_PER_SEC;

    /* do final solution check and output */
    
    o8fin();
    
    return(status);
}

/* **************************************************************************** */
/*        initialization program , standard parameter settings done here        */
/* **************************************************************************** */
void o8st(void) {

    void o8msg    (INTEGER num);
    void esh      (INTEGER i,DOUBLE x[],DOUBLE *hxi);
    void esg      (INTEGER i,DOUBLE x[],DOUBLE *gxi);
    void user_eval(DOUBLE xvar[],INTEGER mode);

    static INTEGER  i,j,k,iz;
    static DOUBLE   tol1 ,xnull[NX+1],bd0,infiny,gxi,hxi,term;
    static char     fil[13],xxx[9] = "xxxxxxxx",name1[9];
    static time_t   tim;
    
    epsmac = pow(two,-20);
    
    L100:

    epsmac = epsmac/two;
    term   = one+epsmac;
    
    if ( term != one ) goto L100;
    
    epsmac = epsmac+epsmac;
    tolmac = epsmac;
    
    L200:

    tol1   = tolmac;
    tolmac = tolmac/twop4;
    
    if ( tolmac != zero ) goto L200;
    
    tolmac = tol1;
    
    /* epsmac machine precision, tolmac smallest machine number */
    /* larger than zero (approximately , base 16 for exponent   */
    /* therefore division by 16 assumed)                        */

    /*                        ***** warning *****                        */
    /* on some machines the computation of tolmac may result in an error */
    /* because underflow is not accepted as zero as is assumed here      */

    if ( n > NX ) {
        fprintf(stderr,"donlp2: n too large/recompile\n");
        exit(1);
    }
    if ( nh+ng > NRESM ) {
        fprintf(stderr,"donlp2:nh or ng too large/recompile\n");
        exit(1);
    }
    if ( tau0 == zero ) tau0 = one;
    if ( del0 == zero ) del0 = tau0*p5;
    
    /* if ( del0 > tau0 ) del0 = tau0*p5 */
    
    if ( nreset > n  ) nreset = n;
    if ( nreset <= 4 ) nreset = 4;
    
    /* standard initialization */
    
    lastch = 0;
    lastdw = 0;
    lastup = 0;
    level  = one;
    tau    = tm1;
    iterma = MAXIT;
    epsx   = tm5;
    sigsm  = sqrt(epsmac);
    smalld = tm1;
    
    /* formerly tm2. smalld has much influence on the Maratos-effect */
    
    smallw = exp(two*log(epsmac)/three);
    rho    = tm6;
    rho1   = tm10;
    del01  = del0/tp1;
    delmin = patMin(del01,patMax(tm6*del0,smallw));
    if ( ! analyt ) delmin = patMin(del01,patMax(epsdif,delmin));
    c1d    = tm2;
    scfmax = tp4;
    taufac = tp1;
    taumax = pow(scfmax,2);
    updmy0 = tm1;
    
    /* take del0 and tau0 from block data or setup0 in function definition */
    /* may be modified by subsequent call of setup                         */
    
    j = 0;
    while ( name[j] == ' ' ) {
        j = j+1;
    }
    if ( name[j] == '\0' ) {
        strcpy(name1,xxx);
    } else {
        k = j+1;
        while ( name[k] != ' ' && name[k] != '\0' && k-j < 8 ) {
            k = k+1;
        }
        strncpy(name1,&name[j],k-j);
        name1[k-j] = '\0';
        
        for (i = 0 ; i <= k-j-1 ; i++) {
            iz = name1[i];
            if ( iz < 48 || ( iz > 57 && iz < 65 )
            || ( iz > 90 && iz < 97 ) || iz > 122 ) name1[i] = 'x';
        }
        if ( k - j < 8 ) strncat(name1,xxx,8-k+j);
    }
    if ( ! silent ) {
        strcpy(fil,name1);
        meu  = fopen(strcat(fil,".mes"),"w");
        strcpy(fil,name1);
        prou = fopen(strcat(fil,".pro"),"w");
    }
    
    infiny = epsmac/tolmac;
    fx     = zero;
    b2n    = zero;
    b2n0   = zero;
    nres   = ng+nh;
    if ( cold ) {
        for (i = 1 ; i <= NX ; i++) {
            for (j = 1 ; j <= NX ; j++) {
                a[j][i] = zero;
            }
            a[i][i]  = one;
            diag0[i] = one;
        }
    }
    for (i = 1 ; i <= NRESM ; i++) {
        diag[i] = zero;
        for (j = 1 ; j <= NX ; j++) {
            qr[j][i]   = zero;
            gres[j][i] = zero;
        }
    }
    for (i = 1 ; i <= NX ; i++) {
        xnull[i] = zero;
        ug[i]    = -infiny;
        og[i]    = infiny;
        llow[i]  = FALSE;
        lup[i]   = FALSE;
    }
    for (i = 1 ; i <= nh ; i++) {
        delfac[i] = one;
    }
    
    if ( bloc ) user_eval(xnull,0);
    
    for (i = nh+1 ; i <= nres ; i++) {
        delfac[i] = one;
        
        /* scan for real lower or upper bounds */
        
        if ( gunit[1][i] == 1 ) {
        
            esg(i-nh,xnull,&gxi);
            
            if ( gunit[3][i] > 0 ) {
                llow[gunit[2][i]] = TRUE;
                ug[gunit[2][i]]   = -gxi/gunit[3][i];
            } else {
                og[gunit[2][i]]   = -gxi/gunit[3][i];
                lup[gunit[2][i]]  = TRUE;
            }
        }
    }
    for (i = nh+1 ; i <= nres ; i++) {
    
        /* modify del0, such that lower and upper bound never become binding */
        /* simultaneously                                                    */
        
        if ( gunit[1][i] == 1 ) {
            j = gunit[2][i];
            if ( og[j] < infiny && ug[j] > -infiny ) {
                del0 = patMin(del0,(og[j]-ug[j])*tm1*myfabs(gunit[3][i]));
            }
        }
    }
    for (i = nh+1 ; i <= nres ; i++) {
    
        /* delfac corresponds to an indirect primal scaling */
        
        if ( gunit[1][i] == 1 ) {
            j = gunit[2][i];
            if ( gunit[3][i] > 0 ) {
                delfac[i] = patMax(delfac[i],myfabs(ug[j])*tm1);
                if ( og[j] < infiny )
                delfac[i] = patMin(delfac[i],(og[j]-ug[j])/(tp1*del0));
            } else {
                delfac[i] = patMax(delfac[i],myfabs(og[j])*tm1);
                if ( ug[j] > -infiny )
                delfac[i] = patMin(delfac[i],(og[j]-ug[j])/(tp1*del0));
            }
        }
    }
    bd0 = infiny;
    for (i = 1 ; i <= n ; i++) {
        if ( ug[i] > zero ) bd0 = patMin(bd0,og[i]);
        if ( og[i] < zero ) bd0 = patMin(bd0,-ug[i]);
    }
    /* change x if necessary , such that bounds not violated */
    
    corr = FALSE;
    
    /* evaluate gradients of supersimple functions only once              */
    /* a function is said to be supersimple if it is of the form a*x[j]+b */
    
    if ( gunit[1][0] == 1 ) {
        gconst[0] = TRUE;
        val[0]    = TRUE;
        for (i = 1 ; i <= n ; i++) {
            gradf[i] = zero;
        }
        gradf[gunit[2][0]] = gunit[3][0]*xsc[gunit[2][0]];
        gfn = myfabs(gunit[3][0]);
    } else {
        val[0] = FALSE;
        for (i = 1 ; i <= n ; i++) {
            gradf[i] = zero;
        }
    }
    for (i = 1 ; i <= nh ; i++) {
        if ( gunit[1][i] == 1 ) {
        
            /* a fixed variable. corrected if necessary */
            
            val[i]    = TRUE;
            gconst[i] = TRUE;
            gres[gunit[2][i]][i] = gunit[3][i]*xsc[gunit[2][i]];
            gresn[i] = myfabs(gunit[3][i])*xsc[gunit[2][i]];
            if ( gresn[i] == zero ) {
                if ( ! silent ) 
                fprintf(meu,"%i fixed variable/zero gradient\n", gunit[2][i]);
                fclose(meu);
                fclose(prou);
                exit(1);
            }
            esh(i,xnull,&hxi);
            
            term = -hxi/gunit[3][i];
            if ( term != x[gunit[2][i]] ) corr = TRUE;
            x[gunit[2][i]] = term;
        }
    }
    if ( bloc ) user_eval(x,0);
    
    for (i = nh+1 ; i <= nres ; i++) {
        if ( gunit[1][i] == 1 ) {
            if ( gunit[3][i] == 0 ) {
                if ( ! silent ) 
                fprintf(meu,"%i bounded variable, zero gradient\n",gunit[2][i]);
                fclose(meu);
                fclose(prou);
                exit(1);
            }
            esg(i-nh,x,&gxi);
            
            gxi = two*delmin-gxi;
            if ( gxi > zero ) {
                corr = TRUE;
                x[gunit[2][i]] = x[gunit[2][i]]+gxi/gunit[3][i];
            }
            gres[gunit[2][i]][i] = gunit[3][i]*xsc[gunit[2][i]];
            gresn[i]  = myfabs(gunit[3][i])*xsc[gunit[2][i]];
            val[i]    = TRUE;
            gconst[i] = TRUE;
        }
    }
    if ( corr && ! silent ) o8msg(13);

    /* remember to reevaluate the functions if corr = TRUE and bloc = TRUE */

    if ( bloc ) user_eval(x,1);

    for (i = 1 ; i <= nres ; i ++) {
        bind[i]  = 0;
        bind0[i] = 0;
        u[i]     = zero;
        u0[i]    = zero;
        cres[i]  = 0;
        cgres[i] = 0;
        
        /* initial weights of the penalty-term */
        
        if ( cold ) o8w[i] = one;
        DONLP2sort[i] = i;
    }   
    clow = one;
    ny   = two;

    /* scf  = weight factor for objective function    */
    /* scf0 = damping factor for tangential direction */
    
    scf    = one;
    scf0   = one;
    sigla  = twop11;
    beta   = four;  /* formerly two */
    alpha  = tm1;
    delta1 = p9;
    delta  = tm3;   /* delta = tm2  formerly */
    theta  = p9;    /* theta = 0.99 formerly */
    icf    = 0;
    icgf   = 0;
    if ( ! silent ) {
        fprintf(prou,"donlp2, v3, 05/29/98, copyright P. Spellucci\n");
        
        time(&tim);
        
        fprintf(prou,"%s",ctime(&tim));
        fprintf(prou,"%s\n", name);
        
        fprintf(meu, "donlp2, v3, 05/29/98, copyright P. Spellucci\n");
        fprintf(meu, "%s",ctime(&tim));
        fprintf(meu, "%s\n", name);
    }
    return;
}

/* **************************************************************************** */
/*                        do final solution check and output                    */
/* **************************************************************************** */
void o8fin(void) {

    void    solchk(void);

    static INTEGER  i,j,k,ih1,ih2,ih3,ih5,ih6,ih7,ih8,ih9;
    static DOUBLE   umin,term;
    static INTEGER  nsing,crtot,cgrtot,nupreg,nbfgs1,nresta;
    patString line ;

    /* termination reason + 11 = message number */
    
    vector<patString> messag(19) ;
    messag[0] = "" ;
    messag[1] = "constraint evaluation returns error with current point" ;
    messag[2] = "objective evaluation returns error with current point" ;
    messag[3] = "QPsolver: working set singular in dual extended QP " ;
    messag[4] = "QPsolver: extended QP-problem seemingly infeasible " ;
    messag[5] = "QPsolver: no descent direction from QP for tau=tau_max" ;
    messag[6] = "QPsolver: on exit correction small, infeasible point";
    messag[7] = "stepsizeselection: computed d not a direction of descent" ;
    messag[8] = "more than MAXIT iteration steps" ;
    messag[9] = "stepsizeselection: no acceptable stepsize in [sigsm,sigla]" ;
    messag[10] = "stepsizeselection: directional deriv. very small, infeasible" ;
    messag[11] = "KT-conditions satisfied, no further correction computed" ;
    messag[12] = "KT-conditions satisfied, computed correction small" ;
    messag[13] = "stepsizeselection: x (almost) feasible, dir. deriv. very small" ;
    messag[14] = "KT-conditions (relaxed) satisfied, singular point" ;
    messag[15] = "very slow primal progress, singular or illconditoned problem" ;
    messag[16] = "very slow progress in x, singular problem" ;
    messag[17] = "correction very small, almost feasible but singular point" ;
    messag[18] = "patMax(n,10) small differences in penalty function ;terminate" ;

    if ( scf != zero ) {
        for (i = 1 ; i <= nres ; i++) {
            u[i] = u[i]/scf;
        }
    }
    /* in solchk the user might add some additional checks and/or output */
    
    solchk();

    if ( silent && ! intakt ) return;

    if ( intakt && ! silent ) printf("%s\n", name);
    
    if ( ! silent ) {
        if ( intakt ) {
            printf(  "\n     n= %9i    nh= %9i    ng= %9i\n", n,nh,ng);
            printf(  "\n  epsx= %9.3e sigsm= %9.3e\n"       , epsx,sigsm);
            printf(  "\nstartvalue\n");
            for (i = 1 ; i <= n ; i++) {
                printf(  " %14.7e ", xst[i]);
                if ( i % 5 == 0 || i == n ) printf(  "\n");
            }
        }
        fprintf(prou,"\n     n= %9i    nh= %9i    ng= %9i\n", n,nh,ng);
        fprintf(prou,"\n  epsx= %9.3e sigsm= %9.3e\n"       , epsx,sigsm);
        fprintf(prou,"\nstartvalue\n");
        for (i = 1 ; i <= n ; i++) {
            fprintf(prou," %14.7e ", xst[i]);
            if ( i % 5 == 0 || i == n ) fprintf(prou,"\n");
        }
    }
    if ( intakt && ! silent ) {
        printf("\n  eps= %9.2e  tol= %9.2e del0= %9.2e delm= %9.2e tau0= %9.2e\n",
        epsmac,tolmac,del0,delmin,tau0);
        printf(  "  tau= %9.2e   sd= %9.2e   sw= %9.2e  rho= %9.2e rho1= %9.2e\n",
        tau,smalld,smallw,rho,rho1);
    }
    if ( ! silent ) { 
        fprintf(prou,
               "\n  eps= %9.2e  tol= %9.2e del0= %9.2e delm= %9.2e tau0= %9.2e\n",
        epsmac,tolmac,del0,delmin,tau0);
        fprintf(prou,
                 "  tau= %9.2e   sd= %9.2e   sw= %9.2e  rho= %9.2e rho1= %9.2e\n",
        tau,smalld,smallw,rho,rho1);
    }
    if ( ! silent ) {
        fprintf(prou,
        " scfm= %9.2e  c1d= %9.2e epdi= %9.2e\n  nre= %9i anal= %9i\n",
        scfmax,c1d,epsdif,nreset,analyt);
        if ( intakt ) printf(
        " scfm= %9.2e  c1d= %9.2e epdi= %9.2e\n  nre= %9i anal= %9i\n",
        scfmax,c1d,epsdif,nreset,analyt);
    }
    if ( ! silent && ! analyt ) {
        fprintf(prou," vbnd= %9.2e efcn= %9.2e diff=%1i\n"
        , taubnd,epsfcn,difftype);
        if ( intakt ) printf("taubnd= %9.2e epsfcn= %9.2e difftype=%1i\n"
        , taubnd,epsfcn,difftype);
    }
    i    = 0;
    j    = 0;
    umin = zero;
    for (k = 1 ; k <= nres ; k++) {
        i = i+ cres[k];
        j = j+cgres[k];
        if ( k > nh ) umin = patMin(umin,u[k]);
    }
    crtot  = i;
    cgrtot = j;
    nsing  = 0;
    nresta = 0;
    nupreg = 0;
    nbfgs1 = 0;
    for (k = 1 ; k <= itstep ; k++) {
        if ( accinf[k][10] == one ) nsing  = nsing+1;
        if ( accinf[k][27] == one ) nbfgs1 = nbfgs1+1;
        
        /* for the Pantoja Mayne update */
        
        if ( accinf[k][29] == zero && accinf[k][27] == one ) nupreg = nupreg+1;
        if ( accinf[k][27] == -one ) nresta = nresta+1;
    }
    k = (int)optite+11;
    if ( k >= 1 && k <= 18 ) {
        line = messag[k];
    } else {
        line = "variable optite undefined on exit";
    }
    if ( intakt && ! silent ) {
      printf(      "\n termination reason:\n %s\n",line.c_str());
        printf(        " evaluations of f                    %9i\n",    icf);
        printf(        " evaluations of grad f               %9i\n",    icgf);
        printf(        " evaluations of constraints          %9i\n",    crtot);
        printf(        " evaluations of grads of constraints %9i\n",    cgrtot);
        printf(        " final scaling of objective          %13.6e\n", scf);
        printf(        " norm of grad(f)                     %13.6e\n", gfn);
        printf(        " lagrangian violation                %13.6e\n", b2n);
        printf(        " feasibility violation               %13.6e\n", upsi);
        printf(        " dual feasibility violation          %13.6e\n", umin);
        printf(        " optimizer runtime sec's             %13.6e\n", runtim);
    }
    if ( ! silent ) {
      fprintf(prou,"\n termination reason:\n %s\n",line.c_str());
        fprintf(prou,  " evaluations of f                    %9i\n",    icf);
        fprintf(prou,  " evaluations of grad f               %9i\n",    icgf);
        fprintf(prou,  " evaluations of constraints          %9i\n",    crtot);
        fprintf(prou,  " evaluations of grads of constraints %9i\n",    cgrtot);
        fprintf(prou,  " final scaling of objective          %13.6e\n", scf);
        fprintf(prou,  " norm of grad(f)                     %13.6e\n", gfn);
        fprintf(prou,  " lagrangian violation                %13.6e\n", b2n);
        fprintf(prou,  " feasibility violation               %13.6e\n", upsi);
        fprintf(prou,  " dual feasibility violation          %13.6e\n", umin);
        fprintf(prou,  " optimizer runtime sec's             %13.6e\n", runtim);
    }                                         
    if ( intakt && ! silent ) printf("\n\n optimal value of f =  %21.14e\n", fx);
    if ( ! silent ) fprintf(prou,    "\n\n optimal value of f =  %21.14e\n", fx);
    
    if ( intakt && ! silent ) {
        printf(      "\n optimal solution  x =\n");
        for (i = 1 ; i <= n ; i++) {
            printf(      " %21.14e", x[i]);
            if ( i % 3 == 0 || i == n ) printf(     "\n");
        }
    }
    if ( ! silent ) {
        fprintf(prou,"\n optimal solution  x =\n");
        for (i = 1 ; i <= n ; i++) {
            fprintf(prou," %21.14e", x[i]);
            if ( i % 3 == 0 || i == n ) fprintf(prou,"\n");
        }
    }
    if ( nres != 0 && ! silent ) {
        fprintf(prou,"\n  multipliers are relativ to scf=1\n");
        fprintf(prou,  "  nr.    constraint      normgrad (or 1)   multiplier\n");
        for (i = 1 ; i <= nres ; i++) {
            fprintf(prou," %4i  %14.7e   %14.7e   %14.7e\n"
            , i,res[i],gresn[i],u[i]);
        }
        if ( intakt ) {
            printf(  "\n  multipliers are relativ to scf=1\n");
            printf(    "  nr.    constraint      normgrad (or 1)   multiplier\n");
            for (i = 1 ; i <= nres ; i++) {
                printf(  " %4i  %14.7e   %14.7e   %14.7e\n"
                , i,res[i],gresn[i],u[i]);
            }
        }
    }
    if ( nres != 0 && ! silent ) {
        fprintf(prou,"\n evaluations of restrictions and their gradients\n");
        for (i = 1 ; i <= nres ; i++) {
            fprintf(prou," (%6i,%6i)", cres[i],cgres[i]);
            if ( i % 5 == 0 || i == nres ) fprintf(prou,"\n");
        }
        if ( intakt ) {
            printf(  "\n evaluations of restrictions and their gradients\n");
            for (i = 1 ; i <= nres ; i++) {
                printf(  " (%6i,%6i)", cres[i],cgres[i]);
                if ( i % 5 == 0 || i == nres ) printf(  "\n");
            }
        }
    }
    if ( itstep > 1 && optite == 0 ) itstep = itstep-1;
    if ( ! silent )
    fprintf(prou,    "last estimate of condition of active gradients %10.3e\n",
    accinf[itstep][13]);
                    
    term = accinf[itstep][14];
    i    = itstep;
    while ( i > 1 && term == zero ) {
        i    = i-1;
        term = accinf[i][14];
    }
    if ( ! silent ) {
        fprintf(prou,"last estimate of condition of approx. hessian  %10.3e\n",
        term);
        fprintf(prou,"iterative steps total           %5i\n", itstep);
        fprintf(prou,"# of restarts                   %5i\n", nresta);
        fprintf(prou,"# of full regular updates       %5i\n", nupreg);
        fprintf(prou,"# of updates                    %5i\n", nbfgs1);
        fprintf(prou,"# of full regularized SQP-steps %5i\n", nsing);

        if ( intakt ) {
            printf(  "last estimate of cond.nr. of active gradients  %10.3e\n",
            accinf[itstep][13]);
            printf(  "last estimate of cond.nr. of approx.  hessian  %10.3e\n",
            accinf[itstep][14]);
            printf(  "iterative steps total           %5i\n", itstep);
            printf(  "# of restarts                   %5i\n", nresta);
            printf(  "# of full regular updates       %5i\n", nupreg);
            printf(  "# of updates                    %5i\n", nbfgs1);
            printf(  "# of regularized full SQP-steps %5i\n", nsing);
        }
    }
    if ( optite < zero ) te1 = TRUE;
    if ( silent )        te1 = FALSE;
    if ( te1 ) {
        for (i = 1 ; i <= itstep ; i++) {
            ih1 = accinf[i][1];
            ih2 = accinf[i][9];
            ih3 = accinf[i][10];
            ih5 = accinf[i][18];
            ih6 = accinf[i][19];
            ih7 = accinf[i][22];
            ih8 = accinf[i][26];
            ih9 = accinf[i][27];
            
            fprintf(prou, 
             "%4i  fx= %13.6e scf= %13.6e psi= %13.6e ups= %13.6e\n",
            ih1,accinf[i][2],accinf[i][3],accinf[i][4],accinf[i][5]);
            fprintf(prou,
            "     del= %13.6e b20= %13.6e b2n= %13.6e  nr=%5i\n",
            accinf[i][6],accinf[i][7],accinf[i][8],ih2);
            fprintf(prou,
            "      si=%4i            u-= %13.6e c-r= %13.6e c-d= %13.6e\n",
            ih3,accinf[i][11],accinf[i][13],accinf[i][14]);
            fprintf(prou,
            "      xn= %13.6e  dn= %13.6e pha=%4i            cl=%14i\n",
            accinf[i][16],accinf[i][17],ih5,ih6);
            fprintf(prou,
            "     skm= %13.6e sig= %13.6e cf+=%5i          dir= %13.6e\n",
            accinf[i][20],accinf[i][21],ih7,accinf[i][23]);
            fprintf(prou,
            "     dsc= %13.6e cos= %13.6e vio=%5i\n",
            accinf[i][24],accinf[i][25],ih8);
            fprintf(prou,
            "     upd=%5i           tk= %13.6e xsi= %13.6e\n",
            ih9,accinf[i][28],accinf[i][29]);
            
            if ( accinf[i][10] == 1. ) {
                fprintf(prou,"     qpt= %13.0e tqp= %13.6e sla= %13.6e\n",
                accinf[i][30],accinf[i][31],accinf[i][32]);
            }
        }
    }
    /*  accinf a c c u m u l a t e d   i n f o r m a t i o n                    */
    /*  on iteration sequence                                                   */
    /*  1: step-nr                                                              */
    /*  2: f(x-k) current value of objective (zero in feasibility improvement   */
    /*            phase (-1) )                                                  */
    /*  3: scf    internal scaling of objective (zero in phase -1)              */
    /*  4: psi    the weighted penalty-term                                     */
    /*  5: upsi   the unweighted penalty-term (L1-norm of constraint vector)    */
    /*  6: del_k_1 bound for currently active constraints                       */
    /*  7: b2n0   L2-norm of projected gradient, based on constraints in level  */
    /*            delmin and below, measured in the norm induced by the         */
    /*            inverse hessian estimate                                      */
    /*  8: b2n    L2-norm of projected gradient based on del_k_1                */
    /*  9: nr     number of binding constraints                                 */
    /* 10: sing   if 1, the binding constraints don't satisfy the regularity    */
    /*            condition                                                     */
    /* 11: umin   infinity norm of negative part of multiplier                  */
    /* 12: -------------                                                        */
    /* 13: cond_r condition number of diagonal part of QR-decomp. of normalized */
    /*            gradients of binding constraints                              */
    /* 14: cond_h condition number of diagonal of Cholesky-factor of updated    */
    /*            full hessian                                                  */
    /* 15: scf0   the relative damping of tangential component if upsi>tau0/2   */
    /* 16: xnorm  L2-norm of x                                                  */
    /* 17: dnorm  L2-norm of d (correction from QP -subproblem, unscaled)       */
    /* 18: phase  -1 : infeasibility improvement phase, 0: initial optimization */
    /*            1  : binding constraints unchanged , 2: d small, Maratos      */
    /*                 correction in use                                        */
    /* 19: c_k    number of decreases of penalty weights                        */
    /* 20: wmax   infinity norm of weights                                      */
    /* 21: sig_k  stepsize from unidimensional minimization (backtracking)      */
    /* 22: cfincr number of objective evaluations for stepsize-algorithm        */
    /* 23: dirder directional derivative of penalty-function along d (scaled)   */
    /* 24: dscal  scaling factor for d                                          */
    /* 25: cosphi cos of arc between d and d_previous. if larger theta , sig    */
    /*            larger than one (up to sigla) is tried                        */
    /* 26: violis[0] number of constraints not binding at x but hit during      */
    /*            line search                                                   */
    /* 27:        type of update for h: 1 normal P&M-BFGS-update,               */
    /*            0 update suppressed, -1 restart with scaled unit matrix ,     */
    /*            2 standard BFGS, 3 BFGS modified by Powells device            */
    /* 28: ny_k/tk modification factor for damping the projector in BFGS/       */
    /*            Pantoja-Mayne-term respectively                               */
    /* 29: 1-my_k/xsik modification factor for damping the quasi-Newton-        */
    /*            relation in BFGS for unmodified BFGS ny_k should be larger    */
    /*            than updmy0 (near one) and 1-my_k equal one./Pantoja-Mayne    */
    /*            term respectively                                             */
    /* 30: qpterm 0, if sing = -1, termination indicator of QP-solver otherwise */
    /*            1: successful, -1: tau becomes larger than tauqp without      */
    /*            slack-variables becoming sufficiently small                   */
    /*            -3: working set of QP-solver becomes linearly dependent       */
    /*            -2: infeasible QP-problem (theoretically impossible)          */
    /* 31: tauqp  weight  of slack-variables in QP-solver                       */
    /* 32: infeas L1-norm of slack-variables in QP-solver                       */

    if ( ! silent ) fclose(prou);
    if ( ! silent ) fclose(meu);
    
    return;
}
/* **************************************************************************** */
/*                      prints informations if te2 = TRUE                       */
/* **************************************************************************** */
void o8info(INTEGER icase) {

    void o8mdru (DOUBLE a[][NRESM+1],INTEGER ma,INTEGER na,char head[],
                 FILE *lognum,LOGICAL fix);
    void o8mdru_(DOUBLE a[][NX+1],   INTEGER ma,INTEGER na,char head[],
                 FILE *lognum,LOGICAL fix);
            
    static INTEGER  i,j,l,k;
    static DOUBLE   y,phih;
    static char     head[41];

    if(!te2) return;
    
    switch (icase) {
        
    case 1:
        fprintf(prou,"\n\n\n");
        for (i = 1 ; i <= 80 ; i++) fprintf(prou,"=");
        fprintf(prou,"\n          %4i-th iteration step\n", itstep);
        fprintf(prou,"   scf= %11.4e psist= %11.4e   psi= %11.4e  upsi= %11.4e\n", 
        scf,psist,psi,upsi);
        fprintf(prou,"  fxst= %11.4e    fx= %11.4e\n", fxst,fx);
        fprintf(prou,"  x=\n");
        for (i = 1 ; i <= n ; i++) {
            fprintf(prou,"  %11.4e", x[i]);
            if ( i % 6 == 0 || i == n ) fprintf(prou,"\n");
        }
        fprintf(prou," valid permutation of x\n\n");
        
        for (i = 1 ; i <= n ; i++) {
            fprintf(prou,"%3i ", perm[i]);
            if ( i % 20 == 0 || i == n ) fprintf(prou,"\n");
        }
        if ( phase >= 0 && te3 ) {
            strcpy(head,"quasi-Newton-matrix full update");
            
            o8mdru_(a,n,n,head,prou,FALSE);
        }
        if ( intakt ) {
            printf(  "\n\n\n");
            for (i = 1 ; i <= 80 ; i++) printf(  "=");
            printf(  "\n          %4i-th iteration step\n", itstep);
            printf(  "   scf= %11.4e psist= %11.4e   psi= %11.4e  upsi= %11.4e\n", 
            scf,psist,psi,upsi);
            printf(  "  fxst= %11.4e    fx= %11.4e\n", fxst,fx);
            printf(  "  x=\n");
            for (i = 1 ; i <= n ; i++) {
                printf("  %11.4e", x[i]);
                if ( i % 6 == 0 || i == n ) printf(  "\n");
            }
            printf(  " valid permutation of x\n\n");
        
            for (i = 1 ; i <= n ; i++) {
                printf(  "%3i ", perm[i]);
                if ( i % 20 == 0 || i == n ) printf(  "\n");
            }
        }
        return;
        
    case 2:
      fprintf(prou,"\n\n  del= %12.5e  b2n0= %12.5e   b2n= %12.5e   gfn= %12.5e\n", 
      del,b2n0,b2n,gfn);
        
        if ( DONLP2alist[0] != 0 ) {
            fprintf(prou,"\n\n values of restrictions\n ");
            for (i = 1 ; i <= DONLP2alist[0] ; i++) {
                fprintf(prou,"(%4i   %11.4e   %11.4e)  ", 
                DONLP2alist[i],res[DONLP2alist[i]],gresn[DONLP2alist[i]]);
                if ( i % 2 == 0 || i == DONLP2alist[0] ) fprintf(prou,"\n ");
            }
        }
        if ( DONLP2alist[0] != 0 && ! singul ) {
            fprintf(prou,"\n\n   diag[r]=\n");
            for (i = 1 ; i <= DONLP2alist[0] ; i++) {
                fprintf(prou,"  %11.4e", diag[i]);
                if ( i % 6 == 0 || i == DONLP2alist[0] ) fprintf(prou,"\n");
            }
        }
        if ( DONLP2alist[0] != 0 && te3 ) {
            for (i = 1 ; i <= DONLP2alist[0] ; i++) {
                l = DONLP2alist[i];
                fprintf(prou,"\n\n gradient of restriction nr.%4i\n ", l);
                for (j = 1 ; j <= n ; j++) {
                    fprintf(prou," %11.4e  ", gres[j][l]);
                    if ( j % 5 == 0 || j == n ) fprintf(prou,"\n ");
                }
            }
        }
        if ( intakt ) {
            printf("\n\n  del= %12.5e  b2n0= %12.5e   b2n= %12.5e   gfn= %12.5e\n", 
            del,b2n0,b2n,gfn);
        
            if ( DONLP2alist[0] != 0 ) {
                printf(  "\n\n values of restrictions\n ");
                for (i = 1 ; i <= DONLP2alist[0] ; i++) {
                    printf(  "(%4i   %11.4e   %11.4e)  ", 
                    DONLP2alist[i],res[DONLP2alist[i]],gresn[DONLP2alist[i]]);
                    if ( i % 2 == 0 || i == DONLP2alist[0] ) printf(  "\n ");
                }
            }
            if ( DONLP2alist[0] != 0 && ! singul ) {
                printf(  "\n\n   diag[r]=\n");
                for (i = 1 ; i <= DONLP2alist[0] ; i++) {
                    printf(  "  %11.4e", diag[i]);
                    if ( i % 6 == 0 || i == DONLP2alist[0] ) printf(  "\n");
                }
            }
            if ( DONLP2alist[0] != 0 && te3 ) {
                for (i = 1 ; i <= DONLP2alist[0] ; i++) {
                    l = DONLP2alist[i];
                    printf(  "\n\n gradient of restriction nr.%4i\n ", l);
                    for (j = 1 ; j <= n ; j++) {
                        printf(  " %11.4e  ", gres[j][l]);
                        if ( j % 5 == 0 || j == n ) printf(  "\n ");
                    }
                }
            }
        }
        return;
        
    case 3:
        if( ! (nr == 0 || phase == -1) ) {
            fprintf(prou,"\n  multipliers: first estimate\n  u =\n");
            for (k = 1 ; k <= nr ; k++) {
                fprintf(prou," %4i  %11.4e", DONLP2alist[k],u[DONLP2alist[k]]);
                if ( k % 4 == 0 || k == nr ) fprintf(prou,"\n");
            }
        }
        if( ! (nr == 0 || phase == -1) && intakt ) {
            printf(      "\n  multipliers: first estimate\n  u =\n");
            for (k = 1 ; k <= nr ; k++) {
                printf(      " %4i  %11.4e", DONLP2alist[k],u[DONLP2alist[k]]);
                if ( k % 4 == 0 || k == nr ) printf(      "\n");
            }
        }
        return;
        
    case 4:
        if( ! (nr == 0 || phase == -1) ) {
            fprintf(prou,"\n  multipliers: second estimate\n  u =\n");
            for (k = 1 ; k <= nr ; k++) {
                fprintf(prou," %4i  %11.4e", DONLP2alist[k],u[DONLP2alist[k]]);
                if ( k % 4 == 0 || k == nr ) fprintf(prou,"\n");
            }
        }
        if( ! (nr == 0 || phase == -1) && intakt ) {
            printf(      "\n  multipliers: second estimate\n  u =\n");
            for (k = 1 ; k <= nr ; k++) {
                printf(      " %4i  %11.4e", DONLP2alist[k],u[DONLP2alist[k]]);
                if ( k % 4 == 0 || k == nr ) printf(      "\n");
            }
        }
        return;
        
    case 5:
        if( intakt )
        printf(          "  condition number of r     %.15e\n",accinf[itstep][13]);
        fprintf(prou,    "  condition number of r     %.15e\n",accinf[itstep][13]);
        if ( phase == -1 ) {
        
            return;
            
        } else {
            fprintf(prou,"  condition number of a     %.15e\n",accinf[itstep][14]);
            if ( intakt ) {
                printf(  "  condition number of a     %.15e\n",accinf[itstep][14]);
            }
            return;
        }
    case 6:
    
        return;
        
    case 7:
        fprintf(prou,"\n\n  phase=%3i  scf0= %11.4e\n", phase,scf0);
        fprintf(prou,    "  d =\n");
        for (i = 1 ; i <= n ; i++) {
            fprintf(prou,"  %11.4e", d[i]);
            if ( i % 6 == 0 || i == n ) fprintf(prou,"\n");
        }
        if ( phase == 2 ) {
            fprintf(prou,"\n\n  dd=\n");
            for (i = 1 ; i <= n ; i++) {
                fprintf(prou,"  %11.4e", dd[i]);
                if ( i % 6 == 0 || i == n ) fprintf(prou,"\n");
            }
        }
        if ( intakt ) {
            printf(  "\n\n  phase=%3i  scf0= %11.4e\n", phase,scf0);
            printf(      "  d =\n");
            for (i = 1 ; i <= n ; i++) {
                printf(  "  %11.4e", d[i]);
                if ( i % 6 == 0 || i == n ) printf(  "\n");
            }
            if ( phase == 2 ) {
                printf(  "\n\n  dd=\n");
                for (i = 1 ; i <= n ; i++) {
                    printf(  "  %11.4e", dd[i]);
                    if ( i % 6 == 0 || i == n ) printf(  "\n");
                }
            }
        }
        return;
        
    case 8:
        y    = tau0*p5;
        phih = fx*scf+psi;
    
        if ( intakt ) {
            printf(  "\n\n start unimin\n\n");
            printf(  "    phi= %11.4e   dphi= %11.4e    psi= %11.4e tau0/2= %11.4e\n",
            phih,dirder,psi,y);
            printf(  "     fx= %11.4e  dscal= %11.4e    scf= %11.4e   upsi= %11.4e\n", 
            fx,dscal,scf,upsi);
        }
        fprintf(prou,"\n\n start unimin\n\n");
        fprintf(prou,"    phi= %11.4e   dphi= %11.4e    psi= %11.4e tau0/2= %11.4e\n",
        phih,dirder,psi,y);
        fprintf(prou,"     fx= %11.4e  dscal= %11.4e    scf= %11.4e   upsi= %11.4e\n",
        fx,dscal,scf,upsi);
        
        return;
        
    case 9:
        fprintf(prou,"    sig= %11.4e     fx= %11.4e    psi= %11.4e   upsi= %11.4e\n",
        sig,fx1,psi1,upsi1);
        if ( intakt ) 
        printf(      "    sig= %11.4e     fx= %11.4e    psi= %11.4e   upsi= %11.4e\n", 
        sig,fx1,psi1,upsi1);
        
        return;
        
    case 10:
        fprintf(prou,"\n\n end unimin\n");
        fprintf(prou,  "\n sig= %11.4e  num. f-evaluations%2i\n", sig,cfincr);
        fprintf(prou,    " list of inactive hit constraints\n");
        for (i = 1 ; i <= violis[0] ; i++) {
            fprintf(prou,"%4i  ", violis[i]);
            if ( i % 13 == 0 || i == violis[0] ) fprintf(prou,"\n");
        }
        if ( violis[0] == 0 ) fprintf(prou,"none\n");
        if ( intakt ) {
            printf(  "\n\n end unimin\n");
            printf(    "\n sig= %11.4e  num. f-evaluations%2i\n", sig,cfincr);
            printf(      " list of inactive hit constraints\n");
            for (i = 1 ; i <= violis[0] ; i++) {
                printf(  "%4i  ", violis[i]);
                if ( i % 13 == 0 || i == violis[0] ) printf(  "\n");
            }
            if ( violis[0] == 0 ) printf("none\n");
        }
        return;
        
    case 11:
        if ( intakt ) printf("additional increase of eta due to large clow\n");
        fprintf(prou,        "additional increase of eta due to large clow\n");
        
        return;
        
    case 12:
        fprintf(prou,
        "\n\n  current scaling,  scf =  %11.4e clow = %12i eta =  %11.4e\n", 
        scf,clow,eta);
        
        if(nres != 0) {
            fprintf(prou,"\n\n  scalres=\n");
            for (i = 1 ; i <= nres ; i++) {
                fprintf(prou,"  %11.4e", o8w[i]);     
                if ( i % 6 == 0 || i == nres ) fprintf(prou,"\n");
            }
        }
        if ( intakt ) {
            printf(
            "\n\n  current scaling,  scf =  %11.4e clow = %12i eta =  %11.4e\n", 
            scf,clow,eta);
            
            if ( nres != 0 ) {
                printf(  "\n\n  scalres=\n");
                for (i = 1 ; i <= nres ; i++) {
                    printf(  "  %11.4e", o8w[i]);
                    if ( i % 6 == 0 || i == nres ) printf(  "\n");
                }
            }
        }
        return;
        
    case 13:
        if ( accinf[itstep][27] == zero ) {
            if ( intakt ) printf("update suppressed\n");
            fprintf(prou,        "update suppressed\n");
        } else if ( accinf[itstep][27] == -one ) {
            fprintf(prou,        "restart with scaled unit matrix\n");
            if ( intakt ) printf("restart with scaled unit matrix\n");
        } else {
            fprintf(prou,"BFGS-update\n");
            fprintf(prou," type = %14.6e\n", accinf[itstep][27]);
            fprintf(prou,"  ny  = %14.6e\n", accinf[itstep][28]);
            fprintf(prou," thet = %14.6e\n", accinf[itstep][29]);
            if ( intakt ) {
                printf(  "BFGS-update\n");
                printf(  " type = %14.6e\n", accinf[itstep][27]);
                printf(  "  ny  = %14.6e\n", accinf[itstep][28]);
                printf(  " thet = %14.6e\n", accinf[itstep][29]);
            }
        }
        return;
        
    case 14:
        if ( accinf[itstep][27] == zero ) {
            if ( intakt ) printf("update suppressed\n");
            fprintf(prou,        "update suppressed\n");
        } else if ( accinf[itstep][27] == one ) {
            fprintf(prou,"BFGS-update as in Pantoja and Mayne\n");
            fprintf(prou,"  tk  = %14.6e\n", accinf[itstep][28]);
            fprintf(prou," xsik = %14.6e\n", accinf[itstep][29]);
            if ( intakt ) {
                printf(  "BFGS-update\n");
                printf(  "  tk  = %14.6e\n", accinf[itstep][28]);
                printf(  " xsik = %14.6e\n", accinf[itstep][29]);
            }
        } else {
            if ( intakt ) printf("restart with scaled unit matrix\n");
            fprintf(prou,        "restart with scaled unit matrix\n");
        }
        return;
        
    case 15:
        if ( intakt ) printf("\n\n\n singular case : full regularized SQP\n");
        fprintf(prou,        "\n\n\n singular case : full regularized SQP\n");
        if ( intakt ) printf("  del = %.15e\n", del);
        fprintf(prou,        "  del = %.15e\n", del);
        
        if ( intakt ) { 
            printf(  "\n\n  scalres=\n");
            for (i = 1 ; i <= nres ; i++) {
                printf(  "  %11.4e", o8w[i]);
                if ( i % 6 == 0 || i == nres ) printf(  "\n");
            }
        }
        fprintf(prou,"\n\n  scalres=\n");
        for (i = 1 ; i <= nres ; i++) {
            fprintf(prou,"  %11.4e",o8w[i]);
            if ( i % 6 == 0 || i == nres ) fprintf(prou,"\n");
        }
        if ( te3 ) {
            strcpy(head,"gradients of constraints");
            
            o8mdru(gres,n,nres,head,prou,FALSE);
        }
        return;
        
    case 16:
        fprintf(prou,"exit from full SQP\n");
        fprintf(prou,"            termination reason  %7.0e\n", 
        accinf[itstep][30]);
        fprintf(prou,"          final value of tauqp  %10.3e\n",
        accinf[itstep][31]);
        fprintf(prou,"      sum norm of slack vector  %10.3e\n",
        accinf[itstep][32]);
                                
        fprintf(prou,"\n\n  phase=%3i  scf0= %11.4e\n",phase,scf0);
        for (i = 1 ; i <= n ; i++) {
            fprintf(prou,"  %11.4e", d[i]);
            if ( i % 6 == 0 || i == n ) fprintf(prou,"\n");
        }
        if ( nres != 0 ) { 
            fprintf(prou,"\n  multipliers: first estimate\n  u =\n");
            for (k = 1 ; k <= nres ; k++) {
                fprintf(prou," %4i  %11.4e", k,u[k]);       
                if ( k % 4 == 0 || k == nres ) fprintf(prou,"\n");
            }
        }
        if ( intakt ) {
            printf(  "exit from full SQP\n");
            printf(  "            termination reason  %7.0e\n", 
            accinf[itstep][30]);
            printf(  "          final value of tauqp  %10.3e\n",
            accinf[itstep][31]);
            printf(  "      sum norm of slack vector  %10.3e\n",
            accinf[itstep][32]);
                                
            printf("\n\n  phase=%3i  scf0= %11.4e\n",  phase,scf0);
            for (i = 1 ; i <= n ; i++) {
                printf(  "  %11.4e", d[i]);
                if ( i % 6 == 0 || i == n ) printf(  "\n");
            }
            if ( nres != 0 ) {
                printf(  "\n  multipliers: first estimate\n  u =\n");
                for (k = 1 ; k <= nres ; k++) {
                    printf(  " %4i  %11.4e", k,u[k]);
                    if ( k % 4 == 0 || k == nres ) printf(  "\n");
                }
            }
        }
        return;
        
    case 17:
        fprintf(prou,"small directional derivative %.15e: finish\n",dirder);
        if ( intakt )
        printf(      "small directional derivative %.15e: finish\n",dirder);
        
        return;
        
    case 18:
        if ( intakt )
        printf(      "small correction from full regularized SQP,finish\n");
        fprintf(prou,"small correction from full regularized SQP,finish\n");
        
        return;
        
    case 19:
        fprintf(prou,        "QP-solver terminated unsuccessfully\n");
        if ( intakt ) printf("QP-solver terminated unsuccessfully\n");
        
        return;
        
    case 20:
        if ( intakt ) printf("restart with scaled unit matrix\n");
        fprintf(prou,        "restart with scaled unit matrix\n");
        
        return;
        
    case 21:
    
        return;
        
    case 22:
    
        return;
    }
}

/* **************************************************************************** */
/*         computation of new scaling factors for L1-penalty-function           */
/* **************************************************************************** */
void o8sce(void) {

    void    o8info(INTEGER icase);
    
    static INTEGER  i;
    static DOUBLE   term,s1,s2,diff0;
    static LOGICAL  wlow;

    wlow = FALSE;
    for (i = 1 ; i <= nres ; i++) {
    
        /* w1 tentative new weights */
        
        term = ny*myfabs(u[i])+tau;
        if ( term > o8w[i] ) {
            w1[i] = term+tau;
        } else {
            w1[i] = o8w[i];
            if ( term < o8w[i]*p5 && bind[i] == 1 ) w1[i] = (term+o8w[i])*p5;
        }
        if ( w1[i] < o8w[i]) wlow = TRUE;
    }
    /* wlow equals TRUE if one tentative weight at least has been decreased */
    
    s1 = zero;
    s2 = zero;
    for (i = 1 ; i <= nres ; i++) {
        if ( i <= nh ) {
            s1 = s1+w1[i]*myfabs(resst[i]);
            s2 = s2+w1[i]*myfabs(res[i]);
        } else {
            s1 = s1-patMin(zero,resst[i])*w1[i];
            s2 = s2-patMin(zero,res[i]  )*w1[i];
        }
    }
    diff0 = (fxst-fx)*scf+(s1-s2);
    if ( wlow && diff0 >= eta*clow && itstep-lastdw > patMax(5,patMin(20,n/10)) ) {
    
        /* accept new (diminished ) weights */
        
        if ( clow > itstep/10 ) {
            eta = onep3*eta;
            if ( ! silent ) o8info(11);
        }
        lastch = itstep;
        lastdw = itstep;
        level  = diff0/iterma;
        psist  = s1;
        psi    = s2;
        for (i = 1 ; i <= nres ; i++) {
            o8w[i] = w1[i];
        }
        clow = clow+one;
    } else {
    
        /* increase individual weights if necessary. let weigths unchanged */
        /* otherwise                                                       */
        
        s1 = zero;
        s2 = zero;
        for (i = 1 ; i <= nres ; i++) {
            if ( w1[i] > o8w[i] ) {
                lastup = itstep;
                lastch = itstep;
            }
            o8w[i] = patMax(o8w[i],w1[i]);
            if ( i <= nh ) {
                s1 = s1+o8w[i]*myfabs(resst[i]);
                s2 = s2+o8w[i]*myfabs(res[i]);
            } else {
                s1 = s1-o8w[i]*patMin(zero,resst[i]);
                s2 = s2-o8w[i]*patMin(zero,res[i]);
            }
        }
        psist = s1;
        psi   = s2;
    }   
    term = zero;
    if ( nres >= 1 ) term = o8w[1];
    for (i = 2 ; i <= nres ; i++) {
        term = patMax(term,o8w[i]);
    }
    accinf[itstep][20] = term;
    
    /* maximum of weights */
    
    accinf[itstep][19] = clow;
    
    if ( ! silent ) o8info(12);
    
    return;
}

/* **************************************************************************** */
/*          computation of the Pantoja-Mayne BFGS-update of hessian             */
/* **************************************************************************** */
void o8bfgs(void) {

    void    o8msg (INTEGER num);
    void    o8inim(void);
    DOUBLE  o8sc1 (INTEGER i,INTEGER j,DOUBLE a[],DOUBLE b[]);
    DOUBLE  o8sc2 (INTEGER n,INTEGER m,INTEGER j,DOUBLE a[][NX+1],   DOUBLE b[]);
    DOUBLE  o8sc3 (INTEGER n,INTEGER m,INTEGER j,DOUBLE a[][NRESM+1],DOUBLE b[]);
    DOUBLE  o8sc3_(INTEGER n,INTEGER m,INTEGER j,DOUBLE a[][NX+1],   DOUBLE b[]);
    void    o8upd(DOUBLE r[][NX+1],DOUBLE z[],DOUBLE y[],INTEGER n,LOGICAL *fail);
    DOUBLE  o8vecn(INTEGER nl,INTEGER nm,DOUBLE x[]);
    
    static INTEGER  i,j;
    static DOUBLE   dg[NX+1],adx[NX+1],den1,den2,den3,
                    th,tk,xsik,ltdx[NX+1],gtdx[NRESM+1],updx[NX+1],updz[NX+1],
                    term,term1,anorm,acond,ndx,ngtdx,den21;
    static LOGICAL  fail;

    for ( i = 1 ;  i <= n ;  i++) {
    
        /* multiply dx = (s in the usual notation) by Cholesky-factor */
        /* stored in the upper half of a                              */
        
         ltdx[i] = o8sc2(i,n,i,a,difx);
         dg[i]   = gphi1[i]-gphi0[i];
    }
    if ( o8vecn(1,n,dg) == zero ) {
    
        /* suppress update */
        
        accinf[itstep][27] = zero;
        accinf[itstep][28] = zero;
        accinf[itstep][29] = zero;
        
        if ( ! silent ) o8msg(21);
        
        return;
    }
    for (i = 1 ; i <= n ; i++) {
        adx[i] = o8sc3_(1,i,i,a,ltdx);
    }
    /* adx = a * ( x-x0), x-x0 = difx */
    
    for (i = 1 ; i <= DONLP2alist[0] ; i++) {
        gtdx[i] = o8sc3(1,n,DONLP2alist[i],gres,difx);
        gtdx[i] = gtdx[i]/gresn[DONLP2alist[i]];
    }
    /* gtdx = grad(res)(transp)*(x-x0) */
    
    ndx   = o8vecn(1,n,difx);
    tk    = patMin(p5,pow(dnorm,2));
    anorm = zero;
    term1 = myfabs(a[1][1]);
    anorm = zero;
    for (i = 1 ; i <= n ; i++) {
        for (j = i ; j <= n ; j++) {
            anorm = anorm+pow(a[i][j],2);
        }
        term1 = patMin(term1,myfabs(a[i][i]));
    }
    if ( term1 != zero ) {
        acond = anorm/pow(term1,2);
    } else {
        acond = epsmac/tolmac;
    }
    den1 = pow(o8vecn(1,n,ltdx),2);
    den2 = o8sc1(1,n,dg,difx);
    if ( den1 <= rho1*anorm*pow(ndx,2) || acond >= one/rho1 ) {
    
        /* take a restart step */
        
        o8inim();
        
        return;
    }
    if ( nres == 0 ) {
    
        /* in the unconstrained case we take the Powell update */
        
        th = one;
        if ( den2 < p2*den1 ) {
            th = p8*den1/(den1-den2);
            for (i = 1 ; i <= n ; i++) {
                dg[i] = th*dg[i]+(one-th)*adx[i];
            }
            den2 = o8sc1(1,n,dg,difx);
        }
        term = one/sqrt(den2);
        for (i = 1 ; i <= n ; i++) {
            dg[i]   = dg[i]*term;
            updz[i] = dg[i];
        }
        term = one/sqrt(den1);
        for (i = 1 ; i <= n ; i++) {
            updx[i] = adx[i]*term;
        }
        accinf[itstep][28] = den2/den1;
        accinf[itstep][29] = th;
        accinf[itstep][27] = two;
        if ( th != one ) accinf[itstep][27] = three;
    } else {
        ngtdx = o8vecn(1,DONLP2alist[0],gtdx);
        term  = one/sqrt(den1);
        for (i = 1 ; i <= n ; i++) {
            updx[i] = adx[i]*term;
        }
        if ( den2 >= rho1*o8sc1(1,n,dg,dg)
        && o8vecn(1,n,dg) >= sqrt(epsmac)*ndx ) {
            xsik = zero;
            for (i = 1 ; i <= n ; i++) {
                updz[i] = dg[i];
            }
            den21 = den2;
        } else {
        
            /* try Pantoja-Mayne modification */
            
            den3 = tk*pow(ndx,2)+pow(ngtdx,2);
            if ( den2 >= rho1*o8sc1(1,n,dg,dg) ) {
                xsik = one;
            } else {
                xsik = one+(tk*pow(ndx,2)+myfabs(den2))/den3;
            }
            for (i = 1 ; i <= n ; i++) {
                term = zero;
                for (j = 1 ; j <= DONLP2alist[0] ; j++) {
                    term1 = gres[i][DONLP2alist[j]]*gtdx[j];
                    term1 = term1/gresn[DONLP2alist[j]];
                    term  = term+term1;
                }
                updz[i] = dg[i]+xsik*(tk*difx[i]+term);
            }
            den21 = o8sc1(1,n,updz,difx);
        }
        term = one/sqrt(den21);
        for (i = 1 ; i <= n ; i++) {
            updz[i] = updz[i]*term;
        }
        th = one;
        if ( den2 < p2*den1 ) {
            th = p8*den1/(den1-den2);
            for (i = 1 ; i <= n ; i++) {
                dg[i] = th*dg[i]+(one-th)*adx[i];
            }
            den2 = o8sc1(1,n,dg,difx);
        }
        term = one/sqrt(den2);
        for (i = 1 ; i <= n ; i++) {
            dg[i] = dg[i]*term;
        }
        if ( o8vecn(1,n,dg) <= tm3*o8vecn(1,n,updz) ) {
        
            /* the Powell update produces a smaller growth */
            
            for (i = 1 ; i <= n ; i++) {
                updz[i] = dg[i];
            }
            accinf[itstep][28] = den2/den1;
            accinf[itstep][29] = th;
            accinf[itstep][27] = two;
            if ( th != one ) accinf[itstep][27] = three;
        } else {
        
            /* no update if strongly irregular */
            
            accinf[itstep][27] = one;
            accinf[itstep][28] = tk;
            accinf[itstep][29] = xsik;
        }
    }
    o8upd(a,updz,updx,n,&fail);
    
    /* check illconditioning after updating */
    
    term  = myfabs(a[1][1]);
    term1 = term;
    i     = 1;
    
    /* in order to overcome a curious error in hp's f77 compiler */
    /* this kind of loop                                         */
    
    while ( i < n ) {
        i     = i+1;
        term  = patMax(term, myfabs(a[i][i]));
        term1 = patMin(term1,myfabs(a[i][i]));
    }
    if ( fail || pow(term1,2) <= rho1*pow(term,2) ) {
    
        /* reset */
        
        o8inim();
    }
    return;
}

/* **************************************************************************** */
/*                      write short information on standard out                 */
/* **************************************************************************** */
void o8shms(void) {

    static DOUBLE   umin;

    if ( te0 && ! silent ) {
        umin = accinf[itstep][11];
        printf(
        "%5i fx= %14.6e upsi= %8.1e b2n= %8.1e umi= %8.1e nr%4i si%2i\n",
        itstep,fx,upsi,b2n,umin,nr,(int)accinf[itstep][10]);
        
        fprintf(prou,
        "%5i fx= %14.6e upsi= %8.1e b2n= %8.1e umi= %8.1e nr%4i si%2i\n",
        itstep,fx,upsi,b2n,umin,nr,(int)accinf[itstep][10]);
    }
    return;
}

/* **************************************************************************** */
/*                      write messages on "special events"                      */
/* **************************************************************************** */
void o8msg(INTEGER num) {

    static INTEGER  i;

    if ( num <= 0 || num  > 25 ) return;
    
    switch (num) {
    
    case 1:
        fprintf(meu,"step=%i\n", itstep);
        fprintf(meu,"rankdeficiency of grad's of active constr.!\n");
        fprintf(meu,"on the basis of delmin!\n");
        
        return;
        
    case 2:
        fprintf(meu,"step=%i\n", itstep);
        fprintf(meu,"rescaling of objective function scf= %.15e\n",scf);
        
        return;
        
    case 3:
        fprintf(meu,"step=%i\n", itstep);
        fprintf(meu,"rankdeficiency of grad's of active constr.!\n");
        fprintf(meu," del= %.15e\n", del);
        
        return;
        
    case 4:
        fprintf(meu,"step=%i\n", itstep);
        fprintf(meu,"rankdeficiency of grad's of active constr.!\n");
        fprintf(meu," del= %.15e\n", del);
        
        return;
        
    case 5:
        fprintf(meu,"qpterm<0: no dir. of. desc., tauqp max\n");
        
        return;
        
    case 6:
        fprintf(meu,"step=%i\n", itstep);
        fprintf(meu,"second order correction suppressed! \n");
        
        return;
        
    case 7:

        fprintf(meu,"step=%i\n", itstep);
        fprintf(meu,"retry next step with a=id \n");
        
        return;
        
    case 8:
        fprintf(meu,"step=%i\n", itstep);
        fprintf(meu,"some constraint has gradient equal to zero \n");
        fprintf(meu,"resulting d may be no direction of descent\n");
        
        return;
        
    case 9:
        fprintf(meu,"step=%i\n", itstep);
        fprintf(meu,"try regularized SQP with increased weights\n");
        fprintf(meu,"since dnorm small or infeasibility large\n");
        
        return;
        
    case 10:
        fprintf(meu,"step=%i\n", itstep);
        fprintf(meu,"QPsolver did not complete, try d anyway, may fail\n");
        
        return;
        
    case 11:
        fprintf(meu,"step=%i\n", itstep);
        fprintf(meu,"direct. deriv. positive! may be due to inaccurate\n");
        fprintf(meu,"gradients or extreme illconditioning\n");
        fprintf(meu,"dirder=  %.15e\n", dirder);
        
        return;
        
    case 12:
        fprintf(meu,"step=%i\n", itstep);
        fprintf(meu,"call of matdru suppressed, mat too large\n");
        
        return;
        
    case 13:
        fprintf(meu,"step=%i\n", itstep);
        fprintf(meu,"startvalue corrected in order to fit bounds\n");
        
        return;
        
    case 14:
        fprintf(meu,"step=%i\n", itstep);
        fprintf(meu,"try full SQP due to slow progress in x \n");
        
        return;
        
    case 15:
        fprintf(meu,"step=%i\n", itstep);
        fprintf(meu,"try full SQP due to small stepsizes while\n");
        fprintf(meu,"infeasibility large\n");
    case 16:
        fprintf(meu,"step=%i\n", itstep);
        fprintf(meu,"on exit from o8qpdu dir. deriv. positive!\n");
        fprintf(meu,"try increased tauqp\n");
        
        return;
        
    case 17:
        fprintf(meu,"step=%i\n", itstep);
        fprintf(meu,"try regularized SQP with increased weights\n");
        fprintf(meu,"no decrease of weights possible\n");
        fprintf(meu,"and incompatibility large\n");
        
        return;
        
    case 18:
        fprintf(meu,"step=%i\n", itstep);
        fprintf(meu,"try regularized SQP with increased weights\n");
        fprintf(meu,"since no direction of descent has been obtained\n");
        
        return;
        
    case 19:
        fprintf(meu,"step=%i\n", itstep);
        fprintf(meu,"degeneracy in dual QP\n");
        fprintf(meu,"restr. %i to be added\n",iptr);
        fprintf(meu,"new rii= %.15e\n",riitr);
        fprintf(meu,"length of current working set=%i\n",iqtr);
        fprintf(meu,"working set\n");
        for (i = 1 ; i <= iqtr ; i++) {
            fprintf(meu,"%4i ",aitr[i]);
            if ( i % 15 == 0 || i == iqtr ) fprintf(meu,"\n");
        }
        fprintf(meu,"primal feasibility violation is= %.15e\n",sstr);
        
        return;
        
    case 20:
        fprintf(meu,"step=%i\n", itstep);
        fprintf(meu,"dual QP seemingly infeasible \n");
        fprintf(meu,"theoretically impossible\n");
        
        return;
        
    case 21:
        fprintf(meu,"step=%i\n", itstep);
        fprintf(meu,"no update since dg=0\n");
        
        return;
        
    case 22:
        fprintf(meu,"step%i\n", itstep);
        fprintf(meu,"function evaluation returns err=true\n");
        
        return;
        
    case 23:
        fprintf(meu,"step%i\n", itstep);
        fprintf(meu,"contraint evaluation returns err=true\n");
        
        return;
        
    case 24:
        fprintf(meu,"step%i\n", itstep);
        fprintf(meu,"current point cannot be changed in current\n");
        fprintf(meu,"direction due to error-message from function\n");
        fprintf(meu,"evaluation\n");
        
        return;
        
    case 25:
        fprintf(meu,"step%i\n", itstep);
        fprintf(meu,"reduce stepsize due to error-indicator set\n");
        fprintf(meu,"by users function evaluation\n");
        
        return;
    }
}

/* **************************************************************************** */
/* equality constrained recursive quadratic programming with                    */
/* multiple inactivation and superlinearly convergent projected                 */
/* BFGS-update (version 12/93 Spellucci )                                       */
/* **************************************************************************** */
patString o8opti(donlp2Parameters theParameters,patIterationBackup* theBackup) {

    void   o8info (INTEGER icase);
    void    o8sce  (void);
    void    o8bfgs (void);
    void    o8shms (void);
    void    o8msg  (INTEGER num);
    void    o8inim (void);
    void    o8dird (void);
    void    o8cutd (void);
    void    o8smax (void);
    void    o8unim (DOUBLE sig1th);
    void    o8egph (DOUBLE gphi[]);
    void    o8dec  (INTEGER nlow,INTEGER nrl);
    void    o8ht   (INTEGER id,INTEGER incr,INTEGER is1,INTEGER is2,INTEGER m,
                    DOUBLE a[][NRESM+1],DOUBLE beta[],DOUBLE b[],DOUBLE c[]);   
    void    o8sol  (INTEGER nlow,INTEGER nup,DOUBLE b[],DOUBLE x[]);
    void    o8solt (INTEGER nlow,INTEGER nup,DOUBLE b[],DOUBLE x[]);
    void    o8rght (DOUBLE a[][NX+1],DOUBLE b[],DOUBLE y[],DOUBLE *yl,INTEGER n);
    void    o8left (DOUBLE a[][NX+1],DOUBLE b[],DOUBLE y[],DOUBLE *yl,INTEGER n);
    DOUBLE  o8vecn (INTEGER nl,INTEGER nm,DOUBLE x[]);
    void    o8qpdu (void);
    void    esf    (DOUBLE x[],DOUBLE *fx);
    void    esgradf(DOUBLE x[],DOUBLE gradf[]);
    void    esh    (INTEGER i,DOUBLE x[],DOUBLE *hxi);
    void    esgradh(INTEGER i,DOUBLE x[],DOUBLE gradhi[]);
    void    esg    (INTEGER i,DOUBLE x[],DOUBLE *gxi);
    void    esgradg(INTEGER i,DOUBLE x[],DOUBLE gradgi[]);
    
    void    user_eval(DOUBLE xvar[],INTEGER mode);

    static INTEGER  l,l0,i,j,k,csssig,csirup,csreg,cschgx;
    static INTEGER  csmdph;
    static DOUBLE   delsig,delx,sum,term;
    static DOUBLE   umin,term1,scfh,unorm;
    static DOUBLE   del1;
    static DOUBLE   qtx[NX+1];
    static DOUBLE   yy[NX+1],yx[NX+1],trvec[NX+1];
    static INTEGER  iumin,rank0,nr0,csdifx,clwold;
    static INTEGER  delist[NRESM+1],nrbas,bindba[NRESM+1];
    static DOUBLE   eps,delold,uminsc,fac,slackn,tauqp0;
    static LOGICAL  nperm,qpnew,etaini;
    /* added feature 25.01.2000 */
    static LOGICAL  viobnd ;
    /* end added feature */

    /* initialization */
    
    /* save starting point for later printing only */
    
    for (i = 1 ; i <= n ; i++) {
        d[i]  = zero;
        d0[i] = zero;
    }
    itstep    = 0;
    DONLP2alist[0]  = nh;
    delist[0] = 0;
    violis[0] = 0;
    upsi      = zero;
    psi       = zero;
    psi0      = zero;
    sig0      = zero;
    d0norm    = one;
    unorm     = one;
    
    /* in order to have cosphi well defined for itstep = 1 */
    
    dnorm  = one;
    del    = del0;
    
    /* count successive regularization steps */
    
    csreg  = 0;
    
    /* count successive small changes in x */
    
    cschgx = 0;
    
    /* count small differences of fx */
    
    csdifx = 0;
    
    /* count irregular quasi-Newton-updates */
    
    csirup = 0;
    
    /* count successive small stepsizes */
    
    csssig = 0;
    
    /* count successive small differences of penalty-function */
    
    csmdph = 0;
    matsc  = one;
    
    /* formerly tauqp = tp2 is worse */
    
    tauqp  = one;
    nperm  = FALSE;
    ident  = FALSE;
    etaini = FALSE;
    if ( n > 100 || nres > 100 ) te3 = FALSE;
    for (i = 1 ; i <= n ; i ++) {
        perm[i]  = i;
        perm1[i] = i;
    }
    for (i = 1 ; i <= nh ; i ++) {
        bind0[i] = 1;
        bind[i]  = 1;
        DONLP2alist[i] = i;
    }
    if ( analyt ) {
        eps = patMin(epsx,sqrt(epsmac));
    } else {
        eps = epsdif;
        if ( epsx < pow(epsdif,2) ) epsx = pow(epsdif,2);
    }
    eps = patMax(epsmac*tp3,patMin(tm3,eps));

    /* calling for external function evaluation necessary only if corr = TRUE */

    /* function and gradient values, from xtr = x*xsc */
    
    for (i = 1 ; i <= nres ; i ++) {
        if ( i <= nh ) {
            cfuerr[i] = FALSE;
            
            esh(i,x,&res[i]);
            
            if ( cfuerr[i] ) {
                if ( ! silent ) o8msg(23);
                optite = -10;
                
                return patString("Constraint evaluation returns error with current point") ;
            }
            term = myfabs(res[i]);
            if ( ! gconst[i] ) {
            
                /* we assume that the gradient can be evaluated whenever */
                /* the function can                                      */
                
                esgradh(i,x,yy);
                
                val[i] = TRUE;
                for (j = 1 ; j <= n ; j++) {
                    gres[j][i] = yy[j];
                }
            }
        } else {
            bind[i]   = 0;
            cfuerr[i] = FALSE;
            
            esg(i-nh,x,&res[i]);
            
            if ( cfuerr[i] ) {
                if ( ! silent ) o8msg(23);
                optite = -10;
                
                return patString("Constraint evaluation returns error with current point") ;

            }
            term = -patMin(zero,res[i]);
            if ( res[i] <= delmin ) {
                bind[i]  = 1;
                DONLP2alist[0] = DONLP2alist[0]+1;
                DONLP2alist[DONLP2alist[0]] = i;
                if ( ! gconst[i] ) {
                    val[i] = TRUE;
                    
                    /* we assume that the gradient can be evaluated whenever */
                    /* the function can                                      */
                    
                    esgradg(i-nh,x,yy);
                    
                    for (j = 1 ; j <= n ; j++) {
                        gres[j][i] = yy[j];
                    }
                }
            }
        }
        upsi = upsi+term;
        psi  = psi+term*o8w[i];
        if ( val[i] && ! gconst[i] ) gresn[i] = patMax(one,o8vecn(1,n,yy));
    }
    L100:

    /* ******************************************** */
    /* obtaining a point feasible within tau0 first */
    /* ******************************************** */
    
    if ( upsi >= tau0 ) {
        scf   = zero;
        phase = -1;
    } else {
        ffuerr = FALSE;
        
        esf(x,&fx);
        
        if ( ffuerr ) {
            if ( ! silent ) o8msg(22);
            optite = -9;
            
                return patString("Objective evaluation returns error with current point") ;
        }
        if ( ! val[0] ) {
        
            /* we assume that the gradient evaluation can be done whenever */
            /* the function has been evaluated                             */
            
            esgradf(x,gradf);
            
            val[0] = TRUE;
        }
        scf    = one;
        phase  = 0;
        fxst   = fx;
        psist  = psi;
        upsist = upsi;
        for (j = 1 ; j <= nres ; j++) {
            resst[j] = res[j];
        }
        eta = zero;
    }
    L200:

    /* *************************************** */
    /* main iteration loop: getting a better x */
    /* *************************************** */
    
    if ( ! ident ) {
      GENERAL_MESSAGE(itstep << '\t' << fx) ;
      if (theBackup != NULL) {
	theBackup->saveCurrentIteration() ;
      }
        if (patFileExists()(theParameters.stopFileName)) {
	  stringstream str ;
	  str << "Iterations interrupted by the user with the file " 
	      << theParameters.stopFileName ;
	  WARNING(str.str()) ;
	  optite = -three;
	  itstep = iterma;
	  b2n    = accinf[itstep][8];
	  b2n0   = accinf[itstep][7];

	  return patString(str.str());
        }
        itstep = itstep+1;
        if ( itstep > iterma ) {
            optite = -three;
            itstep = iterma;
            b2n    = accinf[itstep][8];
            b2n0   = accinf[itstep][7];
            stringstream str ;
	    str << "Maximum number of iterations (" << iterma << ") reached" ;
            return patString(str.str());
        }
        qpnew  = FALSE;
        qpterm = 0;
        delold = del;
        del    = zero;
        b2n0   = -one;
        b2n    = -one;
        singul = FALSE;
        nperm  = FALSE;
        for (i = 1 ; i <= n ; i++) {
            nperm   = nperm || (perm[i] != perm1[i]);
            perm[i] = perm1[i];
            diag[i] = zero;
        }
    }
    /* ********************************************************* */
    /* current valid row permutation for QR-decomposition of     */
    /* matrix of binding gradients in order to obtain continuity */
    /* of the QR-decomposition is given by perm                  */
    /* ********************************************************* */
    
    nr    = DONLP2alist[0];
    nrbas = nr;
    for (i = 1 ; i <= nres ; i++) {
        bindba[i] = bind[i];
    }
    for (j = 1 ; j <= 32 ; j++) {
        accinf[itstep][j] = zero;
    }
    gfn = o8vecn(1,n,gradf);
    
    /* compute new weight of objective function if useful */
    
    if ( nres > 0 && phase >= 0 && ! ident && itstep > 1 && 
        ( (accinf[itstep-1][10] == -1. && scf0 == one )
        || accinf[itstep-1][10] == 1. ) ) {
        
        /* try rescaling the objective function */
    
        term = zero;
        for (i = 1 ; i <= nres ; i++) {
            if ( gunit[1][i] != 1 ) term = patMax(term,gresn[i]);
        }
        scfh = term/patMax(one/scfmax,gfn);
        if ( scfh < one/scfmax ) scfh = one/scfmax;
        if ( scfh > scfmax     ) scfh = scfmax;
        if ( ((fxst-fx)*scfh+scfh/scf*(psist-psi) >= 
            scfh/scf*eta*clow && lastch <= itstep-4 
	      && scfh < tm1*scf) || scfh > tp1*scf ){
            
            /* rescale the objective function if this seems promising and the */
            /* change is significant                                          */
    
            clow  = clow+1;
            term  = scfh/scf;
            psi   = psi*term;
            psist = psist*term;
            for (i = 1 ; i <= nres ; i++) {
                u[i] = u[i]*term;
            }
            unorm  = unorm*term;
            scf    = scfh;
            lastch = itstep;
            term   = sqrt(term);
            for (i = 1 ; i <= n ; i++) {
                diag0[i] = term*diag0[i];
                for (j = 1 ; j <= n ; j++) {
                    a[j][i] = a[j][i]*term;
                }
            }
            matsc = matsc*term;
            if ( ! silent ) o8msg(2);
        }
    }
    accinf[itstep][1] = itstep;
    accinf[itstep][2] = fx;
    accinf[itstep][3] = scf;
    accinf[itstep][4] = psi;
    accinf[itstep][5] = upsi;

    if ( ! silent ) o8info(1);

    /* begin solver */
    
    /* *********************************************** */
    /* QR-decomposition of matrix of binding gradients */
    /* *********************************************** */
    
    if ( nr >= 1 ) {
        o8dec(1,nr);
    } else {
        donlp2rank = 0;
    }
    o8left(a,gradf,yy,&term,n);
    
    for (i = 1 ; i <= n ; i++) {
        qgf[i] = yy[perm[i]];
    }
    o8ht(1,0,1,donlp2rank,n,qr,betaq,qgf,trvec);
    
    for (i = 1 ; i <= n ; i++) {
        qgf[i] = trvec[i];
    }
    if ( donlp2rank != nr && ! silent ) o8msg(1);
    
    /* ******************************************************* */
    /* compute del as function of x (forcing infeasibility and */
    /* the projected gradient to zero)                         */
    /* ******************************************************* */
    
    b2n0 = o8vecn(donlp2rank+1,n,qgf);
    
    sum  = zero;
    for (i = 1 ; i <= nres ; i++) {
        if ( i <= nh ) {
            sum = sum+myfabs(res[i])/gresn[i];
        } else {
            sum = sum-patMin(zero,res[i])/gresn[i];
        }
    }
    if ( itstep > 1 && accinf[itstep-1][8] >= zero 
        && ! etaini && accinf[itstep-1][18] >= 0 ) {
        etaini = TRUE;
        eta    = (accinf[itstep-1][8]/patMax(one,gfn)+sum
                 +patMin(one,slackn)+patMin(one,myfabs(uminsc)))/patMin(30*n,iterma);
        level  = eta;
    }
    delx = delmin;
    term = scf*(fx0-fx)+psi0-psi;
    if ( term > zero && scf != zero ) delx = patMax(delx,exp(p7*p7*log(term)));
    if ( scf == zero ) delx = patMin(del0*tm4,patMax(delx,upsi*tm2));
    delsig = delmin;
    
    /* del should be large enough to include constraints hit in */
    /* step before violis comes from unimin                     */
    
    for (i = 1 ; i <= violis[0] ; i++) {
        j      = violis[i];
        delsig = patMax(delsig,res[j]/gresn[j]/delfac[j]*(one+tm1));
    }
    del = patMin(del0,patMax(patMin(delsig,five*delx),delx));
    if ( violis[0] == 0 ) del = patMin(del,del01);
    
    /* ********************************************* */
    /* if phase = 2 don't loose a binding constraint */
    /* phase = 2 implies delist[0] = 0               */
    /* ********************************************* */
    
    if ( phase == 2 && violis[0] == 0 ) {
        for (i = nh+1 ; i <= nres ; i++) {
            if ( bind0[i] == 1 ) del = patMin(del01,patMax(del,myfabs(res[i])/gresn[i]));
        }
    }
    /* reduce del by a fixed factor ( tm2) if useful, that is if */
    /* delete-list in the previous step was not empty            */
    
    term = del;
    for (i = 1 ; i <= delist[0] ; i++) {
        j     = delist[i];
        term1 = res[j]/gresn[j]*(one-tm2)/delfac[j];
        if ( term1 >= del*tm2 ) term = patMin(term,term1);
    }
    del = term;
    
    /* if delta becomes too large, we may loose complementary slackness */
    
    if ( itstep > 1 && ! ident && scf != zero ) {
        term = zero;
        for (i = nh+1 ; i <= nres ; i++) {
            term = term+patMax(zero,res[i]/gresn[i]-delmin)
                   *patMax(zero,u[i]-smallw)/gresn[i];
        }
        if ( term > zero ) {
            for (i = nh+1 ; i <= nres ; i++) {
                if ( u[i] > smallw && res[i]/gresn[i] > delmin ) {
                    del = patMax(delmin,
                              patMin(del,res[i]/gresn[i]*(one-tm2)/delfac[i]));
                }
            }
        }
    }
    /* if the current step was singular and not successful,try a greater del */
    /* the same, if qpterm in the last step did signal trouble, but          */
    /* stepsize selection was nevertheless successful                        */
    
    if ( itstep > 1 && accinf[itstep-1][30] < 0. ) del = patMin(tp1*delold,del0);
    
    /* ********************************************* */
    /* include nearly binding inequality constraints */
    /* ********************************************* */
    
    for (i = nh+1 ; i <= nres ; i++) {
        term = res[i]/gresn[i];
        if ( bind[i] == 0 && term <= del*delfac[i] ) {
        
            /* it may be useful to include constraints>0 if near its boundary */
            /* but avoid it, if the constraint was in the old delete-list     */
            
            bind[i]         = 1;
            DONLP2alist[0]        = DONLP2alist[0]+1;
            DONLP2alist[DONLP2alist[0]] = i;
            if ( ! val[i] ) {
                val[i] = TRUE;
                
                esgradg(i-nh,x,yy);
                
                for (j = 1 ; j <= n ; j++) {
                    gres[j][i] = yy[j];
                }
                gresn[i] = patMax(one,o8vecn(1,n,yy));
            }
        }
    }
    rank0 = donlp2rank;
    nr0   = nr;
    nr    = DONLP2alist[0] ;
    
    o8dec(nr0+1,nr);
    
    if ( donlp2rank != nr && ! silent ) o8msg(3);
    
    o8ht(1,0,rank0+1,donlp2rank,n,qr,betaq,qgf,trvec);
    
    for (i = 1 ; i <= n ; i++) {
        qgf[i] = trvec[i];
    }
    for (i = 1 ; i <= n ; i++) {
        yy[i] = -qgf[i]*scf;
    }
    for (i = 1 ; i <= n ; i++) {
        yu[i] = zero;
    }
    /* ********************************************* */
    /* first computation of lagrangian multipliers u */
    /* ********************************************* */
    
    o8sol(1,donlp2rank,yy,yu);
    
    umin  = zero;
    unorm = zero;
    for (i = 1 ; i <= nres ; i++) {
        u[i] = zero;
    }
    iumin  = 0;
    uminsc = zero;
    for (i = 1 ; i <= donlp2rank ; i++) {
        unorm = patMax(unorm,myfabs(yu[i]));
        k     = DONLP2alist[colno[i]];
        u[k]  = -yu[i];
        if ( k > nh ) {
            if ( -yu[i]/gresn[k] < uminsc ) {
                iumin  = k;
                uminsc = -yu[i]/gresn[k];
            }
        }
    }
    if ( scf != zero ) {
        for (i = 1 ; i <= n ; i++) {
            yx[i] = scf*gradf[i];
            for (j = 1 ; j <= nres ; j++) {
                yx[i] = yx[i]-gres[i][j]*u[j];
            }
        }
        b2n = o8vecn(1,n,yx)/scf;
    } else {
        b2n = -one;
    }
    if ( ! silent ) o8info(3);
    
    /* compute new delta */
    
    del1 = del;
    if ( b2n >= zero ) del1 = patMax(del,tm1*patMin(del0,
    exp(p7*log( myfabs(b2n)/(gfn+one) + patMax(zero,-uminsc) + sum)) ));
                
    /* exclude constraints which were candidates for inactivating in the */
    /* previous step if useful                                           */
    
    for (i = 1 ; i <= delist[0] ; i++) {
        j     = delist[i];
        term1 = res[j]/gresn[j]*(one-tm2)/delfac[j];
        if ( term1 >= del1*tm2 ) del1 = patMax(delmin,patMin(del1,term1));
    }
    slackn = zero;
    for (i = nh+1 ; i <= nres ; i++) {
        slackn = slackn+
                 patMax(zero,res[i]/gresn[i]-delmin)*patMax(zero,u[i]-smallw)/gresn[i];
    }
    if ( upsi <= delmin && b2n <= epsx*(gfn+one)
        && b2n != -one  && uminsc >= -smallw &&
        slackn <= delmin*smallw*nres ) {
        
        /* sufficient accuracy in Kuhn-Tucker conditions */
    
        optite = zero;
        
        return patString("Convergence reached") ;
    }
    /* include additional constraints if necessary */
    
    l0 = DONLP2alist[0];
    for (i = 1 ; i <= nres ; i++) {
        term = res[i]/gresn[i]/delfac[i];
        if ( term > del && term <= del1 && bind[i] == 0) {
            bind[i]         = 1;
            DONLP2alist[0]        = DONLP2alist[0]+1;
            DONLP2alist[DONLP2alist[0]] = i;
            if ( ! val[i] ) {
                val[i] = TRUE;
                
                esgradg(i-nh,x,yy);
                
                for (j = 1 ; j <= n ; j++) {
                    gres[j][i] = yy[j];
                }
                gresn[i] = patMax(one,o8vecn(1,n,yy));
            }
        }
    }
    del = del1;
    accinf[itstep][6]  = del;
    accinf[itstep][7]  = b2n0;
    accinf[itstep][9]  = DONLP2alist[0];
    accinf[itstep][10] = -1.;
    nr = DONLP2alist[0];
    if ( l0 != nr ) {
        rank0 = donlp2rank;
        
        o8dec(l0+1,nr);
        
        o8ht(1,0,rank0+1,donlp2rank,n,qr,betaq,qgf,trvec);
        
        for (i = 1 ; i <= n ; i++) {
            qgf[i] = trvec[i];
        }
    }
    if ( ! silent ) o8info(2);
    if ( donlp2rank != nr ) {
        if ( ! silent ) o8msg(4);
        
        goto L400;
    }
    /* ******************************************************** */
    /* second solution for multipliers, rank may have changed ! */
    /* ******************************************************** */
    
    for (i = 1 ; i <= n ; i++) {
        yy[i] = -qgf[i]*scf;
    }
    for (i = 1 ; i <= n ; i++) {
        yu[i] = zero;
    }
    o8sol(1,donlp2rank,yy,yu);
    
    /* remember the column interchanges in qr! yu[i] corresponds to */
    /* u[DONLP2alist[colno[i]]]                                           */
    
    umin  = zero;
    unorm = zero;
    for (i = 1 ; i <= nres ; i++) {
        u[i] = zero;
    }
    iumin  = 0;
    uminsc = zero;
    for (i = 1 ; i <= donlp2rank ; i++) {
        unorm = patMax(unorm,myfabs(yu[i]));
        k     = DONLP2alist[colno[i]];
        u[k]  = -yu[i];
        if ( k > nh ) {
            umin = patMin(umin,-yu[i]);
            if ( -yu[i]/gresn[k] < uminsc ) {
                iumin  = k;
                uminsc = -yu[i]/gresn[k];
            }
        }
    }
    if ( scf != zero ) {
        for (i = 1 ; i <= n ; i++) {
            yx[i] = scf*gradf[i];
            for (j = 1 ; j <= nres ; j++) {
                yx[i] = yx[i]-gres[i][j]*u[j];
            }
        }
        b2n = o8vecn(1,n,yx)/scf;
    }
    accinf[itstep][8]  = b2n;
    accinf[itstep][11] = umin;
    
    o8shms();
    
    if ( ! silent ) o8info(4);
    
    delist[0] = 0;
    if ( phase >= 0 && b2n != -one ) {
        if ( myfabs(uminsc) >= patMax(smallw,myfabs(b2n)/(gfn+one)*c1d)) {
            for (i = nh+1 ; i <= nr ; i++) {
                k = DONLP2alist[colno[i]];
                if ( -yu[i]/gresn[k] <= -smallw) {
                    delist[0]         = delist[0]+1;
                    delist[delist[0]] = k;
                }
            }
        }
    }
    /* *********************************************************** */
    /* the new delist doesn't influence the current d but only the */
    /* computation of the next del                                 */
    /* *********************************************************** */
    
    eqres = TRUE;
    for (i = 1 ; i <= nres ; i ++) {
        eqres = eqres && ( bind[i] == bind0[i] );
    }
    /* compute condition number estimators of diag-r and diag of */
    /* Cholesky-decomposition of b                               */
    
    if ( nr > 1 ) {
        term  = zero;
        term1 = one;
        for (i = 1 ; i <= nr ; i++) {
            term  = patMax(term, myfabs(diag[i]));
            term1 = patMin(term1,myfabs(diag[i]));
        }
        accinf[itstep][13] = term/term1;
    } else if ( nr == 1 ) {
        accinf[itstep][13] = 1.;
    } else {
        accinf[itstep][13] = -1.;
    }
    term  = myfabs(a[1][1]);
    term1 = myfabs(a[1][1]);
    i     = 2;
    while ( i <= n ) {
        term  = patMax(term, myfabs(a[i][i]));
        term1 = patMin(term1,myfabs(a[i][i]));
        i     = i+1;
    }
    accinf[itstep][14] = pow(term/term1,2);
    
    if ( ! silent ) o8info(5);
    
    /* since a represents the Cholesky-factor, this square */
    
    slackn = zero;
    for (i = nh+1 ; i <= nres ; i++) {
        slackn = slackn+
                 patMax(zero,res[i]/gresn[i]-delmin)*patMax(zero,u[i]-smallw)/gresn[i];
    }
    if ( umin >= -smallw &&
        slackn <= delmin*smallw*nres &&
        upsi <= nres*delmin && upsi0 <= nres*delmin
        && myfabs(fx-fx0) <= eps*(myfabs(fx)+one) &&
        b2n != -one &&
        b2n <= tp2*epsx*(gfn+one) ) {
        
        csdifx = csdifx+1;
    } else {
        csdifx = 0;
    }
    if ( phase >= 0 && (accinf[itstep][14] > tp3 || ! analyt ) && csdifx > n ) {
        optite = four;
            
        /* to avoid possible slow convergence with singular    */
        /* projected hessian or inaccurate numerical gradients */
            
        return patString("Very slow primal progress, singular or illconditoned ");
    }
    /* compute damping factor for tangential component if upsi>tau0/2 */

    scf0 = one;
    if ( phase >= 0 && upsi > tau0*p5 )
         scf0 = patMax(one/scfmax,(two*(tau0-upsi)/tau0)*upsi*tm1/patMax(one,gfn) )/scf;
    accinf[itstep][15] = scf0;

    /* **************************** */
    /* compute tangential component */
    /* **************************** */
    
    for (i = nr+1 ; i <= n ; i++) {
        qtx[i] = yy[i]*scf0;
    }
    /* qtx[nr+1],..,qtx[n] is s2 */
    
    /* **************************************************** */
    /* compute right hand side and vertical component       */
    /* use damping for inactivation direction if very large */
    /* no indirect inactivation if infeasibility large and  */
    /* we are not almost stationary on the current manifold */
    /* **************************************************** */
    
    fac = one;
    if ( -umin*c1d > b2n+upsi && b2n != -one ) fac = c1d;
    if ( upsi > tau0*p5 ) fac = zero;
    for (i = 1 ; i <= nr ; i++) {
        k    = DONLP2alist[colno[i]];
        term = res[k];
        if ( k > nh && -yu[i] < zero && term > zero ) term = -term;
        if ( k > nh && -yu[i] < zero ) term = term-yu[i]*fac;
        yx[i] = -term;
    }
    o8solt(1,nr,yx,qtx);
    
    /* qtx is transformed direction of descent for phi */
    
    o8ht(-1,0,1,nr,n,qr,betaq,qtx,yx);
    
    for (i = 1 ; i <= n ; i++) {
        qtx[perm[i]] = yx[i];
    }
    
    /* solve l(transp)*d = qtx, l = a = Cholesky-factor of b */
    
    o8rght(a,qtx,d,&term,n);
    
    /* end solver */
    
    /* compute new penalty weights : regular case */
    
    clwold = clow;
    if ( phase >= 0 ) o8sce();
    if ( clow > clwold ) {
    
        /* tau_qp depends on the (new) weights */
        
        term = o8w[1];
        for (i = 1 ; i <= nres ; i++) {
            term = patMax(term,o8w[i]);
        }
        tauqp = patMax(one,patMin(tauqp,term));
    }
    /* compute parameter phase and stopping criterion */

    if ( uminsc < -smallw )     phase = patMin(1,phase);
    if ( ! eqres )              phase = patMin(0,phase);
    if ( eqres && upsi < tau0 ) phase = patMax(1,phase);
    
    /* rescale and project d if appropriate */
    
    o8cutd();
    
    /* compute the directional derivative dirder */
    
    o8dird();
    
    /* terminate if correction is small */
    
    if ( dnorm <= epsx*(xnorm+epsx) && upsi <= delmin
        && b2n != -one && uminsc >= -smallw && b2n <= epsx*(gfn+one) ) {
        
        optite = one;
        
        return patString("Computed correction small, regular case");
    }
    L350:

    /* reenter from the singular case: dirder has been computed already */
    
    accinf[itstep][16] = xnorm;
    accinf[itstep][17] = dnorm;
    accinf[itstep][18] = phase;
    
    /* compute stepsize */
    
    cfincr = icf;

    /* if no descent direction is obtained, check whether restarting the method */
    /* might help                                                               */
    
    if ( dirder >= zero ) {
    
        /* no direction of descent */
        
        if ( ! silent ) o8msg(11);
        stptrm = -two;
        sig    = zero;
        
        goto L360;
    }
    /* if directional derivative correct but very small, terminate */
    /* since no further progress might be possible                 */
    
    if ( -dirder <= epsmac*tp2*(scf*myfabs(fx)+psi+one) ){
        if ( upsi > delmin*nres ) {
            optite = -one;
            stptrm = -one;
        } else {
            optite = two;
            stptrm = one;
        }
        sig = zero;
        
        goto L360;
    }
    /* phase = 2 : we may hope to obtain superlinear convergence */
    /* switch to Maratos-correction is on then                   */
    /* return to phase = 1 if first order correction large       */
    
    if ( phase >= 1 && dnorm <= smalld*(xnorm+smalld) && scf0 == one &&
        uminsc >= -smallw && ! singul) phase = 2;
        
    /* return to phase 1 since correction large again */
    
    if ( phase == 2 && dnorm > (xnorm+smalld) ) phase = 1;

    o8smax();
    
    /* stmaxl is the maximal stepsize such that point on projected */
    /* ray changes with sigma, but sigla at most                   */
    
    for (i = 1 ; i <= n ; i++) {
        dd[i] = zero;
        x1[i] = x[i]+d[i];
    }
    /* compute second order correction of infeasibility if useful */
    viobnd = FALSE ;     
    for ( i = 1 ; i<= n ; i++ ) {
      if ( (llow[i] && x1[i] < ug[i]-taubnd) || (lup[i] && x1[i] > og[i]+taubnd) )
        { viobnd = TRUE ; }
    }

    if ( phase == 2 && dnorm > xnorm*sqrt(epsmac)  && ! singul 
         && ! viobnd && nres > 0 ) {
        if ( bloc ) user_eval(x1,-1);
        
        /* only function values, from xtr = xsc*x1 */
        
        for (i = 1 ; i <= DONLP2alist[0] ; i++) {
            yx[i] = zero;
            if ( i <= nh && ! gconst[DONLP2alist[i]] ) {
                cfuerr[i] = FALSE;
                
                esh(i,x1,&yx[i]);
                
                if ( cfuerr[i] ) goto L355;
            } else {
                if ( ! gconst[DONLP2alist[i]] ) {
                    cfuerr[DONLP2alist[i]] = FALSE;
                
                    esg(DONLP2alist[i]-nh,x1,&yx[i]);
                
                    if ( cfuerr[DONLP2alist[i]] ) goto L355;
                }
            }
            yx[i] = -yx[i];
        }
        for (i = 1 ; i <= DONLP2alist[0] ; i++) {
            yy[i] = yx[colno[i]];
        }
        o8solt(1,nr,yy,dd);
        
        for (i = nr+1 ; i <= n ; i++) {
            dd[i] = zero;
        }
        o8ht(-1,0,1,nr,n,qr,betaq,dd,yx);
        
        for (i = 1 ; i <= n ; i++) {
            dd[perm[i]] = yx[i];
        }
        o8rght(a,dd,dd,&term,n);
        
        if ( sqrt(term) > p5*dnorm ) {
        
            /* second order correction almost as large as first order one: */
            /* not useful                                                  */
            
            for (i = 1 ; i <= n ; i++) {
                dd[i] = zero;
            }
            if ( ! silent ) o8msg(6);
        }
    }
    L355:

    if ( ! silent ) o8info(7);

    sig = patMin(one,stmaxl);
    
    o8unim(sig);

    L360:

    cfincr = icf-cfincr;
    if ( ! silent ) o8info(10);
    
    /* count successive small steps */
    
    term = scf*(fx0-fx)+psi0-psi;
    if ( myfabs(term) <= epsmac*tp3*(scf*myfabs(fx)+psi) ) {
        csmdph = csmdph+1;
    } else {
        csmdph = 0;
    }
    /* csmdph counts contiguous small differences of penalty function phi */
    
    if ( csmdph > patMax(n,10) ) {
        optite = seven;
        
        return patString("Small differences of penalty function");
    }
    if ( sig <= five*tm2 ) {
        if ( sig0 <= five*tm2 ) csssig = csssig+1;
    } else {
        csssig = 0;
    }
    /* csssig counts the number of successive small sig's */

    accinf[itstep][21] = sig;
    accinf[itstep][22] = cfincr;
    accinf[itstep][23] = dirder;
    accinf[itstep][24] = dscal;
    accinf[itstep][25] = cosphi;
    accinf[itstep][26] = violis[0];
    
    if ( sig == zero && stptrm == one && optite == two ) {
    
        /* no further significant progress possible */
        
        if ( ! silent ) o8info(17);
        
        return patString("No further significant progress possible");
    }
    if ( stptrm == one && sig <= tm4 && accinf[itstep][13] > tp4 && ! singul
        && nres > 0) {
        
        /* try a regularized step, hopefully this will give a better d */
        /* and larger sig                                              */
        
        if ( accinf[itstep][14] > tp4 ) o8inim();
        ident  = TRUE;
        singul = TRUE;
        
        goto L400;
    }
    if ( stptrm < zero ) {
    
        /* stepsize selection failed */
        
        if ( ! ident ) {
        
            /* try restart with a = identity scaled */
            
            if ( ! silent ) o8msg(7);
            ident     = TRUE;
            delist[0] = 0;
            violis[0] = 0;
            csreg     = 0;
            csssig    = 0;
            csirup    = 0;
            
            o8inim();
            
            DONLP2alist[0] = nrbas;
            for (i = 1 ; i <= nres ; i++) {
                bind[i] = bindba[i];
            }
            if ( upsi >= tau0 ) {
            
                goto L100;
                
            } else {
            
                goto L200;
            }
        }
        if ( ! singul && ident && accinf[itstep][13] > tp4 && nres > 0 ) {
        
            /* try the full SQP-direction               */
            /* this may be the third try for this point */
            
            singul = TRUE;
            ident  = TRUE;
            
            goto L400;
        }
        if ( stptrm == -two ) {
            optite = -four;
            
            return patString("Fails in finding a descent direction");
        }
        if ( sig == zero && optite == -one ) return patString("Failed to converge");
        
        /* unidimensional search unsuccessfully terminated */
    
        optite = -two;
        
        return patString("Unsuccessful linesearch");
    }
    if ( singul && itstep > n && myfabs(fx-fx0) <= eps*(myfabs(fx)+one) && phase >= 0
        && upsi <= nres*delmin && upsi0 <= nres*delmin
        && slackn <= delmin*smallw*nres && infeas <= upsi && ! ident ) {
        
        /* since multipliers may be incorrect for infeas != zero be careful */
    
        optite = four;
        
        /* avoid slow progress in case of singular constraints */
    
        return patString("Very slow primal progress, singular or illconditoned problem");
    }
    /* relaxed termination criteria in the singular case */
    
    if ( singul && upsi <= delmin*nres && upsi0 <= delmin*nres
         && b2n != -one && b2n <= (gfn+one)*epsx*tp2 && phase >= 0
         && slackn <= delmin*smallw*nres && infeas <= upsi ) {
         
        /* since multipliers may be incorrect for infeas != zero be careful */
        
        optite = three;
        
        return patString("KKT-conditions (relaxed) satisfied, singular point") ;
    }
    k = 0;
    for (i = 1 ; i <= n ; i++) {
        if ( myfabs(difx[i]) >= epsx*(myfabs(x[i])+tm2) ) k = 1;
    }
    if ( k == 0 ) {
        cschgx = cschgx+1;
    } else {
        cschgx = 0;
    }
    if ( cschgx > nreset && singul ) {
    
        /* very slow progress in x in the singular case. terminate */
        
        optite = five;
        
        return patString("Too many tiny corrections");
    }
    /* new value of x has been accepted */
    
    xnorm = o8vecn(1,n,x);
    ident = FALSE;
    
    o8egph(gphi0);
    
    for (i = 0 ; i <= nres ; i++) {
        if ( ! gconst[i] ) val[i] = FALSE;
    }
    if ( bloc ) user_eval(x,2);
    
    /* evaluate gradients only, since function values are already */
    /* valid from unidimensional minimization                     */
    /* argument is xtr = xsc*x                                    */
    
    if ( phase >= 0 && ! gconst[0] ) {
        val[0] = TRUE;
        
        esgradf(x,gradf);
    }
    for (i = 1 ; i <= DONLP2alist[0] ; i++) {
        l = DONLP2alist[i];
        if ( ! val[l] ) {
            val[l] = TRUE;
            if ( l <= nh ) {
            
                esgradh(l,x,yx);
                
            } else {
            
                esgradg(l-nh,x,yx);
            }
            for (j = 1 ; j <= n ; j++) {
                gres[j][l] = yx[j];
            }
            gresn[l] = patMax(one,o8vecn(1,n,yx));
        }
    }
    o8egph(gphi1);
    
    for (i = 1 ; i <= n ; i++) {
        yx[i] = x[i]-x0[i];
        yy[i] = gphi1[i]-gphi0[i];
    }
    /* since a represents the Cholesky-factor, this sqrt */
    
    term = sqrt(o8vecn(1,n,yy)/o8vecn(1,n,yx));
    if ( term != zero && phase >= 0 ) matsc = patMax(one/scfmax,patMin(scfmax,term/2));
    
    /* current scaling of identity in case of restart */

    DONLP2alist[0] = 0;
    for (i = 1 ; i <= nres ; i++) {
        u0[i]    = u[i];
        bind0[i] = bind[i];
        res0[i]  = res[i];
        if ( i <= nh ) {
            bind[i] = 1;
        } else {
            bind[i] = 0;
            if ( res[i]/gresn[i] <= delmin ) bind[i] = 1;
            if ( ! val[i] && bind[i] == 1 ) {
                val[i] = TRUE;
                
                esgradg(i-nh,x,yx);
                
                for (j = 1 ; j <= n ; j++) {
                    gres[j][i] = yx[j];
                }
                gresn[i] = patMax(one,o8vecn(1,n,yx));
            }
        }
        if ( bind[i] == 1 ) {
            DONLP2alist[0]        = DONLP2alist[0]+1;
            DONLP2alist[DONLP2alist[0]] = i;
        }
    }
    /* **** define evaluation sequence for constraints:                         */
    /* **** nonbinding inequalities in ascending order of value                 */
    /* **** before the binding ones (commented here, too expensive for ng large)*/
    /* sort k = 0;                                                              */
    /* sort for (i = 1 ; i <= nres ; i++) {                                     */
    /* sort     if ( bind[i] == 0 ) {                                           */
    /* sort         k       = k+1;                                              */
    /* sort         DONLP2sort[k] = i;                                                */
    /* sort         work[k] = res[i];                                           */
    /* sort     }                                                               */
    /* sort }                                                                   */
    /* sort for (i = 1 ; i <= DONLP2alist[0] ; i++) {                                 */
    /* sort     k       = k+1;                                                  */
    /* sort     DONLP2sort[k] = DONLP2alist[i];                                             */
    /* sort     work[k] = res[DONLP2alist[i]];                                        */
    /* sort }                                                                   */
    /* sort for (i = 1 ; i <= nres-DONLP2alist[0] ; i++) {                            */
    /* sort     term = work[i];                                                 */
    /* sort     j    = i;                                                       */
    /* sort     for (k = i+1 ; k <= nres-DONLP2alist[0] ; k++) {                      */
    /* sort         if ( work[k] < term ) {                                     */
    /* sort             term = work[k];                                         */
    /* sort             j    = k;                                               */
    /* sort         }                                                           */
    /* sort     }                                                               */
    /* sort     k       = DONLP2sort[i];                                              */
    /* sort     DONLP2sort[i] = DONLP2sort[j];                                              */
    /* sort     DONLP2sort[j] = k;                                                    */
    /* sort     term    = work[i];                                              */
    /* sort     work[i] = work[j];                                              */
    /* sort     work[j] = term;                                                 */
    /* sort }                                                                   */
    /* **** bind now corresponds to the state of the new point                  */
    /* **** but there may be gradients evaluated at the new point               */
    /* **** not yet put to bind                                                 */
    /* ****                                                                     */
    /* **** update the unprojected quasi-Newton-matrix anyway                   */
    /* ****                                                                     */

    if ( scf != zero ) {
        if ( csirup > nreset || csssig > nreset || csreg  > nreset ) {
            csreg  = 0;
            csssig = 0;
            csirup = 0;
         
            o8inim();
            
        } else {
        
            o8bfgs();
         
            /* for projected update:  if ( ! silent ) o8info(13) */
            /* for Pantoja&Mayne update:                         */
        
            if ( ! silent ) o8info(14);
        }
    }
    /* proceed */

    if ( accinf[itstep][27] == one ) {
        if ( itstep > 1 && 
            accinf[itstep-1][29] != zero && accinf[itstep][29] != zero ) {

            /* c ount s uccessive ir regular up dates */
            
            csirup = csirup+1;
        } else {
            csirup = 0;
        } 
    }
    /* accinf[itstep][27] = 1   update Pantoja&Mayne                   */
    /*                    = 0   noupdate                               */
    /*                    = -1  restart                                */
    /*                    = 2   normal BFGS (nr = 0)                   */
    /*                    = 3   Powell's modified BFGS (nr = 0)        */
    /* accinf[itstep][28] = modification term tk/ den2/den1  resp.     */
    /* accinf[itstep][29] = modification factor th/xsik resp.          */
    /* csirup counts the number of successive irregular updating steps */
    
    if ( phase == -1 ) {
    
        goto L100;
        
    } else {
    
        goto L200;
    }
    L400:

    singul = TRUE;
    phase  = patMin(phase,0);
    accinf[itstep][10] = one;
    
    /* try to compute a descent direction using                       */
    /* an extended quadratic program with                             */
    /* individual slack variable for any constraint                   */
    /* compute damping factor for tangential component if upsi>tau0/2 */
    /* by rescaling f if possible                                     */
    
    scf0 = one;
    if ( phase >= 0 && upsi > tau0*p5 ) {
        scfh = patMax(one/scfmax,
               patMin(scfmax,(two*(tau0-upsi)/tau0)*upsi*tau/patMax(one,gfn) ) );
        if ( ((fxst-fx)*scfh+scfh/scf*(psist-psi) >= 
             scfh/scf*eta*clow && lastch <= itstep-4
		&& scfh < tm1*scf) || scfh > tp1*scf ){
             
            /* rescale the objective function if this seems promising and */
            /* the change is significant                                  */
            
            clow  = clow+1;
            term  = scfh/scf;
            scf0  = term ;
            psi   = psi*term;
            psist = psist*term;
            for (i = 1 ; i <= nres ; i++) {
                u[i] = u[i]*term;
            }
            unorm  = unorm*term;
            scf    = scfh;
            lastch = itstep;
            accinf[itstep][15] = scf;
            term = sqrt(term);
            for ( i = 1 ; i <= n ; i++) {
                diag0[i] = term*diag0[i];
                for ( j = 1 ; j <= n ; j++) {
                    a[j][i] = a[j][i]*term;
                }
            }
            matsc = matsc*term;
            if ( ! silent ) o8msg(2);
        }
    }
    /* slack is upsi at most */
    
    accinf[itstep][32] = upsi;
    if ( ! silent ) o8info(15);

    accinf[itstep][13] = -one;
    term  = myfabs(a[1][1]);
    term1 = term;
    i     = 2;
    while ( i <= n ) {
        term  = patMax(term ,myfabs(a[i][i]));
        term1 = patMin(term1,myfabs(a[i][i]));
        i     = i+1;
    }
    accinf[itstep][14] = pow(term/term1,2);
    if ( ! silent ) o8info(5);
    clwold = clow;
    
    /* save for restart */
    
    tauqp0 = tauqp;
    for (i = 1 ; i <= nres ; i++) {
        u[i] = zero;
    }
    for (i = 1 ; i <= n ; i++) {
        dd[i] = zero;
    }
    o8qpdu();
    
    if ( dnorm == zero && qpterm == 1 && optite == three ) return patString("KKT-conditions (relaxed) satisfied, singular point");

    if ( dnorm <= epsx*(patMin(xnorm,one)+epsx) && qpterm < 0 ) {
    
        /* may be it failed because of illconditioning */
    
        if ( upsi >= nres*delmin && qpnew ) {
        
            /* restarting the method has been done already: game is over */
            
            optite = -one;
            if ( ! silent ) o8info(18);
            
            return patString("Small correction from QP, infeasible point");
        }
        if ( qpnew ) {
            optite = qpterm-five;
            if ( ! silent ) o8info(18);
            
            return patString("Error while solving QP");
        }
        /* try a = id */
        
        qpnew = TRUE;
        for (i = 1 ; i <= nres ; i++) {
            o8w[i] = one;
        }
        lastch = itstep;
        
        o8inim();
            
        ident = TRUE;
        
        /* formerly tauqp = tauqp0 */
        
        tauqp    = one;
        DONLP2alist[0] = nrbas;
        for (i = 1 ; i <= nres ; i++) {
            bind[i] = bindba[i];
        }
        if ( scf == zero ) {
        
            goto L100;
            
        } else {
        
            goto L200;
        }
    }
    accinf[itstep][11] = zero;
    delist[0] = 0;
    umin      = zero;
    
    /* b2n is defined also internally in o8qpdu */
    
    o8shms();
    
    if ( (qpterm >= 0 || qpterm == -3) && scf != zero ) {
        unorm = myfabs(u[1]);
        for (i = 2 ; i <= nres ; i++) {
            unorm = patMax(unorm,myfabs(u[i]));
        }
        for (i = 1 ; i <= n ; i++) {
            yx[i] = scf*gradf[i];
            for (j = 1 ; j <= nres ; j++) {
                yx[i] = yx[i]-gres[i][j]*u[j];
            }
        }
        b2n = o8vecn(1,n,yx)/scf;
    } else {
        b2n = -one;
        
        /* signals "undefined" here */
    }
    if ( ! silent ) o8info(2);
    if ( ! silent ) o8info(16);
    
    /* to avoid termination */
    
    if ( b2n == -one ) b2n = epsmac/tolmac;
    if ( qpterm >= 0 && dnorm <= tm2*epsx*(epsx+patMin(one,xnorm))) {
      patString status ;
        if ( upsi <= nres*delmin ) {
            optite = six;
	    status = patString("Correction from QP very small, almost feasible, singular") ;
        } else {
            optite = -five;
	    status = patString("QP solver: on exit correction small, infeasible point") ;
        }
        return status;
    }
    /* check whether QPsolver terminated unsuccessfully */
    
    if ( qpterm < 0 ) {
    
        /* we have a unfeasible solution for QP. try it */
        /* but don't complain if it fails               */
        
        if ( ! silent ) o8msg(10);
    }
    if ( clow > clwold ) {
        term = one;
        for (i = 1 ; i <= nres ; i++) {
            term = patMax(term,o8w[i]);
        }
        tauqp = patMax(one,term);
    }
    if ( tauqp > pow(taufac,3)*tauqp0 ) tauqp = pow(taufac,3)*tauqp0;
    
    /* no change of tauqp otherwise */
    
    b2n0   = b2n;
    umin   = zero;
    uminsc = zero;
    if ( qpterm >= 1 ) {
        slackn = zero;
        for (i = nh+1 ; i <= nres ; i++) {
            slackn = slackn+
            patMax(zero,res[i]/gresn[i]-delmin)*patMax(zero,u[i]-smallw)/gresn[i];
        }
    } else {
    
        /* slack is undefined, since multipliers are undefined */
        /* use this value to prevent premature termination     */
        
        slackn = one;
    }
    goto L350;
}

/* **************************************************************************** */
/*      initialize the quasi Newton update with a multiple of the identity      */
/* **************************************************************************** */
void o8inim(void) {

    void o8info(INTEGER icase);
    
    static INTEGER i,j;

    for (i = 1 ; i <= n ; i++) {
        for (j = 1 ; j <= n ; j++) {
            a[j][i] = zero;
        }
        a[i][i]  = matsc;
        diag0[i] = matsc;
    }
    accinf[itstep][27] = -one;
    accinf[itstep][14] = one;
    if ( ! silent ) o8info(20);
    
    return;
}

/* **************************************************************************** */
/*          compute the directional derivative of the L1-penalty-function       */
/* **************************************************************************** */
void o8dird(void) {
    
    DOUBLE o8sc1(INTEGER i,INTEGER j,DOUBLE a[],DOUBLE b[]);
    DOUBLE o8sc3(INTEGER n,INTEGER m,INTEGER j,DOUBLE a[][NRESM+1],DOUBLE b[]);
    
    static INTEGER  i;
    static DOUBLE   term,term1;

    /* compute directional derivative of Zangwill function */
    
    //    L100:

    dirder = o8sc1(1,n,gradf,d)*scf;
    
    for (i = 1 ; i <= nres ; i++) {
        term  = o8sc3(1,n,i,gres,d)*o8w[i];
        term1 = res[i];
        if ( i <= nh ) {
            if ( term1/gresn[i] <= -tp3*epsmac ) {
                dirder = dirder-term;
            } else {
                if ( term1/gresn[i] > tp3*epsmac ) {
                    dirder = dirder+term;
                } else {
                    dirder = dirder+myfabs(term);
                }
            }
        } else {
            if ( bind[i] == 1 ) {
                if ( myfabs(term1)/gresn[i] <= tp3*epsmac ) {
                    dirder = dirder-patMin(zero,term);
                } else {
                    if ( term1/gresn[i] < -tp3*epsmac ) {
                        if ( term > zero ) term = patMin(term,-res[i]*o8w[i]);
                        
                        /* only negative values of the constraints */
                        /* contribute to the Zangwill function     */
                        
                        dirder = dirder-term;
                    }
                }
            }
        }
    }
    return;
}

/* **************************************************************************** */
/*                      cut d if appropriate and rescale                        */
/* **************************************************************************** */
void o8cutd(void) {

    DOUBLE o8sc1 (INTEGER i,INTEGER j,DOUBLE a[],DOUBLE b[]);
    DOUBLE o8vecn(INTEGER nl,INTEGER nm,DOUBLE x[]);

    static INTEGER  i;
    static DOUBLE   term,term1;

    xnorm  = o8vecn(1,n,x);
    term   = beta*(xnorm+one);
    dnorm  = o8vecn(1,n,d);
    d0norm = o8vecn(1,n,d0);
    dscal  = one;
    if ( dnorm*d0norm != zero ) {
        cosphi = o8sc1(1,n,d,d0)/(d0norm*dnorm);
    } else {
        cosphi = zero;
    }
    if ( dnorm > term ) {
    
        /* d too long: rescale */
        
        term1 = term/dnorm;
        dnorm = term;
        dscal = term1;
        for (i = 1 ; i <= n ; i++) {
            d[i]  = d[i]*term1;
            dd[i] = dd[i]*pow(term1,2);
        }
    }
    /* since we project the ray with respect to the bounds, be sure */
    /* to compute the directional derivative correctly              */
    /* therefore correct d and dd appropriately                     */
    
    for (i = 1 ; i <= n ; i++) {
        if ( llow[i] && x[i]+sigsm*d[i] <= ug[i] ) {
            d[i]  = zero;
            dd[i] = patMax(zero,dd[i]);
        }
        if ( lup[i] && x[i]+sigsm*d[i] >= og[i] ) {
            d[i]  = zero;
            dd[i] = patMin(zero,dd[i]);
        }
    }
    dnorm = o8vecn(1,n,d);
    
    return;
}

/* **************************************************************************** */
/* compute maximum stepsize stmaxl such that projection on the box of lower     */
/* and upper bounds changes for sig in [0,stmaxl],if such exists                */
/* **************************************************************************** */
void o8smax(void) {
    
    static INTEGER  i;
    static LOGICAL  exis;

    exis = TRUE;
    for (i = 1 ; i <= n ; i++) {
        exis = exis &&( ( d[i] == zero )
            || ( lup[i]  && d[i] > zero )
            || ( llow[i] && d[i] < zero ) );
    }
    if ( exis ) {
        stmaxl = sigsm;
        for (i = 1 ; i <= n ; i++) {
            if ( llow[i] && d[i] < zero ) {
                if ( -d[i]*sigla >= x[i]-ug[i] ) {
                    stmaxl = patMax(stmaxl,(x[i]-ug[i])/(-d[i]));
                } else {
                    stmaxl = sigla;
                }
            }
            if ( lup[i] && d[i] > zero ) {
                if ( d[i]*sigla >= og[i]-x[i] ) {
                    stmaxl = patMax(stmaxl,(og[i]-x[i])/d[i]);
                } else {
                    stmaxl = sigla;
                }
            }
        }
    } else {
        stmaxl = sigla;
    }
    /* but never use stepsize larger than sigla */
    
    stmaxl = patMin(sigla,stmaxl);
    
    return;
}

/* **************************************************************************** */
/*        restore the best point found so far to be the current new point       */
/* **************************************************************************** */
void o8rest(void) {
    
    static INTEGER  j;

    phi1  = phimin;
    psi1  = psimin;
    upsi1 = upsim;
    sig   = sigmin;
    fx1   = DONLPfmin;
    for (j = 1 ; j <= n ; j++) {
        x1[j] = xmin[j];
    }
    for (j = 1 ; j <= nres ; j++) {
        res1[j] = resmin[j];
    }
    return;
}

/* **************************************************************************** */
/*             save the best point found so far in ...min variables             */
/* **************************************************************************** */
void o8save(void) {
    
    static INTEGER  i;
    
    phimin = phi1;
    upsim  = upsi1;
    psimin = psi1;
    DONLPfmin   = fx1;
    sigmin = sig;
    for (i = 1 ; i <= n ; i++) {
        xmin[i] = x1[i];
    }
    for (i = 1 ; i <= nres ; i++) {
        resmin[i] = res1[i];
    }
    return;
}

/* **************************************************************************** */
/*                    evaluate the functions at the new point                   */
/* **************************************************************************** */
void o8eval(DOUBLE sigact,DOUBLE *sigres,LOGICAL *reject,LOGICAL *error) {

    void    esf(DOUBLE x[],DOUBLE *fx);
    void    esh(INTEGER i,DOUBLE x[],DOUBLE *hxi);
    void    esg(INTEGER i,DOUBLE x[],DOUBLE *gxi);
    
    void    user_eval(DOUBLE xvar[],INTEGER mode);
    
    static INTEGER  i,j;
    static DOUBLE   term;
    
    sig = sigact;
    for (i = 1 ; i <= n ; i ++) {
        x1[i] = x[i]+sig*(d[i]+sig*dd[i]);
        
        /* project with respect to the box-constraints */
        
        if ( llow[i] ) x1[i] = patMax(x1[i],ug[i]);
        if ( lup[i]  ) x1[i] = patMin(x1[i],og[i]);
    }
    *reject = FALSE;
    *error  = FALSE;
    *sigres = sig;
    upsi1   = zero;
    psi1    = zero;
    if ( bloc ) user_eval(x1,-1);
    
    /* only function values, from xtr = x1*xsc */
    
    for (j = 1 ; j <= nres ; j++) {
        i = DONLP2sort[j];
        if ( i <= nh ) {
            cfuerr[i] = FALSE;
            
            esh(i,x1,&res1[i]);
            
            if ( cfuerr[i] ) {
                *error = TRUE;
                
                return;
            }
            term = myfabs(res1[i]);
        } else {
            cfuerr[i] = FALSE;
            
            esg(i-nh,x1,&res1[i]);
            
            if ( cfuerr[i] ) {
                *error = TRUE;
                
                return;
            }
            term = -patMin(zero,res1[i]);
            if ( res1[i] < -delmin && bind[i] == 0) {
                violis[0]         = violis[0]+1;
                violis[violis[0]] = i;
            }
        }
        /* violis is the list of inequality-contraints currently        */
        /* not binding which have been hit during unidimensional search */
        /* sigres is the smallest zero of secants through constraints   */
        /* which change sign along [x,x+d]                              */

        upsi1 = upsi1+term;
        if ( upsi1 > tau0 && phase != -1 ) {
            *reject = TRUE;
            
            return;
        }
        psi1 = psi1+term*o8w[i];
        if ( res1[i]*res[i] < zero && sig <= one && ( bind[i] == 0 ||
            (bind[i] == 1 && (myfabs(res[i])/gresn[i] >= tp3*epsmac 
            || myfabs(res1[i])/gresn[i] >= tp3*epsmac ) ) ) )
            *sigres = patMin(*sigres,sig*res[i]/(res[i]-res1[i]));
    }
    if ( phase != -1 ) {
        ffuerr = FALSE;
        
        esf(x1,&fx1);
        
        if ( ffuerr ) {
            *error = TRUE;
            
            return;
        }
    } else {
        fx1 = zero;
    }
    phi1 = scf*fx1+psi1;
    
    return;
}

/* **************************************************************************** */
/*          determination of stepsize by an Armijo-like test for descent        */
/* **************************************************************************** */
void o8unim(DOUBLE sig1th) {

    /* sig1th the first proposed stepsize for searching on the arc              */
    /* if sig = one did'nt work                                                 */

    /* n = number of variables                                                  */
    /* x = current point                                                        */
    /* d = direction of descent. dd = second order correction                   */
    /* x,d = input                                                              */
    /* x0,d0 etc. information from previous step                                */

    /* xnorm,dnorm = euclidean length of x and d                                */
    /* stptrm = 1 on success , = -1 or = -2 otherwise                           */
    /* sig    =  computed stepsize                                              */
    /* it is assumed that one is asymptotically optimal                         */
    /* sigsm = smallest acceptable stepsize                                     */
    /* sigla = largest  acceptable stepsize                                     */
    /* alpha = smallest feasible reduction factor for stepsize                  */
    /* delta = multiplier for derivative should be smaller than .25             */
    /* beta  = maximum feasible increase of x-norm for sig = 1                  */
    /* theta = bound for cos(angle(current direction, previous direction))      */
    /*         if overridden, stepsize larger than one is tried                 */

    /* ************************************************************************ */

    void    o8info(INTEGER icase);
    void    o8rest(void);
    void    o8save(void);
    void    o8eval(DOUBLE sigact,DOUBLE *sigres,LOGICAL *reject,LOGICAL *error);

    static INTEGER  i,l,j;
    static DOUBLE   term,maxphi;
    static DOUBLE   sigres ,diff;
    static DOUBLE   step[NSTEP+1];
    static LOGICAL  desc,descre,sminfe,lainfe,reject,error;

    step[1] = p5;
    step[2] = twom2;
    for (i = 3 ; i <= NSTEP ; i++) {
        step[i] = tm1;
    }
    /* projection of d, rescaling and computing dirder has been done already */
    
    l         = 0;
    error     = FALSE;
    phi       = scf*fx+psi;
    sig       = sig1th;
    violis[0] = 0;
    if ( ! silent ) o8info(8);
    
    L100:

    l = l+1;
    if ( l > NSTEP ) {
        if ( error && ! silent ) o8msg(24);
        stptrm = -one;
        sig    = zero;
        
        return;
    }
    /* compute a new x and test for descent */
    
    o8eval(sig,&sigres,&reject,&error);
    
    if ( error ) {
        if ( sig > one ) {
        
            o8rest();
            
            goto L200;
            
        } else {
            if ( ! silent ) o8msg(25);
            sig = step[l]*sig;
            
            goto L100;
        }
    }
    if ( reject ) {
        if ( sig > one ) {
        
            o8rest();
            
            goto L200;
            
        } else {
            sig = step[l]*sig;
            
            goto L100;
        }
    }
    if ( ! silent ) o8info(9);
    
    /* new function value */
    
    if ( sig > one ) {
        if ( phi1 >= phimin ) {
        
            /* phi does'nt decrease further */
            
            o8rest();
            
            goto L200;
            
        } else {
            if ( sig < stmaxl ) {
            
                o8save();
                
                sig = patMin(stmaxl,sig+sig);
                
                goto L100;
                
            } else {
            
                goto L200;
            }
        }
    }
    if ( lastch >= itstep-3 || phase != 2 || singul ) {
    
        /* require monotonic behaviour */
    
        diff = phi-phi1;
    } else {
        maxphi = phi;
        for (j = 1 ; j <= 3 ; j++) {
            maxphi = patMax(scf*accinf[itstep-j][2]+accinf[itstep-j][4],maxphi);
        }
        diff = maxphi-phi1;
    }
    desc   = diff >= patMin(-sig*delta*dirder,level);
    descre = upsi - upsi1 >= sig*pow(delta,2)*upsi/tauqp;
    sminfe = upsi <= tau0*p5 && upsi1 <= tau0;
    lainfe = upsi > tau0*p5;
    
    if ( desc && ( sminfe || ( lainfe && descre ) ) ) {
    
        /* Goldstein-Armijo descent test satisfied */
        
        if ( sig == one && ( (cosphi >= theta && sig0 >= one 
             && (phase+1)*(phase-2) != 0
             && ! singul) || diff >= -sig*delta1*dirder )
             && stmaxl > one && upsi < tau0*p5 ) {
             
            /* 1 >= delta1 >> delta > 0 */
             
            /* try stepsize larger than one */
            /* save the current point as the best one */
            
            o8save();
            
            sig = patMin(stmaxl,sig+sig);
            
            goto L100;
        }
        if ( sig <= one && upsi > tau0*p5 && upsi1 > upsi ) goto L300;
        
        goto L200;
        
    } else {
    
        goto L300;
    }
    L200:

    /* accept new x, save old values */
    
    fx0    = fx;
    fx     = fx1;
    upsi0  = upsi;
    upsi   = upsi1;
    psi0   = psi;
    psi    = psi1;
    stptrm = one;
    sig0   = sig;
    for (i = 1 ; i <= n ; i++) {
        x0[i]   = x[i];
        d0[i]   = d[i];
        x[i]    = x1[i];
        difx[i] = x[i]-x0[i];
    }
    d0norm = dnorm;
    x0norm = xnorm;
    for (i = 1 ; i <= nres ; i++) {
        res[i] = res1[i];
    }
    return;
    
    /* continue reducing sig */
    
    L300:

    if ( sigres < sig ) {
        sig = patMin(p5*sig,patMax(step[l]*sig,sigres));
    } else {
        term = (diff-dirder*sig)*two;
        if ( term > epsmac*(scf*myfabs(fx)+psi) ) {
            sig = patMin(p5*sig,patMax(step[l]*sig,-dirder*pow(sig,2)/term));
        } else {
            sig = step[l]*sig;
        }
    }
    /* test for (i = 1 ; i <= n ; i++) {                            */
    /* test     if ( myfabs(x1[i]-x[i]) >= sigsm*(myfabs(x[i])+tm2)) {  */
    /* test         goto L100;                                      */
    /* test     }                                                   */
    /* test }                                                       */
    
    if ( sig*patMax(one,dnorm) >= sigsm ) goto L100;
    stptrm = -one;
    sig    = zero;
    
    return;
}

/* **************************************************************************** */
/*                scalar product of two vectors or parts of vectors             */
/* **************************************************************************** */
DOUBLE o8sc1(INTEGER i,INTEGER j,DOUBLE a[],DOUBLE b[]) {

    /* mulpily two vectors */
    
    static INTEGER  k;
    static DOUBLE   s;

    if ( i > j ) {
        return (zero);
    } else {
        s = zero;
        for (k = i ; k <= j ; k++) {
            s = s+a[k]*b[k];
        }
        return (s);
    }
}

/* **************************************************************************** */
/*                    multiply row j of matrix a with vector b                  */
/* **************************************************************************** */
DOUBLE o8sc2(INTEGER n,INTEGER m,INTEGER j,DOUBLE a[][NX+1],DOUBLE b[]) {
    
    static DOUBLE   s;
    static INTEGER  i;

    s = zero;
    for (i = n ; i <= m ; i++) {
        s = s+a[j][i]*b[i];
    }
    return (s);
}

/* **************************************************************************** */
/*          multiply column j section (n to m) of matrix a with vector b        */
/* **************************************************************************** */
DOUBLE o8sc3(INTEGER n,INTEGER m,INTEGER j,DOUBLE a[][NRESM+1],DOUBLE b[]) {

    static DOUBLE   s;
    static INTEGER  i;

    s = zero;
    for (i = n ; i <= m ; i++) {
        s = s+a[i][j]*b[i];
    }
    return (s);
}

DOUBLE o8sc3_(INTEGER n,INTEGER m,INTEGER j,DOUBLE a[][NX+1],DOUBLE b[]) {

    static DOUBLE   s;
    static INTEGER  i;

    s = zero;
    for (i = n ; i <= m ; i++) {
        s = s+a[i][j]*b[i];
    }
    return (s);
}

/* **************************************************************************** */
/* subprogram for structured output of a submatrix a[ma][na]                    */
/* on channel lognum in fix or float format with heading "head".                */
/* uses a fixed format string with 70 print columns                             */
/* **************************************************************************** */
void o8mdru(DOUBLE a[][NRESM+1],INTEGER ma,INTEGER na,char head[],
            FILE *lognum,LOGICAL fix) {

    static INTEGER  i,j,jo,ju;
    static INTEGER  anz;

    fprintf(lognum,"\n%40s\n", head);
    anz = 4;
    jo  = 0;
    while ( jo < na ) {
        ju = jo+1;
        jo = patMin(ju+anz-1,na);
        
        fprintf(lognum,"\nrow/column");
        for (j = ju ; j <= jo ; j++) {
            fprintf(lognum,"      %3i      ", j);
            if ( (j-ju+1) % 4 == 0 || j == jo ) fprintf(lognum,"\n");
        }
        for (i = 1 ; i <= ma ; i++) {
            if ( fix ) {
                fprintf(lognum,"   %4i   ",i);
                for (j = ju ; j <= jo ; j++) {
                    fprintf(lognum,"%14.7f ", a[i][j]);
                    if ( (j-ju+1) % 4 == 0 || j == jo ) fprintf(lognum,"\n");
                }
            } else {
                fprintf(lognum,"   %4i   ",i);
                for (j = ju ; j <= jo ; j++) {
                    fprintf(lognum," %13.6e ", a[i][j]);
                    if ( (j-ju+1) % 4 == 0 || j == jo ) fprintf(lognum,"\n");
                }
            }
        }
    }
    return;
}
void o8mdru_(DOUBLE a[][NX+1],INTEGER ma,INTEGER na,char head[],
             FILE *lognum,LOGICAL fix) {

    static INTEGER  i,j,jo,ju;
    static INTEGER  anz;

    fprintf(lognum,"\n%40s\n", head);
    anz = 4;
    jo  = 0;
    while ( jo < na ) {
        ju = jo+1;
        jo = patMin(ju+anz-1,na);
        
        fprintf(lognum,"\nrow/column");
        for (j = ju ; j <= jo ; j++) {
            fprintf(lognum,"      %3i      ", j);
            if ( (j-ju+1) % 4 == 0 || j == jo ) fprintf(lognum,"\n");
        }
        for (i = 1 ; i <= ma ; i++) {
            if ( fix ) {
                fprintf(lognum,"   %4i   ",i);
                for (j = ju ; j <= jo ; j++) {
                    fprintf(lognum,"%14.7f ", a[i][j]);
                    if ( (j-ju+1) % 4 == 0 || j == jo ) fprintf(lognum,"\n");
                }
            } else {
                fprintf(lognum,"   %4i   ",i);
                for (j = ju ; j <= jo ; j++) {
                    fprintf(lognum," %13.6e ", a[i][j]);
                    if ( (j-ju+1) % 4 == 0 || j == jo ) fprintf(lognum,"\n");
                }
            }
        }
    }
    return;
}

/* **************************************************************************** */
/*                        compute gradient of lagrangian                        */
/* **************************************************************************** */
void o8egph(DOUBLE gphi[]) {

    static INTEGER  i,j,l;

    for (i = 1 ; i <= n ; i++) {
         gphi[i] = gradf[i] * scf;
         for ( j = 1 ; j <= nh ; j++) {
             gphi[i] = gphi[i]-u[j]*gres[i][j];
         }
         for (j = nh+1 ; j <= DONLP2alist[0] ; j++) {
             l = DONLP2alist[j];
             if( u[l] > zero ) gphi[i] = gphi[i]-gres[i][l]*u[l];
             
            /* include constraints, whose multipliers are of correct sign only */
         }
    }
    return;
}

/* **************************************************************************** */
/* QR-decomposition of matrix of gradients of binding constraints               */
/* this set may be expanded using multiple calls to o8dec.                      */
/* No exit on singular r-factor here. Information on                            */
/* the decompostion is stored in betaq and in and below the                     */
/* diagonal of qr. r-factor is stored in diag (diagonal ) and                   */
/* above the diagonal of qr. cscal is the column scaling of the                 */
/* original matrix. Column pivoting is done here and is stored in colno         */
/* **************************************************************************** */
void o8dec(INTEGER nlow,INTEGER nrl) {
    
    DOUBLE  o8sc3 (INTEGER n,INTEGER m,INTEGER j,DOUBLE a[][NRESM+1],DOUBLE b[]);
    void    o8ht  (INTEGER id,INTEGER incr,INTEGER is1,INTEGER is2,INTEGER m,
                   DOUBLE a[][NRESM+1],DOUBLE beta[],DOUBLE b[],DOUBLE c[]);    
    void    o8left(DOUBLE a[][NX+1],DOUBLE b[],DOUBLE y[],DOUBLE *yl,INTEGER n);
    DOUBLE  o8vecn(INTEGER nl,INTEGER nm,DOUBLE x[]);
    
    static INTEGER  n1,n2;
    static INTEGER  i,j,k,l,i1,i2,ipiv;
    static DOUBLE   sum,term,dalpha,dbeta,qrii,dl;
    static DOUBLE   qri[NX+1],qri0[NX+1],curle;

    if ( nlow > nrl ) return;
    if ( nlow == 1 ) donlp2rank = 0;
    dl = one/(n+n+n);
    for ( i = nlow ; i <= nrl ; i++) {
        diag[i]  = zero;
        betaq[i] = zero;
        colno[i] = i;
        for (j = 1 ; j <= n ; j++) {
            qri[j] = gres[j][DONLP2alist[i]];
        }
        o8left(a,qri,qri0,&sum,n);
        
        if ( sum == zero ) {
            cscal[i] = one;
            colle[i] = zero;
            for (j = 1 ; j <= n ; j++) {
                qr[j][i] = zero;
            }
        } else {
            for (j = 1 ; j <= n ; j++) {
                qri[j] = qri0[perm[j]];
            }
            term     = one/sqrt(patMax(sum,pow(rho,2)));
            cscal[i] = term;
            if ( nlow > 1 ) {
            
                o8ht(1,0,1,donlp2rank,n,qr,betaq,qri,qri0);
                
                for (j = 1 ; j <= n ; j++) {
                    qri[j] = qri0[j];
                }
            }
            for (j = 1 ; j <= n ; j++) {
                qr[j][i] = qri[j]*term;
            }
            /* colle : length of remaining column squared */
    
            colle[i] = pow(o8vecn(donlp2rank+1,n,qri)*term,2);
        }
    }
    if ( nlow > 1 && donlp2rank < nlow-1 ) {
    
        /* shift zero block to the right */
        
        i1 = nlow-1-donlp2rank;
        i2 = nrl-nlow+1;
        for (i = 1 ; i <= patMin(i1,i2) ; i++) {
            ipiv        = donlp2rank+i;
            k           = nrl-i+1;
            term        = betaq[k];
            betaq[k]    = betaq[ipiv];
            betaq[ipiv] = term;
            j           = colno[k];
            colno[k]    = colno[ipiv];
            colno[ipiv] = j;
            term        = colle[k];
            colle[k]    = colle[ipiv];
            colle[ipiv] = term;
            for (j = 1 ; j <= n ; j++) {
                term        = qr[j][k];
                qr[j][k]    = qr[j][ipiv];
                qr[j][ipiv] = term;
             }
         }
    }
    if ( nlow > 1 ) {
        n1 = donlp2rank+1;
        n2 = n1+nrl-nlow;
    } else {
        n1 = nlow;
        n2 = nrl;
    }
    for ( i = n1 ; i <= n2 ; i++) {
    
        /* search for pivot column */
        
        ipiv  = i;
        curle = colle[i];
        for (j = i+1 ; j <= n2 ; j++) {
            if ( colle[j] > curle ) {
                curle = colle[j];
            }
        }
        for (j = n2 ; j >= i ; j--) {
            if ( colle[j] >= curle/three ) ipiv = j;
        }
        /* interchange columns explicitly, if necessary   */
        /* make interchanges continuous with respect to x */
        
        if ( ipiv != i ) {
            j           = colno[i];
            colno[i]    = colno[ipiv];
            colno[ipiv] = j;
            term        = colle[i];
            colle[i]    = colle[ipiv];
            colle[ipiv] = term;
            for (k = 1 ; k <= n ; k++) {
                term        = qr[k][i];
                qr[k][i]    = qr[k][ipiv];
                qr[k][ipiv] = term;
            }
        }
        sum = zero;
        for (j = i ;  j <= n ;  j++) {
            term   = qr[j][i];
            qri[j] = term;
            sum    = term*term+sum;
        }
        if ( sum <= pow(rho,2) ) {
        
            /* set tiny values to zero */
            
            for (j = i ; j <= n2 ; j++) {
                colle[j] = zero;
                for (k = i ; k <= n ; k++) {
                    qr[k][j] = zero;
                }
            }
            donlp2rank = i-1;
            
            return;
        }
        qrii   = qri[i];
        dalpha = -sqrt(sum);
        if ( myfabs(qrii) <= -dalpha*dl ) {
            term = zero;
            for (j = i+1 ; j <= n ; j++) {
                if ( myfabs(qri[j]) > term ) {
                    term = myfabs(qri[j]);
                    l    = j;
                }
            }
            k        = perm1[i];
            perm1[i] = perm1[l];
            perm1[l] = k;
        }
        if ( qrii < zero ) dalpha = -dalpha;
        dbeta    = one/(sum-qrii*dalpha);
        diag[i]  = dalpha;
        betaq[i] = dbeta;
        qri[i]   = qrii-dalpha;
        qr[i][i] = qri[i];
        donlp2rank     = i;
        i1       = i+1;
        for (j = i1 ; j <= n2 ; j++) {
            sum = dbeta*o8sc3(i,n,j,qr,qri);
            for (k = i ; k <= n ; k++) {
                qr[k][j] = qr[k][j]-sum*qri[k];
            }
            colle[j] = colle[j]-pow(qr[i][j],2);
        }
    }
    return;
}

/* **************************************************************************** */
/*                  application of Householder transformations                  */
/* **************************************************************************** */
void o8ht(INTEGER id,INTEGER incr,INTEGER is1,INTEGER is2,INTEGER m,
          DOUBLE a[][NRESM+1],DOUBLE beta[],DOUBLE b[],DOUBLE c[]) {

    /* application of Householder transformations stored     */
    /* in the lower or strict lower (if incr = 0 or 1 resp.) */
    /* triangle of a and in the vector beta on b giving c.   */
    /* only columns is1 to is2 are used in forward manner    */
    /* if id > 0,backwards otherwise.                        */
    /* rows is1 to m of c are changed only                   */
    
    DOUBLE o8sc3(INTEGER n,INTEGER m,INTEGER j,DOUBLE a[][NRESM+1],DOUBLE b[]);

    static DOUBLE   sum;
    static INTEGER  i,j,k,it;

    for (i = 1 ; i <= m ; i++) {
        c[i] = b[i];
    }
    if(is1 > m)   return;
    if(is2 < is1) return;
    for (i = is1 ; i <= is2 ; i++) {
        it = i;
        if(id < 0) it = is2-it+is1;
        
        /* it = index of transformation */
        
        j   = it+incr;
        sum = beta[it]*o8sc3(j,m,it,a,c);
        for (k = j ; k <= m ; k++) {
            c[k] = c[k]-sum*a[k][it];
        }
    }
    return;
}
    
/* **************************************************************************** */
/* solve triangular system r*x = b, r defined by Householder-QR-                */
/* decomposition decomp (with column scaling)                                   */
/* **************************************************************************** */
void o8sol(INTEGER nlow,INTEGER nup,DOUBLE b[],DOUBLE x[]) {
    
    static DOUBLE   xl[NX+1];
    static DOUBLE   sum;
    static INTEGER  i,j;

    for (i = nup ; i >= nlow ; i--) {
        sum = zero;
        for (j = i+1 ; j <= nup ; j++) {
            sum = sum+qr[i][j]*xl[j];
        }
        xl[i] = (b[i]-sum)/diag[i];
    }
    for (i = nlow ; i <= nup ; i++) {
        x[i] = xl[i]*cscal[colno[i]];
    }
    /* there must follow interchange of x as given by colno */
    /* e.g. xx(colno[i]) = x[i]                             */

    return;
}

/* **************************************************************************** */
/* solve triangular system r(transpose)*x = b, r defined by                     */
/* Householder-QR-decomposition decomp (with column scaling]                    */
/* **************************************************************************** */
void o8solt(INTEGER nlow,INTEGER nup,DOUBLE b[],DOUBLE x[]) {

    static INTEGER  i,j;
    static DOUBLE   sum;

    for (i = nlow ; i <= nup ; i++) {
    
        /* b has been permuted already ! */
        
        x[i] = b[i]*cscal[colno[i]];
    }
    for (i = nlow ; i <= nup ; i++) {
        sum = zero;
        for (j = nlow ; j <= i-1 ; j++) {
            sum = sum+qr[j][i]*x[j];
        }
        x[i] = (x[i]-sum)/diag[i];
    }
    return;
}
    
/* **************************************************************************** */
/* lenght of vector (a,b). numerically stable version with                      */
/* overflow / underflow saveguard                                               */
/* **************************************************************************** */
DOUBLE o8dsq1(DOUBLE a,DOUBLE b) {

    static DOUBLE   a1,b1;
    static DOUBLE   res;

    a1 = myfabs(a);
    b1 = myfabs(b);
    if ( a1 > b1 ) {
        res = a1*sqrt(one+pow(b1/a1,2));
    } else {
        if ( b1 > a1 ) {
            res = b1*sqrt(one+pow(a1/b1,2));
        } else {
            res = a1*sqrt(two);
        }
    }
    return (res);
}

/* **************************************************************************** */
/*                  computes the upper triangular Cholesky-factor               */
/* **************************************************************************** */
void o8upd(DOUBLE r[][NX+1],DOUBLE z[],DOUBLE y[],INTEGER n,LOGICAL *fail) {
    
    DOUBLE  o8dsq1(DOUBLE a,DOUBLE b);
    void    o8left(DOUBLE a[][NX+1],DOUBLE b[],DOUBLE y[],DOUBLE *yl,INTEGER n);

    static INTEGER  i,j,i1;
    static DOUBLE   sdiag[NX+1],rn1[NX+1],o8w[NX+1];
    static DOUBLE   yl,zl,wl,wn1,ai,bi,h;

    /* o8upd computes the upper triangular Cholesky-factor               */
    /* r1 of r(transpose)*r+z*z(transpose)-y*y(transpose)                */
    /* and restores it in r. The strict lower triangle of r              */
    /* remains unchanged.                                                */
    /* fail = TRUE if the decomposition does'nt exist, stop on dimension */
    /* error, fail = false on success.                                   */
    
    if ( n > NX ) {
        fprintf(stderr,"o8upd wrong call\n");
        exit(1);
    }
    *fail = FALSE;
    
    /* save subdiagonal */
    
    for (i = 1 ; i <= n-1 ; i++) {
        sdiag[i]  = r[i+1][i];
        r[i+1][i] = zero;
    }
    /* step one: include z*z(transpose) */
    
    zl = zero;
    for (i = 1 ; i <= n ; i++) {
        zl = zl + pow(z[i],2);
    }
    if ( zl != zero ) {
    
        /* solve w(transpose)*r = z(transpose) */
        
        o8left(r,z,o8w,&wl,n);
        
        wl = sqrt(wl+one);
        
        /* u[2]*u[3]*...*u[n]*w = ( norm(w),0,..,0)(transpose) */
        /* u[i] rotations                                      */
        
        for (i = n ; i >= 2 ; i--) {
            if ( o8w[i] != zero ) {
                i1        = i-1;
                ai        = o8w[i1];
                bi        = o8w[i];
                o8w[i1]     = o8dsq1(ai,bi);
                ai        = ai/o8w[i1];
                bi        = -bi/o8w[i1];
                r[i][i1]  = bi*r[i1][i1];
                r[i1][i1] = ai*r[i1][i1];
                for (j = i ; j <= n ; j++) {
                    h        = ai*r[i1][j] - bi*r[i][j];
                    r[i][j]  = bi*r[i1][j] + ai*r[i][j];
                    r[i1][j] = h;
                }
            }
        }
        /* r = d*r, d = diag(wl,1,...,1), r now Hessenberg */
        
        for (i = 1 ; i <= n ; i++) {
            r[1][i] = r[1][i]*wl;
        }
        /* r = u[n-1]*...*u[1]*r now upper triangular, */
        /* u[i] givens-rotations                       */
        
        for (i = 1 ; i <= n-1 ; i++) {
            i1 = i+1;
            ai = r[i][i];
            bi = -r[i1][i];
            h  = o8dsq1(ai,bi);
            if ( h != zero ) {
                ai       = ai/h;
                bi       = bi/h;
                r[i][i]  = h;
                r[i1][i] = zero;
                for (j = i+1 ; j <= n ; j++) {
                    h        = ai*r[i][j] - bi*r[i1][j];
                    r[i1][j] = bi*r[i][j] + ai*r[i1][j];
                    r[i][j]  = h;
                }
            }
        }
    }
    /* step two : include -y*y(transpose) */
    
    yl = zero;
    for (i = 1 ; i <= n ; i++) {
        yl = yl + pow(y[i],2);
    }
    if ( yl != zero ) {
    
        o8left(r,y,o8w,&wl,n);
        
        if ( wl >= one ) {
            *fail = TRUE;
        } else {
            wl  = sqrt(myfabs(one-wl));
            wn1 = wl;
            
            /* ************************************************ */
            /*   ( r(new), 0 )                  (    r  , w )   */
            /*   (-----------)  =  u[1]*...u[n]*(-----------)   */
            /*   (y(transp),1)                  ((0,..,0),wl)   */
            /* ************************************************ */
            
            for (i = n ; i >= 1 ; i--) {
                ai  = wn1;
                bi  = o8w[i];
                wn1 = o8dsq1(ai,bi);
                if ( wn1 != zero ) {
                    ai      = ai/wn1;
                    bi      = bi/wn1;
                    rn1[i]  = bi*r[i][i];
                    r[i][i] = ai*r[i][i];
                    for (j = i+1 ; j <= n ; j++) {
                        h       = ai*r[i][j] - bi*rn1[j];
                        rn1[j]  = bi*r[i][j] + ai*rn1[j];
                        r[i][j] = h;
                    }
                }
            }
        }
    }
    /* restore subdiagonal */
    
    for ( i = 1 ;  i <= n-1 ;  i++) {
        r[i+1][i] = sdiag[i];
    }
    return;
}

/* **************************************************************************** */
/*                                    o8rght                                    */
/* **************************************************************************** */
void o8rght(DOUBLE a[][NX+1],DOUBLE b[],DOUBLE y[],DOUBLE *yl,INTEGER n) {
    
    static INTEGER  i,j;
    static DOUBLE   h;

    /* o8rght assumes that the Cholesky-factor of a         */
    /* a = r(transpose)*r is stored in the upper half of a. */
    /* b is a right hand side. o8rght solves                */
    /* r*y = b                                              */
    /* yl = pow(norm(y),2)                                  */
    
    *yl = zero;
    for (i = n ; i >= 1 ; i--) {
        h = b[i];
        for (j = i+1 ; j <= n ; j++) {
            h = h - a[i][j]*y[j];
        }
        h    = h/a[i][i];
        y[i] = h;
        *yl  = pow(h,2) + *yl;
    }
    return;
}

/* **************************************************************************** */
/*                                  o8left                                      */
/* **************************************************************************** */
void o8left(DOUBLE a[][NX+1],DOUBLE b[],DOUBLE y[],DOUBLE *yl,INTEGER n) {
    
    static INTEGER  i,j;
    static DOUBLE   h;

    /* o8left assumes that the Cholesky-factor of a         */
    /* a = r(transpose)*r is stored in the upper half of a. */
    /* b is a right hand side. o8left solves                */
    /* r(transpose)*y = b                                   */
    /* yl = pow(norm(y),2)                                  */
    
    *yl = zero;
    for (i = 1 ; i <= n ; i++) {
        h = b[i];
        for (j = 1 ; j <= i-1 ; j++) {
            h = h - a[j][i]*y[j];
        }
        h    = h/a[i][i];
        y[i] = h;
        *yl  = pow(h,2) + *yl;
    }
    return;
}

/* **************************************************************************** */
/*                      euclidean norm of x , avoid overflow                    */
/* **************************************************************************** */
DOUBLE o8vecn(INTEGER nl,INTEGER nm,DOUBLE x[]) {

    static INTEGER  i;
    static DOUBLE   xm,h;

    if ( nm < nl ) {
    
        return (zero);
    }
    xm = myfabs(x[nl]);
    for (i = nl+1 ; i <= nm ; i++) {
        xm = patMax(xm,myfabs(x[i]));
    }
    if ( xm == zero ) {
    
        return (zero);
        
    } else {
        h = zero;
        for (i = nl ; i <= nm ; i++) {
            h = h+pow(x[i]/xm,2);
        }
        return (xm*sqrt(h));
    }
}
    
/* **************************************************************************** */
/*  solution of extended quadratic program                                      */
/*  scf*gradf(x)*d+(1/2)*d*a*d+summe( tauqp*d[i] +( my/2)*pow(d[i],2) )         */
/*  minimal subject to                                                          */
/*  d[i] >= 0, i = 1,...,nr                                                     */
/*  (gres[.][j]*d+res[j])+vz*d[j] = 0, j = 1,...,nh vz = -sign(res[j])          */
/*  (gres[.][DONLP2alist[j]])*d+res[DONLP2alist[j]])+d[j] >= 0, j = nh+1,....,nr            */
/*  the weight tauqp is adapted during solution                                 */
/*  the quasi-Newton-matrix a is taken from o8comm.h                            */
/*  a is regularized if not sufficiently well conditioned                       */
/*  the resulting d[1+nr],...,d[n+nr] is a direction of descent for             */
/*  the Zangwill function of the corresponding nonlinear                        */
/*  optimization problem                                                        */
/*  f(x) = min, res[j] = 0, j = 1,..nh, res[j] >= 0 , j = nh+1,nres             */
/*  at the currrent point x if the weight tauqp is chosen appropriately         */
/*  the quadratic programming problem is solved using the method                */
/*  of Goldfarb and Idnani                                                      */
/*  variables are stored in xd (solution) and ud (multipliers)                  */
/*  in the following order xd = ( d[j], j = 1,nr; d = direct. of desc.)         */
/*  ud = ( multipliers for d[i] >= 0 , i = 1,..,nr ;                            */
/*  multipliers for the equality constraints ,                                  */
/*  multipliers for the general inequality constraints )                        */
/* **************************************************************************** */
void o8qpdu(void) {
    
    void    o8info(INTEGER icase);
    void    o8msg (INTEGER num);
    void    o8dird(void);
    void    o8cutd(void);
    DOUBLE  o8sc1 (INTEGER i,INTEGER j,DOUBLE a[],DOUBLE b[]);
    DOUBLE  o8sc3 (INTEGER n,INTEGER m,INTEGER j,DOUBLE a[][NRESM+1],DOUBLE b[]);
    DOUBLE  o8vecn(INTEGER nl,INTEGER nm,DOUBLE x[]);
    void    o8zup (DOUBLE z[]);
    void    o8rup (DOUBLE rv[]);
    void    o8dlcd(INTEGER ai[],INTEGER l);
    void    o8adcd(void);
    void    o8rinv(INTEGER n,DOUBLE a[][NX+1],INTEGER ndual,DOUBLE x[][NDUALM+1]);
    
    static DOUBLE   g0[NDUALM+1],ci0[MDUALM+1],cii[MDUALM+1],
                    cei[NDUALM+1],
                    xd[NDUALM+1],s[MDUALM+1],z[NDUALM+1],vr[MDUALM+1];
    static DOUBLE   xdold[NDUALM+1],udold[MDUALM+1];
    static INTEGER  ai[MDUALM+1],iai[MDUALM+1],iaexcl[MDUALM+1],aiold[MDUALM+1];
    static DOUBLE   y[NDUALM+1],mult[NRESM+1];
    
    static DOUBLE   infe1,s1,s2,tiny,
                    my,zz,ss,su,t,t1,t2,f,psid,c1,c2,cdiag,term,
                    su1,su2,condr,infiny,term1,term2,
                    diff0;
    static INTEGER  i,j,k,ip,l,qpdel[MDUALM+1],incr;
    static LOGICAL  wlow;

    ndual = n+nr;
    
    /* number of equality constraints in QP-problem */
    
    me = nh;
    
    /* QP inequality constraints = active constraints - */
    /* equality-constraints + slack's                   */
    
    mi     = 2*nr-nh;
    infiny = epsmac/tolmac;
    if ( analyt ) {
        tiny = 2*nr*epsmac*tp3;
    } else {
        tiny = 2*nr*patMax(epsdif,epsmac*tp3);
    }
    qpterm = 0;
    for (i = 1 ; i <= nr ; i++) {
    
        /* check gradients of active constraints against zero */
        
        for (j = 1 ; j <= n ; j++) {
            y[j] = gres[j][DONLP2alist[i]];
        }
        mult[i] = one;
        if ( o8vecn(1,n,y) == zero ) {
            mult[i] = zero;
            if ( ! silent ) o8msg(8);
        }
    }
    /* restart point in case of increase of tauqp */
    
    L10:

    /* initialize matrices j and r */
    
    for (i = 1 ; i <= ndual ; i++) {
        ddual[i] = zero;
        for (j = 1 ; j <= ndual ; j++) {
            r[j][i]  = zero;
            xj[j][i] = zero;
        }
    }
    rnorm = one;
    rlow  = one;
    term1 = zero;
    for (i = 1 ; i <= nres ; i++) {
        u[i] = zero;
        if ( o8w[i] > term1 ) term1 = o8w[i];
    }
    accinf[itstep][19] = clow;
    accinf[itstep][20] = term1;
    accinf[itstep][31] = tauqp;
    for (i = 1 ; i <= me+mi ; i++) {
        ud[i] = zero;
    }
    c1 = myfabs(a[1][1]);
    for (i = 1 ; i <= n ; i++) {
        c1 = patMax(c1,myfabs(a[i][i]));
    }
    c1 = c1*tp1;
    
    /* we require much more regularity of a in the singular case */
    
    for (i = 1 ; i <= n ; i++) {
        if ( myfabs(a[i][i]) < sqrt(rho1)*c1 ) a[i][i] = sqrt(rho1)*c1;
    }
    /* invert the Cholesky-factor and store in xj ( Idnanis j-matrix) */
    
    o8rinv(n,a,ndual,xj);
    
    c1   = myfabs(a[1][1]);
    incr = nr;
    c2   = myfabs(xj[1+incr][1+incr]);
    for (i = 1 ; i <= n ; i++) {
        c1 = patMax(c1,myfabs(a[i][i]));
        c2 = patMax(c2,myfabs(xj[i+incr][i+incr]));
    }
    my = zero;
    for (i = 1 ; i <= n ; i++) {
        for (j = i ; j <= n ; j++) {
             my = my+pow(a[i][j],2);
        }
    }
    my    = my/n;
    cdiag = one/sqrt(my);
    for ( i = 1 ; i <= incr ; i++) {
        xj[i][i] = cdiag;
    }
    for (i = 1 ; i <= ndual ; i++) {
        if ( i >  incr ) g0[i] = gradf[i-incr]*scf;
        if ( i <= incr ) g0[i] = tauqp;
    }
    /* compute unconstrained solution */
    /* the Cholesky-factor of a is stored in the upper triangle */
    
    for (i = 1 ; i <= n ; i++) {
        su = zero;
        for (j = 1 ; j <= i-1 ; j++) {
            su = su+a[j][i]*y[j+incr];
        }
        y[i+incr] = (g0[i+incr]-su)/a[i][i];
    }
    for (i = n ; i >= 1 ; i--) {
        su = zero;
        for (j = i+1 ; j <= n ; j++) {
            su = su+a[i][j]*y[j+incr];
        }
        y[i+incr] = (y[i+incr]-su)/a[i][i];
    }
    for (i = 1 ; i <= incr ; i++) {
    
        /* initially assume the slacks being zero */
        
        y[i] = zero;
    }
    for (i = 1 ; i <= ndual ; i++) {
        xd[i] = -y[i];
    }
    /* unconstrained minimizer of the QP: slacks come first */
    
    f = p5*o8sc1(1,ndual,g0,xd);
    
    /* define the initial working set: all slacks are at their */
    /* lower bounds                                            */
    
    iq = nr;
    for (i = 1 ; i <= iq ; i++) {
        ai[i]   = i;
        r[i][i] = one;
        ud[i]   = tauqp;
    }
    rnorm = one;
    rlow  = one;
    
    /* introduction of equality constraints */
    
    for (i = 1 ; i <= me ; i++) {
        for ( j = 1 ; j <= iq ; j ++) {
            ud1[j] = ud[j];
        }
        L20:

        ud1[iq+1] = zero;
        
        /* an equality constraint is indicated by the negative index */
    
        ai[iq+1] = -nr-i;
        for (j = 1 ; j <= n ; j++) {
        
            /* me = nh and DONLP2alist[i] = i for i = 1,..,nh */
            
            cei[j+incr] = gres[j][i];
        }
        for (j = 1 ; j <= incr ; j++) {
            cei[j] = zero;
        }
        cei[i] = one;
        if ( res[i] > zero ) cei[i] = -one;
        for (j = 1 ; j <= ndual ; j++) {
            np[j] = cei[j];
        }
        o8zup(z);
        
        if ( iq != 0 ) o8rup(vr);
        
        /* |z| = 0? */
        
        /* old      zz = o8sc1(1,ndual,z,z)     */
            
        zz   = o8vecn(1,ndual,z);
        term = o8sc1(1,ndual,z,np);
        
        /* old      if (zz != zero && term > zero ) { */
        
        if ( zz >= tiny*rnorm && term > zero ) {
            t2 = (-o8sc1(1,ndual,np,xd)-res[i])/term;
        } else if ( (-o8sc1(1,ndual,np,xd)-res[i]) >= zero ) {
            t2 = infiny;
        } else {
            t2 = -infiny;
        }
        /* in addition of an equality constraint, t2 may be positive or negative*/
        
        if ( iq != 0 ) o8rup(vr);
        l = 0;
        if ( t2 > zero ) {
            t1 = infiny;
            for (k = 1 ; k <= iq ; k++) {
                if( vr[k] > zero && ai[k] > 0 ) {
                    if( ud1[k]/vr[k] < t1 ) {
                        t1 = ud1[k]/vr[k];
                    }
                }
            }
            t = patMin(t1,t2);
        } else {
            t1 = infiny;
            for ( k = 1 ; k <= iq ; k++) {
                if( vr[k] < zero && ai[k] > 0 ) {
                    if( ud1[k]/myfabs(vr[k]) < t1 ) {
                        t1 = ud1[k]/myfabs(vr[k]);
                    }
                }
            }
            t1 = -t1;
            t  = patMax(t1,t2);
            
            /* t now negative */
        }
        /* add constraint , otherwise we must first delete some */
        /* inequality constraint with negative multiplier       */
        /* first delete then add!                               */
        
        if ( myfabs(t)  >= infiny ) goto L2000;
        
        if ( myfabs(t2) >= infiny ) {
        
            /* purely dual step */
            
            for (k = 1 ; k <= iq ; k++) {
                ud1[k] = ud1[k]+t*(-vr[k]);
                if ( ud1[k] < zero && ai[k] > 0 ) ud1[k] = zero;
            }
            ud1[iq+1] = ud1[iq+1]+t;
            qpdel[0]  = 0;
            for (j = 1 ; j <= iq ; j++) {
                if ( ud1[j] <= tiny && ai[j] > 0 ) {
                    qpdel[0]        = qpdel[0]+1;
                    qpdel[qpdel[0]] = ai[j];
                }
            }
            for (k = 1 ; k <= qpdel[0] ; k++) {
                l      = qpdel[k];
                iai[l] = l;
                
                o8dlcd(ai,l);
            }
            goto L20;
        }
        for (k = 1 ; k <= ndual ; k++) {
            xd[k] = xd[k]+t*z[k];
        }
        for (k = 1 ; k <= iq ; k++) {
            ud1[k] = ud1[k]+t*(-vr[k]);
            if ( ud1[k] < zero && ai[k] > 0 ) ud1[k] = zero;
        }
        ud1[iq+1] = ud1[iq+1]+t;
        
        f = f+t*o8sc1(1,ndual,z,np)*(p5*t+ud1[iq+1]);
        
        if ( myfabs(t2-t1) <= tiny ) {
            qpdel[0] = 0;
            for (j = 1 ; j <= iq ; j++) {
                if ( ud1[j] <= tiny && ai[j] > 0 ) {
                    qpdel[0]        = qpdel[0]+1;
                    qpdel[qpdel[0]] = ai[j];
                }
            }
            for (k = 1 ; k <= qpdel[0] ; k++) {
                l      = qpdel[k];
                iai[l] = l;
                
                o8dlcd(ai,l);
            }
            ai[iq+1] = -i-nr;
            
            o8adcd();
            
        } else if ( t == t2 ) {
            ai[iq+1] = -i-nr;
            
            o8adcd();
            
        } else {
            qpdel[0] = 0;
            for (j = 1 ; j <= iq ; j++) {
                if ( ud1[j] <= tiny && ai[j] > 0 ) {
                    qpdel[0]        = qpdel[0]+1;
                    qpdel[qpdel[0]] = ai[j];
                }
            }
            for (k = 1 ; k <= qpdel[0] ; k++) {
                l      = qpdel[k];
                iai[l] = l;
                
                o8dlcd(ai,l);
            }
            goto L20;
        }
        for (j = 1 ; j <= iq ; j++) {
            ud[j] = ud1[j];
        }  
    }
    /* set iai = k\ai */
    
    for (i = 1 ; i <= mi ; i++) {
        iai[i] = i;
    }
    /* step 1 */
    
    L50:

    /* ai = QP - working set , iai[i] = 0 if i in ai */
    
    for (i = 1 ; i <= iq ; i++) {
        ip = ai[i];
        if ( ip > 0 ) iai[ip] = 0;
    }
    /* s[xd] = ci(trans)*xd+ci0 >= 0 ? */
    
    psid = zero;
    
    /* psid : the measure of infeasibility */
    
    for (i = 1 ; i <= mi ; i++) {
    
        /* iaexcl: if = 0, exclude from addition in this cycle */
        
        iaexcl[i] = 1;
        su        = zero;
        
        /* numbers of inequality constraints:                           */
        /* i = 1,...,nr corresponds to the constraints v >= 0, u_a >= 0 */
        /* i = nr+1,....,mi to the regularized general inequalities     */
        
        if ( i > nr ) {
            k = DONLP2alist[i+nh-incr];
            for (j = 1 ; j <= n ; j++) {
                cii[j+incr] = gres[j][k];
            }
            for (j = 1 ; j <= incr ; j++) {
                cii[j] = zero;
            }
            cii[nh+i-incr] = one;
            ci0[i]         = res[k];
        } else {
            for (j = 1 ; j <= ndual ; j++) {
                cii[j] = zero;
            }
            ci0[i] = zero;
            cii[i] = one;
        }
        su   = o8sc1(1,ndual,cii,xd)+ci0[i];
        s[i] = su;
        psid = psid+patMin(zero,su);
    }
    for (i = 1 ; i <= iq ; i++) {
        udold[i] = ud[i];
        aiold[i] = ai[i];
    }
    for (i = 1 ; i <= ndual ; i++) {
        xdold[i] = xd[i];
    }
    L60:

    ss = zero;
    ip = 0;
    
    /* introduce most violated inequality constraint */
    
    for (i = 1 ; i <= mi ; i++) {
        if( s[i] < ss && iai[i] != 0 && iaexcl[i] != 0 ) {
            ss = s[i];
            ip = i;
        }
    }
    if ( iq > 1 ) {
        condr = rnorm/rlow;
    } else {
        condr = one;
    }
    if ( myfabs(psid) <= tiny*(c1*c2+condr) || ip == 0) {
    
        /* successful termination of QP-solver for current tauqp */
        
        if ( myfabs(psid) > tiny*(c1*c2+condr) && ! silent ) o8msg(10);
        qpterm = 1;
        accinf[itstep][30] = one;
        accinf[itstep][13] = condr;
        accinf[itstep][14] = c1*c2;
        for (i = 1 ; i <= n ; i++) {
            d[i] = xd[i+incr];
        }
        /* new : dnorm added */
        
        dnorm  = o8vecn(1,n,d);
        infeas = zero;
        
        for (i = 1 ; i <= incr ; i++) {
            infeas = infeas+myfabs(xd[i]);
        }
        /* L1-norm of slack variables */
        
        accinf[itstep][31] = tauqp;
        accinf[itstep][32] = infeas;
        wlow = FALSE;
        su1  = zero;
        su2  = zero;
        for (i = 1 ; i <= iq ; i++) {
            if ( ai[i] < 0 ) {
                u[-(ai[i]+nr)] = ud[i];
            } else {
                if ( ai[i] > nr ) u[DONLP2alist[ai[i]+nh-nr]] = ud[i];
            }
        }
        term1 = zero;
        for (j = 1 ; j <= n ; j++) {
            np[j] = gradf[j]*scf;
        }
        for (i = 1 ; i <= nres ; i++) {
            for (j = 1 ; j <= n ; j++) {
                np[j] = np[j]-gres[j][i]*u[i];
            }
        }
        b2n = o8vecn(1,n,np);
        
        if ( scf != zero ) b2n = b2n/scf;
        
        /* correction in the original variables */
            
        infe1 = zero;
        for (i = 1 ; i <= nr ; i++) {
            infe1 = infe1+myfabs(xd[i])*mult[i];
        }
/*        if ( upsi <= delmin*nres && upsi0 <= delmin*nres   */ /*old*/
        if ( upsi <= delmin*nres 
            && b2n <= (gfn+one)*epsx*tp2 && phase >= 0
            && infeas <= delmin*nres ) {
         
            /* since multipliers may be incorrect for infeas != zero be careful */
            /* we consider the problem as successfully solved with reduced      */
            /* requirements                                                     */

            for (i = 1 ; i <= n ; i++) {
                d[i] = zero;
            }
            dirder = zero;
            dnorm  = zero;
            optite = three;
            qpterm = 1 ;   /* new */
            return ;
        }
        /* there may be an additional increase of tauqp necessary again */

        if ( infe1 >= (one-delta1/tauqp)*upsi && 
            ( o8vecn(1,n,d) <= patMin(infe1,pow(infe1,2))*tp1
            || upsi > tau0*p5 ) ) {
            
            /* further increase tauqp ! */
            
            for (i = 1 ; i <= nres ; i++) {
                u[i]     = zero;
                slack[i] = zero;
            }
            if ( ! silent ) o8msg(17);
            if ( tauqp*taufac > taumax ) {
                if ( ! silent ) o8msg(5);
                qpterm = -1;
                accinf[itstep][30] = qpterm;
                accinf[itstep][31] = tauqp;
                accinf[itstep][32] = infeas;
                for (i = 1 ; i <= n ; i++) {
                    d[i] = zero;
                }
                dnorm = zero;
                dirder = zero ;
                return;
                
            } else {
                tauqp = tauqp*taufac;
                
                goto L10;
            }
        }
        /* compute new weights for the penalty-function */
        
	//        L500:

        for (i = 1 ; i <= nres ; i++) {
            slack[i] = zero;
        }
        for (i = 1 ; i <= nr ; i++) {
            slack[DONLP2alist[i]] = xd[i];
        }
        wlow = FALSE;
        for (i = 1 ; i <= nres ; i++) {
            w1[i] = o8w[i];
            if ( i <= nh ) {
                if ( myfabs(slack[i]) > myfabs(res[i])+tiny ) {
                    w1[i] = myfabs(u[i]);
                } else {
                    w1[i] = ny*myfabs(u[i])+tau;
                }
            } else {
                if ( bind[i] == 0 ) {
                    w1[i] = patMax(o8w[i]*p8,tau);
                } else {
                    if ( res[i] >= zero && slack[i] <= tiny )
                    w1[i] = patMax(ny*myfabs(u[i])+tau,(myfabs(u[i])+w1[i])*p5);
                    if ( res[i] >= zero && slack[i] > tiny )
                    w1[i] = myfabs(u[i]);
                    if ( res[i] < zero && slack[i] <= -res[i]+tiny )
                    w1[i] = patMax(ny*myfabs(u[i])+tau,(myfabs(u[i])+w1[i])*p5);
                    if ( res[i] < zero && slack[i] >  -res[i]+tiny )
                    w1[i] = myfabs(u[i]);
                }
            } 
            if ( w1[i] < o8w[i] ) wlow = TRUE;
        }
        if ( wlow ) {
            s1 = zero;
            s2 = zero;
            for (i = 1 ; i <= nres ; i++) {
                if ( i <= nh ) {
                    s1 = s1+w1[i]*myfabs(resst[i]);
                    s2 = s2+w1[i]*myfabs(res[i]);
                } else {
                    s1 = s1-patMin(zero,resst[i])*w1[i];
                    s2 = s2-patMin(zero,res[i]  )*w1[i];
                }
            }
            diff0 = (fxst-fx)*scf+(s1-s2);
            if ( diff0 >= eta*clow && itstep-lastdw >= patMax(5,patMin(n/10,20)) ) {
            
                /* accept new (diminished ) weights */
            
                lastdw = itstep;
                lastch = itstep;
                level  = diff0/iterma;
                psist  = s1;
                psi    = s2;
                for (i = 1 ; i <= nres ; i++) {
                    if ( w1[i] != o8w[i] ) lastch = itstep;
                    o8w[i] = w1[i];
                }
                clow = clow+one;
                if ( clow > itstep/10 ) {
                
                    /* additional increase of eta */
                    
                    eta = eta*onep3;
                    if ( ! silent ) o8info(11);
                }
                if ( ! silent ) o8info(12);
                
                goto L1000;
            }
        }
        /* we cannot accept new weights */
        /* reset weights                */
        
        for (i = 1 ; i <= nres ; i++) {
            w1[i] = o8w[i];
            if ( i <= nh ) {
                if ( slack[i] > myfabs(res[i]) ) w1[i] = myfabs(u[i]);
                if ( slack[i] <= myfabs(res[i]) ) {
                    if ( o8w[i] <= myfabs(u[i]) && myfabs(u[i]) <= o8w[i]+tau ) {
                        w1[i] = o8w[i]+two*tau;
                    } else {
                        w1[i] = patMax(o8w[i],ny*myfabs(u[i])+tau);
                    }
                }
            } else {
                if ( slack[i] > -patMin(-tiny,res[i]) && bind[i] == 1 ) {
                    w1[i] = myfabs(u[i]);
                } else if ( bind[i] == 1 && slack[i] <= -patMin(-tiny,res[i])
                            && u[i] <= o8w[i]+tau && o8w[i] >= u[i] ) {
                    
                    w1[i] = o8w[i]+two*tau;
                } else if ( bind[i] == 1 ) {
                    w1[i] = patMax(o8w[i],ny*myfabs(u[i])+tau);
                }
            }
        }
        term1 = zero;
        for (i = 1 ; i <= nres ; i++) {
            if ( w1[i] > o8w[i] || w1[i] < o8w[i] ) lastch = itstep;
            if ( w1[i] > o8w[i] ) lastup = itstep;
            if ( w1[i] < o8w[i] ) lastdw = itstep;
            o8w[i]  = w1[i];
            term1 = patMax(term1,o8w[i]);
        }
        s1 = zero;
        s2 = zero;
        for (i = 1 ; i <= nres ; i++) {
            if ( i <= nh ) {
                s1 = s1+o8w[i]*myfabs(resst[i]);
                s2 = s2+o8w[i]*myfabs(res[i]);
            } else {
                s1 = s1-o8w[i]*patMin(zero,resst[i]);
                s2 = s2-o8w[i]*patMin(zero,res[i]);
            }
        }
        psist = s1;
        psi   = s2;
        if ( ! silent ) o8info(12);
        accinf[itstep][20] = term1;
        accinf[itstep][19] = clow;
        
        goto L1000;
    }
    if ( ip > nr ) {
        k = DONLP2alist[ip+nh-nr];
        for (j = 1 ; j <= n ; j++) {
            cii[j+incr] = gres[j][k];
        }
        for (j = 1 ; j <= incr ; j++) {
            cii[j] = zero;
        }
        cii[nh+ip-nr] = one;
        ci0[ip]       = res[k];
    } else {
        for (j = 1 ; j <= ndual ; j++) {
            cii[j] = zero;
        }
        ci0[ip] = zero;
        cii[ip] = one;
    }
    for (i = 1 ; i <= ndual ; i++) {
        np[i] = cii[i];
    }
    for (i = 1 ; i <= iq ; i++) {
        ud1[i] = ud[i];
    }
    ud1[iq+1] = zero;
    ai[iq+1]  = ip;
    
    L100:

    /* step 2a */
    
    o8zup(z);
    
    if ( iq != 0 ) o8rup(vr);
    
    l  = 0;
    t1 = infiny;
    for (k = 1 ; k <= iq ; k++) {
        if(ai[k] > 0 && vr[k] > zero) {
            if ( ud1[k]/vr[k] < t1 ) {
                t1 = ud1[k]/vr[k];
            }
        }
    }
    /* |z| = 0? */
    
    /* old      zz = o8sc1(1,ndual,z,z) */
    
    zz   = o8vecn(1,ndual,z);
    term = o8sc1(1,ndual,z,np);
    
    /* old      if (zz != zero && term > zero ) { */
    
    if ( zz >= tiny*rnorm && term > zero ) {
        t2 = -s[ip]/term;
    } else {
        t2 = infiny;
    }
    t = patMin(t1,t2);
    if ( t >= infiny ) goto L2000;

    if( t2 >= infiny ) {
        for (k = 1 ; k <= iq ; k++) {
            ud1[k] = ud1[k]+t*(-vr[k]);
            if ( ud1[k] < zero && ai[k] > 0 ) ud1[k] = zero;
        }
        ud1[iq+1] = ud1[iq+1]+t;
        qpdel[0]  = 0;
        for (i = 1 ; i <= iq ; i++) {
            if ( ud1[i] <= tiny && ai[i] > 0 ) {
                qpdel[0]        = qpdel[0]+1;
                qpdel[qpdel[0]] = ai[i];
            }
        }
        for (k = 1 ; k <= qpdel[0] ; k++) {
            l      = qpdel[k];
            iai[l] = l;
            
            o8dlcd(ai,l);
        }
        goto L100;
    }
    for (k = 1 ; k <= ndual ; k++) {
        xd[k] = xd[k]+t*z[k];
    }
    for (k = 1 ; k <= iq ; k++) {
        ud1[k] = ud1[k]+t*(-vr[k]);
        if ( ud1[k] < zero && ai[k] > 0 ) ud1[k] = zero;
    }
    ud1[iq+1] = ud1[iq+1]+t;
    
    f = f+t*o8sc1(1,ndual,z,np)*(p5*t+ud1[iq+1]);
    
    if ( t2 <= t1-tiny ) {
    
        /* ddual is computed by o8zup */
        
        if ( o8vecn(iq+1,ndual,ddual) < epsmac*rnorm ) {
        
            /* degeneracy: adding this constraint gives a singular working set */
            /* theoretically impossible, but due to roundoff this may occur.   */
            /* mark this constraint and try to add another one                 */
            
            iptr = ip;
            iqtr = iq;
            for (i = 1 ; i <= iq ; i++) {
                aitr[i] = ai[i];
            }
            sstr  = ss;
            riitr = o8vecn(iq+1,ndual,ddual);
            
            if ( ! silent ) o8msg(19);
            iaexcl[ip] = 0;
            for (i = 1 ; i <= mi ; i++) {
                iai[i] = i;
            }
            for (i = 1 ; i <= iq ; i++) {
                ai[i] = aiold[i];
                if (ai[i] > 0 ) iai[ai[i]] = 0;
                ud1[i] = udold[i];
            }
            for (i = 1 ; i <= ndual ; i++) {
                xd[i] = xdold[i];
            }
            goto L60;
        }
        /* add constraint, l-pair */
        
        o8adcd();
        
        iai[ip] = 0;
        for (i = 1 ; i <= iq ; i++) {
            ud[i] = ud1[i];
        }
        goto L50;
    }
    su = zero;
    if ( ip > nr ) {
    
        /* a general linear inequality constraint */
        
        k = DONLP2alist[ip+nh-nr];
        for (j = 1 ; j <= n ; j++) {
            cii[j+incr] = gres[j][k];
        }
        for (j = 1 ; j <= incr ; j++) {
            cii[j] = zero;
        }
        cii[nh+ip-nr] = one;
        ci0[ip]       = res[k];
        
        s[ip] = o8sc1(1,ndual,cii,xd)+ci0[ip];
    } else {
    
        /* a slack constraint */
        
        s[ip] = xd[ip];
    }
    /* now t = t1 */
    
    qpdel[0] = 0;
    for (i = 1 ; i <= iq ; i++) {
        if ( ud1[i] <= tiny && ai[i] > 0 ) {
            qpdel[0]        = qpdel[0]+1;
            qpdel[qpdel[0]] = ai[i];
        }
    }
    for (k = 1 ; k <= qpdel[0] ; k++) {
        l      = qpdel[k];
        iai[l] = l;
        
        o8dlcd(ai,l);
    }
    if ( t2 <= t1+tiny ) {
        if ( o8vecn(iq+1,ndual,ddual) < epsmac*rnorm ) {
        
            /* degeneracy */
            
            iptr = ip;
            iqtr = iq;
            for (i = 1 ; i <= iq ; i++) {
                aitr[i] = ai[i];
            }
            sstr  = ss;
            riitr = o8vecn(iq+1,ndual,ddual);
            
            if ( ! silent ) o8msg(19);
            iaexcl[ip] = 0;
            for (i = 1 ; i <= mi ; i++) {
                iai[i] = i;
            }
            for (i = 1 ; i <= iq ; i++) {
                ai[i] = aiold[i];
                if ( ai[i] > 0 ) iai[ai[i]] = 0;
                ud1[i] = udold[i];
            }
            for (i = 1 ; i <= ndual ; i++) {
                xd[i] = xdold[i];
            }
            goto L60;
        }
        /* add constraint, l-pair */
        
        o8adcd();
        
        iai[ip] = 0;
        for (i = 1 ; i <= iq ; i++) {
            ud[i] = ud1[i];
        }
        goto L50;
        
    } else {
    
        goto L100;
    }
    /* this is the exit point of o8qpdu                                     */
    /* we either may have successful or unsuccessful termination here       */
    /* the latter with qpterm = -2 or -3, in which case it may nevertheless */
    /* be possible to use the computed d. -2 or -3 exit is theoretically    */
    /* impossible but may occur due to roundoff effects.                    */
    /* we check the directional derivative of the penalty-function now      */
    
    L1000:

    /* cut and rescale d if appropriate */
    
    o8cutd();
    
    /* compute the directional derivative dirder */
    
    o8dird();
    
    if ( dirder >= zero || ( -dirder <= epsmac*tp2*(scf*myfabs(fx)+psi+one) &&
        infeas > patMax(upsi,nres*delmin) ) ) {
        if ( ! silent ) o8msg(18);
        if ( tauqp <= taumax/taufac ) {
            tauqp = tauqp*taufac;
            
            goto L10;
            
        } else {
            if ( ! silent ) o8msg(5);
            qpterm = -1;
            accinf[itstep][30] = qpterm;
            accinf[itstep][31] = tauqp;
            accinf[itstep][32] = infeas;
            for (i = 1 ; i <= n ; i++) {
                d[i] = zero;
            }
            dnorm = zero;
        }
    }
    return;
    
    L2000:

    /* QP infeasible ( in this application impossible , theoretically) */
    
    if ( ! silent ) o8msg(20);
    qpterm = -2;
    accinf[itstep][30] = -two;
    accinf[itstep][13] = condr;
    accinf[itstep][14] = c1*c2;
    for (i = 1 ; i <= n ; i++) {
        d[i] = xd[i+incr];
    }
    dnorm = o8vecn(1,n,d);
    su1   = zero;
    for (i = 1 ; i <= incr ; i++) {
        su1 = su1+myfabs(xd[i]);
    }
    /* L1-norm of slack variables */
    infeas = su1 ;
    accinf[itstep][32] = su1;
    wlow = FALSE;
    su1  = zero;
    su2  = zero;
    for (i = 1 ; i <= iq ; i++) {
        if ( ai[i] < 0 ) {
            u[-(ai[i]+nr)] = ud[i];
        } else {
            if ( ai[i] > nr ) u[DONLP2alist[ai[i]-nr+nh]] = ud[i];
        }
    }
    term1 = zero;
    for (j = 1 ; j <= n ; j++) {
        np[j] = gradf[j]*scf;
    }
    for (i = 1 ; i <= nres ; i++) {
        for (j = 1 ; j <= n ; j++) {
            np[j] = np[j]-gres[j][i]*u[i];
        }
        term  = o8sc3(1,n,i,gres,d);
        term2 = res[i]/patMax(one,gresn[i]);
        w1[i] = o8w[i];
        if ( i <= nh ) {
            if ( slack[i] > myfabs(res[i]) ) w1[i] = myfabs(u[i]);
            if ( slack[i] <= myfabs(res[i]) ) {
                if ( o8w[i] <= myfabs(u[i]) && myfabs(u[i]) <= o8w[i]+tau ) {
                    w1[i] = o8w[i]+two*tau;
                } else {
                    w1[i] = patMax(o8w[i],ny*myfabs(u[i])+tau);
                }
            }
            su1 = su1+myfabs(res[i]  )*w1[i];
            su2 = su2+myfabs(resst[i])*w1[i];
        } else {
            if ( slack[i] > -patMin(-tiny,res[i]) && bind[i] == 1 ) {
                w1[i] = myfabs(u[i]);
            } else if ( bind[i] == 1 && slack[i] <= -patMin(-tiny,res[i])
                        && u[i] <= o8w[i]+tau && o8w[i] >= u[i] ) {
                w1[i] = o8w[i]+two*tau;
            } else if ( bind[i] == 1 ) {
                w1[i] = patMax(o8w[i],ny*myfabs(u[i])+tau);
            }
            su1 = su1-w1[i]*patMin(zero,res[i]);
            su2 = su2-w1[i]*patMin(zero,resst[i]);
        }
        if ( o8w[i] != w1[i] ) lastch = itstep;
        o8w[i]  = w1[i];
        term1 = patMax(term1,o8w[i]);
    }
    psist = su2;
    psi   = su1;
    
    b2n   = sqrt(o8sc1(1,n,np,np));
    
    if ( scf != zero ) 
    { 
      b2n = b2n/scf;
    }
    else
    {
      b2n = -one;
    }
    if ( wlow ) {
        clow   = clow+one;
        lastch = itstep;
        lastdw = itstep;
    }
    if ( ! silent ) o8info(12);
    accinf[itstep][19] = clow;
    accinf[itstep][20] = term1;
    accinf[itstep][31] = tauqp;
    if ( upsi <= delmin*nres 
            && b2n <= (gfn+one)*epsx*tp2 && phase >= 0
            && infeas <= delmin*nres ) 
    {
         
      /* since multipliers may be incorrect for infeas != zero be careful */
      /* we consider the problem as successfully solved with reduced      */
      /* requirements                                                     */

            for (i = 1 ; i <= n ; i++) {
                d[i] = zero;
            }
            dnorm  = zero;
            optite = three;
            qpterm = 1 ;   /* new */
            return ;
    }

    
    goto L1000;
}

/* **************************************************************************** */
/*                  compute updated projected gradient (primal)                 */
/* **************************************************************************** */
void o8zup(DOUBLE z[]) {

    static INTEGER  i,j;
    static DOUBLE   su;

    /* d = j(trans) *np */
    
    for (i = 1 ; i <= ndual ; i++) {
        su = zero;
        for (j = 1 ; j <= ndual ; j++) {
            su = su+xj[j][i]*np[j];
        }
        ddual[i] = su;
    }
    /* computation of z */
    
    for (i = 1 ; i <= ndual ; i++) {
        z[i] = zero;
        for (j = iq+1 ; j <= ndual ; j++) {
            z[i] = z[i]+xj[i][j]*ddual[j];
        }
    }
    return;
}

/* **************************************************************************** */
/*                    compute correction of dual multipliers                    */
/* **************************************************************************** */
void o8rup(DOUBLE rv[]) {

    static DOUBLE   s;
    static INTEGER  i,j;

    for (i = iq ; i >= 1 ; i--) {
        s = zero;
        for (j = i+1 ; j <= iq ; j++) {
            s = s+r[i][j]*rv[j];
        }
        rv[i] = (ddual[i]-s)/r[i][i];
    }
    return;
}

/* **************************************************************************** */
/*                          delete constraint nr. l                             */
/* **************************************************************************** */
void o8dlcd(INTEGER ai[],INTEGER l) {
    
    DOUBLE o8dsq1(DOUBLE a,DOUBLE b);
    
    static INTEGER  qq,i,j,k;
    static DOUBLE   t1,t2,cc,ss,h,c1,s1,xny;

    for (i = 1 ; i <= iq ; i++) {
        if ( ai[i] == l ) {
            qq = i;
            
            goto L10;
         }
    }
    L10:

    for (i = qq ; i <= iq-1 ; i++) {
        ai[i]  = ai[i+1];
        ud1[i] = ud1[i+1];
        for (j = 1 ; j <= ndual ; j++) {
            r[j][i] = r[j][i+1];
        }
    }
    //    L20:

    ai[iq]    = ai[iq+1];
    ud1[iq]   = ud1[iq+1];
    ai[iq+1]  = 0;
    ud1[iq+1] = zero;
    for (j = 1 ; j <= iq ; j++) {
        r[j][iq] = zero;
    }
    iq = iq-1;
    
    if ( iq == 0 ) goto L100;
    
    for (j = qq ; j <= iq ; j++) {
        cc = r[j][j];
        ss = r[j+1][j];
        h  = o8dsq1(cc,ss);
        
        if ( h == zero ) goto L90;
        
        c1 = cc/h;
        s1 = ss/h;
        r[j+1][j] = zero;
        if ( c1 < zero ) {
            r[j][j] = -h;
            c1      = -c1;
            s1      = -s1;
        } else {
            r[j][j] = h;
        }
        xny = s1/(one+c1);
        for (k = j+1 ; k <= iq ; k++) {
            t1        = r[j][k];
            t2        = r[j+1][k];
            r[j][k]   = t1*c1+t2*s1;
            r[j+1][k] = xny*(t1+r[j][k])-t2;
        }
        for (k = 1 ; k <= ndual ; k++) {
            t1         = xj[k][j];
            t2         = xj[k][j+1];
            xj[k][j]   = t1*c1+t2*s1;
            xj[k][j+1] = xny*(xj[k][j]+t1)-t2;
        }
        L90:;
    }
    L100:

    rnorm = one;
    rlow  = one;
    
    /* in order to avoid a compiler error of hp in +op3 mode */
    
    if ( iq >= 1 ) {
        rnorm = myfabs(r[1][1]);
        rlow  = myfabs(r[1][1]);
        i     = 1;
        while ( i < iq ) {
            i     = i+1;
            rnorm = patMax(rnorm,myfabs(r[i][i]));
            rlow  = patMin(rlow, myfabs(r[i][i]));
        }
    }
    return;
}

/* **************************************************************************** */
/*                  add constraint whose gradient is given by np                */
/* **************************************************************************** */
void o8adcd(void) {
    
    DOUBLE o8dsq1(DOUBLE a,DOUBLE b);
                    
    static INTEGER  i,j,k;
    static DOUBLE   cc,ss,h,s1,c1,t1,t2,xny;

    for (j = ndual ; j >= iq+2 ; j--) {
        cc = ddual[j-1];
        ss = ddual[j];
        h  = o8dsq1(cc,ss);
        
        if ( h == zero ) goto L20;
        
        ddual[j] = zero;
        s1       = ss/h;
        c1       = cc/h;
        if ( c1 < zero ) {
            c1         = -c1;
            s1         = -s1;
            ddual[j-1] = -h;
        } else {
            ddual[j-1] = h;
        }
        xny = s1/(one+c1);
        for (k = 1 ; k <= ndual ; k++) {
            t1         = xj[k][j-1];
            t2         = xj[k][j];
            xj[k][j-1] = t1*c1+t2*s1;
            xj[k][j]   = xny*(t1+xj[k][j-1])-t2;
        }
        L20:;
    }
    iq = iq+1;
    for (i = 1 ; i <= iq ; i++) {
        r[i][iq] = ddual[i];
    }
    rnorm = one;
    rlow  = one;
    
    /* in order to avoid a compiler error of hp in +op3 mode */
    
    if ( iq >= 1 ) {
        rnorm = myfabs(r[1][1]);
        rlow  = myfabs(r[1][1]);
        i     = 1;
        while ( i < iq ) {
            i     = i+1;
            rnorm = patMax(rnorm,myfabs(r[i][i]));
            rlow  = patMin(rlow, myfabs(r[i][i]));
        }
    }
    return;
}

/* **************************************************************************** */
/*          computes the inverse of the upper triangular matrix part            */
/* **************************************************************************** */
void o8rinv(INTEGER n,DOUBLE a[][NX+1],INTEGER ndual,DOUBLE x[][NDUALM+1]) {
            
    /* computes the inverse of the upper triangular matrix part */
    /* of a and stores it in the upper triangle of the          */
    /* right lower minor of x                                   */
    /* actual dimension of a is n and of x ndual                */
    
    static INTEGER  l,j,k,incr;
    static DOUBLE   su;

    incr = ndual-n;
    
    /* incr = nr */
    
    for (j = n ; j >= 1 ; j--) {
    
        /* we assume a being sufficiently regular here. given in this */
        /* application. See top of o8qpdu                             */
        
        x[j+incr][j+incr] = one/a[j][j];
        for (k = j-1 ; k >= 1 ; k--) {
            su = zero;
            for (l = k+1 ; l <= j ; l++) {
                su = su+a[k][l]*x[l+incr][j+incr];
            }
            x[k+incr][j+incr] = -su/a[k][k];
        }
    }
    return;
}

/* **************************************************************************** */
/* suite of function interfaces, for performing scaling and                     */
/* external evaluations                                                         */
/* if bloc = TRUE then it is assumed that prior to calls to es<xyz>             */
/* valid function information is computed (via a call of user_eval)             */
/* and stored in fu and fugrad , setting valid = TRUE afterwards                */
/* **************************************************************************** */

/* **************************************************************************** */
/*                              objective function                              */
/* **************************************************************************** */
void esf(DOUBLE x[],DOUBLE *fx) {


    void ef(DOUBLE x[],DOUBLE *fx);
    
    static INTEGER  i;
    
    if ( bloc ) {
        if ( valid ) {
            *fx = fu[0];
            icf = icf+1;
        } else {
            fprintf(stderr,"donlp2: bloc-call, function info invalid\n");
            exit(1);
        }
    } else {
        for (i = 1 ; i <= n ; i++) {
            xtr[i] = x[i]*xsc[i];
        }
        ef(xtr,fx);
    }
    return;
}

/* **************************************************************************** */
/*                        gradient of objective function                        */
/* **************************************************************************** */
void esgradf(DOUBLE x[],DOUBLE gradf[]) {

    void    ef    (DOUBLE x[],DOUBLE *fx);
    void    egradf(DOUBLE x[],DOUBLE gradf[]);
    
    static INTEGER  j;
    static DOUBLE   d1,d2,d3,sd1,sd2,sd3,fhelp,fhelp1,fhelp2,
                    fhelp3,fhelp4,fhelp5,fhelp6,xincr,xhelp,floc;
                            
    if ( gunit[1][0] == 1 ) {
        for (j = 1 ; j <= n ; j++) {
            gradf[j] = zero;
        }
        gradf[gunit[2][0]] = gunit[3][0]*xsc[gunit[2][0]];
        
        return;
    }
    if ( bloc ) {
        if ( valid ) {
            icgf = icgf+1;
            for (j = 1 ; j <= n ; j++) {
                gradf[j] = xsc[j]*fugrad[j][0];
            }
            return;
        } else {
            fprintf(stderr,"donlp2: bloc call with function info invalid\n");
            exit(1);
        }
    } else {
        for (j = 1 ; j <= n ; j++) {
            xtr[j] = xsc[j]*x[j];
        }
        if ( analyt ) {
        
             egradf(xtr,gradf);
             
        } else {
            if ( difftype == 1 ) {
                deldif = patMin(tm1*sqrt(epsfcn),tm2);
                
                ef(xtr,&floc);
                
                for (j = 1 ; j <= n ; j++) {
                    xhelp = xtr[j];
                    xincr = patMin(patMin(tm2,deldif*myfabs(xhelp)+deldif),taubnd);
                    if ( xhelp >= zero ) {
                        xtr[j] = xhelp+xincr;
                    } else {
                        xtr[j] = xhelp-xincr;
                    }
                    ef(xtr,&fhelp);
                    
                    gradf[j] = (fhelp-floc)/(xtr[j]-xhelp);
                    xtr[j]   = xhelp;
                }
            } else if ( difftype == 2 ) {
                deldif = patMin(tm1*pow(epsfcn,one/three),tm2);
                for (j = 1 ; j <= n ; j++) {
                    xhelp  = xtr[j];
                    xincr  = patMin(patMin(tm2,deldif*myfabs(xhelp)+deldif),taubnd);
                    xtr[j] = xhelp+xincr;
                    
                    ef(xtr,&fhelp1);
                    
                    xtr[j] = xhelp-xincr;
                    
                    ef(xtr,&fhelp2);
                    
                    gradf[j] = (fhelp1-fhelp2)/(xincr+xincr);
                    xtr[j]   = xhelp;
                }
            } else {
                deldif = patMin(tm1*pow(epsfcn,one/seven),tm2);
                for (j = 1 ; j <= n ; j++) {
                    xhelp  = xtr[j];
                    xincr  = patMin(patMin(tm2,deldif*myfabs(xhelp)+deldif),taubnd/four);
                    xtr[j] = xhelp-xincr;
                    
                    ef(xtr,&fhelp1);
                    
                    xtr[j] = xhelp+xincr;
                    
                    ef(xtr,&fhelp2);
                    
                    xincr  = xincr+xincr;
                    d1     = xincr;
                    xtr[j] = xhelp-xincr;
                    
                    ef(xtr,&fhelp3);
                    
                    xtr[j] = xhelp+xincr;
                    
                    ef(xtr,&fhelp4);
                    
                    xincr  = xincr+xincr;
                    d2     = xincr;
                    xtr[j] = xhelp-xincr;
                    
                    ef(xtr,&fhelp5);
                    
                    xtr[j] = xhelp+xincr;
                    
                    ef(xtr,&fhelp6);
                    
                    xtr[j]   = xhelp;
                    d3       = xincr+xincr;
                    sd1      = (fhelp2-fhelp1)/d1;
                    sd2      = (fhelp4-fhelp3)/d2;
                    sd3      = (fhelp6-fhelp5)/d3;
                    sd3      = sd2-sd3;
                    sd2      = sd1-sd2;
                    sd3      = sd2-sd3;
                    gradf[j] = sd1+p4*sd2+sd3/c45;
                }
            }
        }
        for (j = 1 ; j <= n ; j++) {
            gradf[j] = xsc[j]*gradf[j];
        }
    }
    return;
}

/* **************************************************************************** */
/*              compute the i-th equality constaint, value is hxi               */
/* **************************************************************************** */
void esh(INTEGER i,DOUBLE x[],DOUBLE *hxi) {

    void eh(INTEGER i,DOUBLE x[],DOUBLE *hxi);
    
    static INTEGER  j;
    
    if ( bloc ) {
        if ( valid ) {
            if ( gunit[1][i] != 1 ) cres[i] = cres[i]+1;
            *hxi = fu[i];
            
            return;
            
        } else {
            fprintf(stderr,"donlp2: bloc-call with functoninformation invalid\n");
            exit(1);
        }
    } else {
        for (j = 1 ; j <= n ; j++) {
            xtr[j] = x[j]*xsc[j];
        }
        eh(i,xtr,hxi);
        
        return;
    }
}

/* **************************************************************************** */
/*              compute the gradient of the i-th equality constraint            */
/* **************************************************************************** */
void esgradh(INTEGER i,DOUBLE x[],DOUBLE gradhi[]) {

    void    eh    (INTEGER i,DOUBLE x[],DOUBLE *hxi);
    void    egradh(INTEGER i,DOUBLE x[],DOUBLE gradhi[]);
    
    static INTEGER  j;
    static DOUBLE   d1,d2,d3,sd1,sd2,sd3,fhelp,fhelp1,fhelp2,
                    fhelp3,fhelp4,fhelp5,fhelp6,xincr,xhelp,hixloc;
         
    if ( gunit[1][i] == 1 ) {
        for (j = 1 ; j <= n ; j++) {
            gradhi[j] = zero;
        }
        gradhi[gunit[2][i]] = gunit[3][i]*xsc[gunit[2][i]];
        
        return;
    }
    if ( bloc ) {
        if ( valid ) {
            cgres[i] = cgres[i]+1;
            for (j = 1 ; j <= n ; j++) {
                gradhi[j] = xsc[j]*fugrad[j][i];
            }
            return;
            
        } else {
            fprintf(stderr,"donlp2: bloc call with function info invalid\n");
            exit(1);
        }
    } else {
        for (j = 1 ; j <= n ; j++) {
            xtr[j] = xsc[j]*x[j];
        }
        if ( analyt ) {
        
            egradh(i,xtr,gradhi);
            
        } else {
            if ( difftype == 1 ) {
                deldif = patMin(tm1*sqrt(epsfcn),tm2);
                
                eh(i,xtr,&hixloc);
                
                for (j = 1 ; j <= n ; j++) {
                    xhelp = xtr[j];
                    xincr = patMin(patMin(tm2,deldif*myfabs(xhelp)+deldif),taubnd);
                    if ( xhelp >= zero ) {
                        xtr[j] = xhelp+xincr;
                    } else {
                        xtr[j] = xhelp-xincr;
                    }
                    eh(i,xtr,&fhelp);
                    
                    gradhi[j] = (fhelp-hixloc)/(xtr[j]-xhelp);
                    xtr[j]    = xhelp;
                }
            } else if ( difftype == 2 ) {
                deldif = patMin(tm1*pow(epsfcn,one/three),tm2);
                for (j = 1 ; j <= n ; j++) {
                    xhelp  = xtr[j];
                    xincr  = patMin(patMin(tm2,deldif*myfabs(xhelp)+deldif),taubnd);
                    xtr[j] = xhelp+xincr;
                    
                    eh(i,xtr,&fhelp1);
                    
                    xtr[j] = xhelp-xincr;
                    
                    eh(i,xtr,&fhelp2);
                    
                    gradhi[j] = (fhelp1-fhelp2)/(xincr+xincr);
                    xtr[j]    = xhelp;
                }
            } else {
                deldif = patMin(tm1*pow(epsfcn,one/seven),tm2);
                for (j = 1 ; j <= n ; j++) {
                    xhelp  = xtr[j];
                    xincr  = patMin(patMin(tm2,deldif*myfabs(xhelp)+deldif),taubnd/four);
                    xtr[j] = xhelp-xincr;
                    
                    eh(i,xtr,&fhelp1);
                    
                    xtr[j] = xhelp+xincr;
                    
                    eh(i,xtr,&fhelp2);
                    
                    xincr  = xincr+xincr;
                    d1     = xincr;
                    xtr[j] = xhelp-xincr;
                    
                    eh(i,xtr,&fhelp3);
                    
                    xtr[j] = xhelp+xincr;
                    
                    eh(i,xtr,&fhelp4);
                    
                    xincr  = xincr+xincr;
                    d2     = xincr;
                    xtr[j] = xhelp-xincr;
                    
                    eh(i,xtr,&fhelp5);
                    
                    xtr[j] = xhelp+xincr;
                    
                    eh(i,xtr,&fhelp6);
                    
                    xtr[j]    = xhelp;
                    d3        = xincr+xincr;
                    sd1       = (fhelp2-fhelp1)/d1;
                    sd2       = (fhelp4-fhelp3)/d2;
                    sd3       = (fhelp6-fhelp5)/d3;
                    sd3       = sd2-sd3;
                    sd2       = sd1-sd2;
                    sd3       = sd2-sd3;
                    gradhi[j] = sd1+p4*sd2+sd3/c45;
                }
            }
        }
        for (j = 1 ; j <= n ; j++) {
            gradhi[j] = xsc[j]*gradhi[j];
        }
        return;
    }
}

/* **************************************************************************** */
/*              compute the i-th inequality constaint, bounds included          */
/* **************************************************************************** */
void esg(INTEGER i,DOUBLE x[],DOUBLE *gxi) {

    void eg(INTEGER i,DOUBLE x[],DOUBLE *gxi);

    static INTEGER  j;
    
    if ( bloc ) {
        if ( valid ) {
            if ( gunit[1][i+nh] != 1 ) cres[i+nh] = cres[i+nh]+1;
            *gxi = fu[i+nh];
        } else {
            fprintf(stderr,"donlp2: bloc call with function info invalid\n");
            exit(1);
        }
    } else {
        for (j = 1 ; j <= n ; j++) {
            xtr[j] = x[j]*xsc[j];
        }
        eg(i,xtr,gxi);
    }
    return;
}

/* **************************************************************************** */
/*                  evaluation of gradient of i-th inequality                   */
/* **************************************************************************** */
void esgradg(INTEGER i,DOUBLE x[],DOUBLE gradgi[]) {

    /* evaluation of gradient of i-th inequality                    */
    /* not necessary for bounds, but constant gradients must be set */
    /* here e.g. using dcopy from a data-field                      */

    void    eg    (INTEGER i,DOUBLE x[],DOUBLE *gxi);
    void    egradg(INTEGER i,DOUBLE x[],DOUBLE gradgi[]);
    
    static INTEGER  j;
    static DOUBLE   d1,d2,d3,sd1,sd2,sd3,fhelp,fhelp1,fhelp2,
                    fhelp3,fhelp4,fhelp5,fhelp6,xincr,xhelp,gixloc;
                            
    if ( gunit[1][i+nh] == 1 ) {
        for (j = 1 ; j <= n ; j++) {
            gradgi[j] = zero;
        }
        gradgi[gunit[2][i+nh]] = gunit[3][i+nh]*xsc[gunit[2][i+nh]];
        
        return;
    }
    if ( bloc ) {
        if ( valid ) {
            if ( gunit[1][i+nh] != 1 ) cgres[i+nh] = cgres[i+nh]+1;
            for (j = 1 ; j <= n ; j++) {
                gradgi[j] = xsc[j]*fugrad[j][i+nh];
            }
            return;
            
        } else {
            fprintf(stderr,"donlp2: bloc call with function info invalid\n");
            exit(1);
        }
    } else {
        for (j = 1 ; j <= n ; j++) {
            xtr[j] = x[j]*xsc[j];
        }
        if ( analyt ) {
        
            egradg(i,xtr,gradgi);
            
        } else {
            if ( difftype == 1 ) {
                deldif = patMin(tm1*sqrt(epsfcn),tm2);
                
                eg(i,xtr,&gixloc);
                
                for (j = 1 ; j <= n ; j++) {
                    xhelp = xtr[j];
                    xincr = patMin(patMin(tm2,deldif*myfabs(xhelp)+deldif),taubnd);
                    if ( xhelp >= zero ) {
                        xtr[j] = xhelp+xincr;
                    } else {
                        xtr[j] = xhelp-xincr;
                    }
                    eg(i,xtr,&fhelp);
                    
                    gradgi[j] = (fhelp-gixloc)/(xtr[j]-xhelp);
                    xtr[j]    = xhelp;
                }
            } else if ( difftype == 2 ) {
                deldif = patMin(tm1*pow(epsfcn,one/three),tm2);
                for (j = 1 ; j <= n ; j++) {
                    xhelp  = xtr[j];
                    xincr  = patMin(patMin(tm2,deldif*myfabs(xhelp)+deldif),taubnd);
                    xtr[j] = xhelp+xincr;
                    
                    eg(i,xtr,&fhelp1);
                    
                    xtr[j] = xhelp-xincr;
                    
                    eg(i,xtr,&fhelp2);
                    
                    gradgi[j] = (fhelp1-fhelp2)/(xincr+xincr);
                    xtr[j]    = xhelp;
                }
            } else {
                deldif = patMin(tm1*pow(epsfcn,one/seven),tm2);
                for (j = 1 ; j <= n ; j++) {
                    xhelp  = xtr[j];
                    xincr  = patMin(patMin(tm2,deldif*myfabs(xhelp)+deldif),taubnd/four);
                    xtr[j] = xhelp-xincr;
                    
                    eg(i,xtr,&fhelp1);
                    
                    xtr[j] = xhelp+xincr;
                    
                    eg(i,xtr,&fhelp2);
                    
                    xincr  = xincr+xincr;
                    d1     = xincr;
                    xtr[j] = xhelp-xincr;
                    
                    eg(i,xtr,&fhelp3);
                    
                    xtr[j] = xhelp+xincr;
                    
                    eg(i,xtr,&fhelp4);
                    
                    xincr  = xincr+xincr;
                    d2     = xincr;
                    xtr[j] = xhelp-xincr;
                    
                    eg(i,xtr,&fhelp5);
                    
                    xtr[j] = xhelp+xincr;
                    
                    eg(i,xtr,&fhelp6);
                    
                    xtr[j]    = xhelp;
                    d3        = xincr+xincr;
                    sd1       = (fhelp2-fhelp1)/d1;
                    sd2       = (fhelp4-fhelp3)/d2;
                    sd3       = (fhelp6-fhelp5)/d3;
                    sd3       = sd2-sd3;
                    sd2       = sd1-sd2;
                    sd3       = sd2-sd3;
                    gradgi[j] = sd1+p4*sd2+sd3/c45;
                }
            }
        }
        for (j = 1 ; j <= n ; j++) {
            gradgi[j] = xsc[j]*gradgi[j];
        }
    }
    return;
}


