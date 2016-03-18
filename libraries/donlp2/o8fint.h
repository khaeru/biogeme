/* **************************************************************************** */
/*                                   o8fint.h                                   */
/* **************************************************************************** */

/* if bloc = TRUE then it is assumed that functionevaluation takes place        */      
/* at once in an external routine and                                           */
/* that user_eval has been called before calling for evaluation of functions    */
/* the latter then simply consists in copying data                              */
/* to donlp2's own data                                                         */
/* user_eval must set valid = TRUE, if functionvalues are valid for the         */
/* current xtr                                                                  */
/* corr is set to true by donlp2, if the initial x does not satisfy             */
/* the bound constraints. x is modified in this case                            */
/* difftype = 1,2,3 numerical differentiation by the ordinary forward           */
/* differences, by central differences or by Richardson-extrapolation           */
/* of order 6, requiring n, 2n , 6n additional function evaluations             */
/* respectively                                                                 */
/* epsfcn is the assumed precision of the function evaluation, to be            */
/* set by the user                                                              */
/* taubnd: amount by which bound constraints may be violated during             */
/* finite differencing, set by the user                                         */

X LOGICAL   bloc,valid,corr;
X INTEGER   difftype;
X DOUBLE    xtr[NX+1],xsc[NX+1],fu[NRESM+1],fugrad[NX+1][NRESM+1],
            fud[NRESM+1][7],epsfcn,taubnd;

