//-*-c++-*------------------------------------------------------------
//
// File name : patDonlp2.h
// Author :    \URL[Michel Bierlaire]{http://rosowww.epfl.ch/mbi}
// Date :      Tue Mar 19 17:17:15 2002
//
//--------------------------------------------------------------------

#ifndef patDonlp2_h
#define patDonlp2_h 

/**
@doc This class defines a C++ interface to the Donlp2 routine, by
P. Spellucci. Actually, it is based on a translation of the f77-version of
donlp2.f as of 7/99 to ANSI-C done by Serge Schoeffler. 
@author \URL[Michel Bierlaire]{http://rosowww.epfl.ch/mbi}, EPFL (Tue Mar 19 17:17:15 2002) 

*/


#include "o8para.h"
#include "patError.h"
#include "trNonLinearAlgo.h"
#include "donlp2Parameters.h"
class patNonLinearProblem ;

class patDonlp2 : public trNonLinearAlgo {

public:
  /**
   */
  patDonlp2(donlp2Parameters theParameters, patNonLinearProblem* aProblem = NULL) ;
  /**
   */
  virtual ~patDonlp2() ;
  /**
   */
  void setProblem(patNonLinearProblem* aProblem) ;

  /**
   */
  patNonLinearProblem* getProblem() ;

  /**
   */
  void defineStartingPoint(const patVariables& x0) ;

  /**
   */
  patVariables getStartingPoint() ;
  /**
   */
  patVariables getSolution(patError*& err) ;

  /**
   */
  patReal getValueSolution(patError*& err) ;

  /**
   */
  patVariables getLowerBoundsLambda() ;
  /**
   */
  patVariables getUpperBoundsLambda() ;
  /**
   */
  patVariables getNonLinIneqLambda() ;
  /**
   */
  patVariables getLinIneqLambda() ;
  /**
   */
  patVariables getNonLinEqLambda() ;
  /**
   */
  patVariables getLinEqLambda() ;

  /**
   */
  patString getName() ;

  /**
     @return Diagnostic from cfsqp
   */
  patString run(patError*& err) ;
  /**
     @return number of iterations. If there is any error, 0 is returned. 
   */
  patULong nbrIter() ;


private:

  patVariables startingPoint ;
  patVariables solution ;
  patVariables lowerBoundsLambda ;
  patVariables upperBoundsLambda ;
  patVariables nonLinIneqLambda ;
  patVariables linIneqLambda ;
  patVariables nonLinEqLambda ;
  patVariables linEqLambda ;

  donlp2Parameters theParameters ;
};

/**
 @doc User functions required by donlp2. See User's guide.
 */
void setup(void) ;

/**
 @doc User functions required by donlp2. See User's guide.
 */
void eval_extern(INTEGER mode) ;


/**
 @doc User functions required by donlp2. See User's guide.
 */
void setup0(void) ;


/**
 @doc User functions required by donlp2. See User's guide.
 */
void solchk(void) ;
/**
 @doc User functions required by donlp2. Objective function
 */
void ef(DOUBLE x[],DOUBLE *fx) ;
/**
 @doc User functions required by donlp2. Gradient of objective function
 */
void egradf(DOUBLE x[],DOUBLE gradf[]) ;
/**
 @doc User functions required by donlp2. Compute the i-th equality constaint, value is hxi
 */
void eh(INTEGER i,DOUBLE x[],DOUBLE *hxi) ;

/**
 @doc User functions required by donlp2. Compute the gradient of the i-th equality constraint
 */
void egradh(INTEGER i,DOUBLE x[],DOUBLE gradhi[]) ;

/**
 @doc User functions required by donlp2. Compute the i-th inequality constaint, bounds included.
 */
void eg(INTEGER i,DOUBLE x[],DOUBLE *gxi) ;

/**
 @doc User functions required by donlp2. See User's guide. Compute the gradient of the i-th inequality constraint
 */
void egradg(INTEGER i,DOUBLE x[],DOUBLE gradgi[]) ;







#endif
