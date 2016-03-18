//-*-c++-*------------------------------------------------------------
//
// File name : patDonlp2.cc
// Author :    Michel Bierlaire
// Date :      Tue Mar 19 17:27:03 2002
//
//--------------------------------------------------------------------

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "patDonlp2.h"
#include <sstream>
#include <numeric>
#include <cstring>
#include "trFunction.h"
#include "patMath.h"
#include "patErrMiscError.h"
#include "patErrNullPointer.h"
#include "patNonLinearProblem.h"
#include "donlp2Parameters.h"

/**
   @doc global pointer to the problem for the functions used by solvOpt
*/
patNonLinearProblem* globalDonlp2Problem ;

/**
   @doc global pointer to the starting point
*/
patVariables* globalStartingPoint ;

/**
 */
patString donlp2(donlp2Parameters theParameters,patIterationBackup* theBackup) ;

patDonlp2::patDonlp2(donlp2Parameters p, patNonLinearProblem* aProblem) :
  trNonLinearAlgo(aProblem),
  startingPoint((aProblem==NULL)?0:aProblem->nVariables()),
  solution((aProblem==NULL)?0:aProblem->nVariables()),
  lowerBoundsLambda((aProblem==NULL)?0:aProblem->nVariables()),
  upperBoundsLambda((aProblem==NULL)?0:aProblem->nVariables()),
  nonLinIneqLambda((aProblem==NULL)?0:aProblem->nNonLinearIneq()),
  linIneqLambda((aProblem==NULL)?0:aProblem->nLinearIneq()),
  nonLinEqLambda((aProblem==NULL)?0:aProblem->nNonLinearEq()),
  linEqLambda((aProblem==NULL)?0:aProblem->nLinearEq()),
  theParameters(p)
{
  globalDonlp2Problem = aProblem ;
  globalStartingPoint = &startingPoint ;
}

patDonlp2::~patDonlp2() {

}

void patDonlp2::setProblem(patNonLinearProblem* aProblem) {
  theProblem = aProblem ;
  if (theProblem != NULL) {
    solution.resize(theProblem->nVariables()) ;
    lowerBoundsLambda.resize(theProblem->nVariables()) ;
    upperBoundsLambda.resize(theProblem->nVariables()) ;
    nonLinIneqLambda.resize(theProblem->nNonLinearIneq()) ;
    linIneqLambda.resize(theProblem->nLinearIneq()) ;
    nonLinEqLambda.resize(theProblem->nNonLinearEq()) ;
    linEqLambda.resize(theProblem->nLinearEq()) ;
  }
}

patNonLinearProblem* patDonlp2::getProblem() {
  return theProblem ;
}

void patDonlp2::defineStartingPoint(const patVariables& x0) {
  startingPoint=x0 ;
}

patVariables patDonlp2::getStartingPoint() {
  return startingPoint ;
}

patVariables patDonlp2::getSolution(patError*& err) {
  return solution ;
}

patReal patDonlp2::getValueSolution(patError*& err) {
  patBoolean success = patTRUE ;
  if (theProblem == NULL) {
    err = new patErrNullPointer("patNonLinearProblem") ;
    WARNING(err->describe()) ;
    return patReal();

  }
  trFunction* f = theProblem->getObjective(err) ;
  if (err != NULL) {
    WARNING(err->describe()) ;
    return patReal() ;
  }
  if (f == NULL) {
    err = new patErrNullPointer("trFunction") ;
    WARNING(err->describe()) ;
    return patReal();
  }

  patReal result = f->computeFunction(&solution,&success,err) ;
  if (err != NULL) {
    WARNING(err->describe()) ;
    return patReal() ;
  }
  if (!success) {
    WARNING("Problem in function evaluation") ;
  }
  return result ;
}


patVariables patDonlp2::getLowerBoundsLambda() {
  return lowerBoundsLambda ;
}
patVariables patDonlp2::getUpperBoundsLambda() {
  return upperBoundsLambda ;
}
patVariables patDonlp2::getNonLinIneqLambda() {
  return nonLinIneqLambda ;
}

patVariables patDonlp2::getLinIneqLambda() {
  return linIneqLambda ;
}
patVariables patDonlp2::getNonLinEqLambda() {
  return nonLinEqLambda ;
}
patVariables patDonlp2::getLinEqLambda() {
  return linEqLambda ;
}



patString patDonlp2::run(patError*& err) {
#define  X extern
#include "o8fuco.h"
#include "o8comm.h"
#undef   X

  if (theProblem == NULL) {
    patString diag("No problem has been defined") ;
    return diag ;
  }

  DEBUG_MESSAGE("Call to donlp2") ;
  patString status = donlp2(theParameters,theBackup) ;
  DEBUG_MESSAGE("Done...") ;
  for (unsigned k = 0 ; k < theProblem->nVariables() ; ++k) {
    solution[k] = x[k+1] ;
  }

  theProblem->computeLagrangeMultipliers(solution,err) ;
  if (err != NULL) {
    WARNING(err->describe()) ;
    DELETE_PTR(err) ;
  }
  return status ;
}


patString patDonlp2::getName() {
  return patString("DONLP2") ;
}

void setup(donlp2Parameters theParameters) {
#define  X extern
#include "o8comm.h"
#undef   X


  epsx = theParameters.epsx ;
  delmin = theParameters.delmin ;
  smallw = theParameters.smallw ;
  epsdif = theParameters.epsdif ;
  nreset = theParameters.nreset ;

  DEBUG_MESSAGE("epsx = " << epsx) ;
  DEBUG_MESSAGE("delmin = " << delmin) ;
  DEBUG_MESSAGE("smallw = " << smallw) ;
  DEBUG_MESSAGE("epsdif = " << epsdif) ;
  DEBUG_MESSAGE("nreset = " << nreset) ;

  te2 = TRUE ;
  return ;
}
void setup0(void) {
#define  X extern
#include "o8fuco.h"
#include "o8comm.h"
#undef   X


  
  strncpy(name,globalDonlp2Problem->getProblemName().c_str(),41);
  n    = int(globalDonlp2Problem->nVariables()) ;
  nh   = int(globalDonlp2Problem->nNonLinearEq() + globalDonlp2Problem->nLinearEq()) ;
    ng   = 
      int(globalDonlp2Problem->nNonLinearIneq() + 
	  globalDonlp2Problem->nLinearIneq() + 
	  2*globalDonlp2Problem->nVariables()) ;

    if (globalStartingPoint != NULL) {
      for (int i = 0 ;
	   i < n ;
	   ++i) {
	x[i+1] = (*globalStartingPoint)[i] ;
      }
    }
    else {
      for (int i = 0 ;
	   i < n ;
	   ++i) {
	x[i+1] = 0.0 ;
      }
    }

    
    
    del0 = 1.e0;
    tau0 = 1.e0;
    tau  = .1e0;
    
    silent = FALSE;
    analyt = TRUE;
    epsdif = 0.e0;

    // The vector gunit describes the nature of the inequality constraints.

    // I don't understand why it must start at nh, but it seems to be necessary

    for (unsigned short i = 0 ;
	 i <= nh;
	 ++i) {
        gunit[1][i] = -1;
        gunit[2][i] = 0;
        gunit[3][i] = 0;      
    }

      unsigned short index = nh+1 ;
      for (unsigned short i = 0 ;
  	 i < globalDonlp2Problem->nNonLinearIneq() ;
  	 ++i) {
          gunit[1][index] = -1;
          gunit[2][index] = 0;
          gunit[3][index] = 0;      
  	++index ;
      }
      for (unsigned short i = 0 ;
  	 i < globalDonlp2Problem->nLinearIneq() ;
  	 ++i) {
          gunit[1][index] = -1;
          gunit[2][index] = 0;
          gunit[3][index] = 0;      
  	++index ;
      }
     for (unsigned short i = 0 ;
  	 i < globalDonlp2Problem->nVariables() ;
  	 ++i) {
        // Lower bounds
          gunit[1][index] = 1;
          gunit[2][index] = i+1;
          gunit[3][index] = 1;      
  	++index ;
      }
      for (unsigned short i = 0 ;
  	 i < globalDonlp2Problem->nVariables() ;
  	 ++i) {
        // Upper bounds
          gunit[1][index] = 1;
          gunit[2][index] = i+1;
          gunit[3][index] = -1;      
  	++index ;
      }
 
    return;

}

void solchk(void) {

}

void ef(DOUBLE x[],DOUBLE *fx) {
#define  X extern
#include "o8fuco.h"
#undef   X

  trVector xstl(n) ;
  for (unsigned short i = 0 ; i < n ; ++i) {
    xstl[i] = x[i+1] ;
  }
  patError* err = NULL ;
  patBoolean success = patTRUE ;
  trFunction* f = globalDonlp2Problem->getObjective(err) ;
  if (err != NULL) {
    WARNING(err->describe()) ;
    *fx = patMaxReal ;
    return ;
  }
  if (f == NULL) {
    WARNING("Null pointer");
    *fx = patMaxReal ;
    return ;
  }

  *fx = f->computeFunction(&xstl,&success,err) ;
  if (err != NULL) {
    WARNING(err->describe()) ;
    *fx = patMaxReal ;
    return ;
  }
  if (!success) {
    *fx = patMaxReal ;
    return ;
  }
  return ;


}

void egradf(DOUBLE x[],DOUBLE gradf[]) {
#define  X extern
#include "o8fuco.h"
#undef   X
  trVector xstl(n) ;
  for (unsigned short i = 0 ; i < n ; ++i) {
    xstl[i] = x[i+1] ;
  }
  patError* err = NULL ;
  patBoolean success = patTRUE ;
  trFunction* f = globalDonlp2Problem->getObjective(err) ;
  if (err != NULL) {
    WARNING(err->describe()) ;
    return ;
  }
  if (f == NULL) {
    WARNING("Null pointer");
    return ;
  }

  trVector grad(xstl.size());
  f->computeFunctionAndDerivatives(&xstl,&grad,NULL,&success,err) ;
  if (err != NULL) {
    WARNING(err->describe()) ;
    return ;
  }
  if (!success) {
    WARNING("Error in function computation") ;
    return ;
  }

  for (unsigned short i = 0 ; i < n ; ++i) {
    gradf[i+1] = grad[i] ;
  }
  return ;

}

void eh(INTEGER i,DOUBLE x[],DOUBLE *hxi) {

#define  X extern
#include "o8fuco.h"
#undef   X

  ++cres[i] ;

  trVector xstl(n) ;
  for (unsigned short k = 0 ; k < n ; ++k) {
    xstl[k] = x[k+1] ;
  }
  patError* err = NULL ;
  patBoolean success = patTRUE ;

  // i is numbered from 1

  unsigned short j = i-1 ;
  if (j < globalDonlp2Problem->nNonLinearEq()) {
    // This is a nonlinear constraint
    trFunction* f = globalDonlp2Problem->getNonLinEquality(j,err) ;
    if (err != NULL) {
      WARNING(err->describe()) ;
      *hxi = patMaxReal ;
      return ;
    }
    if (f == NULL) {
      WARNING("Null pointer");
      *hxi = patMaxReal ;
      return ;
    }
    
    *hxi = f->computeFunction(&xstl,&success,err) ;
    if (err != NULL) {
      WARNING(err->describe()) ;
      *hxi = patMaxReal ;
      return ;
    }
    if (!success) {
      *hxi = patMaxReal ;
      return ;
    }
    return ;
    
  }
  else {
    // This is a linear constraint
    j -= globalDonlp2Problem->nNonLinearEq() ;
    pair<patVariables,patReal> f = globalDonlp2Problem->getLinEquality(j,err) ;
    if (err != NULL) {
      WARNING(err->describe()) ;
      *hxi = patMaxReal ;
      return ;
    }
    *hxi = inner_product(f.first.begin(),f.first.end(),xstl.begin(),-f.second) ;
    return ;

  }
}

void egradh(INTEGER i,DOUBLE x[],DOUBLE gradhi[]) {
#define  X extern
#include "o8fuco.h"
#undef   X

  trVector xstl(n) ;
  for (unsigned short k = 0 ; k < n ; ++k) {
    xstl[k] = x[k+1] ;
  }
  patError* err = NULL ;
  patBoolean success = patTRUE ;

  // i is numbered from 1

  unsigned short j = i-1 ;
  if (j < globalDonlp2Problem->nNonLinearEq()) {
    // This is a nonlinear constraint
    trFunction* f = globalDonlp2Problem->getNonLinEquality(j,err) ;
    if (err != NULL) {
      WARNING(err->describe()) ;
      return ;
    }
    if (f == NULL) {
      WARNING("Null pointer");
      return ;
    }
    trVector grad(xstl.size());
    f->computeFunctionAndDerivatives(&xstl,&grad,NULL,&success,err) ;
    if (err != NULL) {
      WARNING(err->describe()) ;
      return ;
    }
    if (!success) {
      return ;
    }
    for (unsigned short k = 0 ; k < n ; ++k) {
      gradhi[k+1] = grad[k] ;
    }
    return ;
  }
  else {
    // This is a linear constraint
    j -= globalDonlp2Problem->nNonLinearEq() ;
    pair<patVariables,patReal> f = globalDonlp2Problem->getLinEquality(j,err) ;
    if (err != NULL) {
      WARNING(err->describe()) ;
      return ;
    }
    for (unsigned short k = 0 ; k < n ; ++k) {
      gradhi[k+1] = f.first[k] ;
    }
    return ;
  }
}

void eg(INTEGER i,DOUBLE x[],DOUBLE *gxi) {

  // WARNING: the patNonLinearProblem defines inequalities as g(x) <= 0, while
  // donlp2 defines iniequalities as g(x) >=0. The sign must therefore be
  // changed.

#define  X extern
#include "o8fuco.h"
#undef   X

  if ( gunit[1][i+nh] == -1 ) ++cres[i+nh] ;
    
  patError* err = NULL ;
  static patVariables lowerBound = globalDonlp2Problem->getLowerBounds(err) ;
  if (err != NULL) {
    WARNING(err->describe()) ;
    *gxi = -patMaxReal ;
    return ;
  }
  static patVariables upperBound = globalDonlp2Problem->getUpperBounds(err) ;
  if (err != NULL) {
    WARNING(err->describe()) ;
    *gxi = -patMaxReal ;
    return ;
  }
  
  trVector xstl(n) ;
  for (unsigned short k = 0 ; k < n ; ++k) {
    xstl[k] = x[k+1] ;
  }

  patBoolean success = patTRUE ;

  // i is numbered from 1

  unsigned short j = i-1 ;
  if (j < globalDonlp2Problem->nNonLinearIneq()) {
    // This is a nonlinear constraint
    trFunction* f = globalDonlp2Problem->getNonLinInequality(j,err) ;
    if (err != NULL) {
      WARNING(err->describe()) ;
      *gxi = -patMaxReal ;
      return ;
    }
    if (f == NULL) {
      WARNING("Null pointer");
      *gxi = -patMaxReal ;
      return ;
    }
    
    *gxi = -f->computeFunction(&xstl,&success,err) ;
    if (err != NULL) {
      WARNING(err->describe()) ;
      *gxi = -patMaxReal ;
      return ;
    }
    if (!success) {
      *gxi = -patMaxReal ;
      return ;
    }
    return ;
  }
  j -= globalDonlp2Problem->nNonLinearIneq() ;
  if (j < globalDonlp2Problem->nLinearIneq()) {
    // This is a linear constraint
    pair<patVariables,patReal> f = globalDonlp2Problem->getLinInequality(j,err) ;
    if (err != NULL) {
      WARNING(err->describe()) ;
      *gxi = -patMaxReal ;
      return ;
     }
    *gxi = -inner_product(f.first.begin(),
			 f.first.end(),
			 xstl.begin(),
			 -f.second) ;
    return ;
  }
  j -= globalDonlp2Problem->nLinearIneq() ;
  
  if ( j < globalDonlp2Problem->nVariables()) {
    // This is a lower bound
    *gxi = -lowerBound[j]+xstl[j] ;
    return ;
  }
  j -= globalDonlp2Problem->nVariables() ;
// This is an upper bound
  *gxi = -xstl[j] + upperBound[j] ;
  return ;
}

void egradg(INTEGER i,DOUBLE x[],DOUBLE gradgi[]) {
#define  X extern
#include "o8fuco.h"
#undef   X

  if ( gunit[1][i+nh] == 1 ) return;
  trVector xstl(n) ;
  for (unsigned short k = 0 ; k < n ; ++k) {
    xstl[k] = x[k+1] ;
  }
  patError* err = NULL ;
  patBoolean success = patTRUE ;

  // i is numbered from 1

  unsigned short j = i-1 ;
  if (j < globalDonlp2Problem->nNonLinearIneq()) {
    // This is a nonlinear constraint
    trFunction* f = globalDonlp2Problem->getNonLinInequality(j,err) ;
    if (err != NULL) {
      WARNING(err->describe()) ;
      return ;
    }
    if (f == NULL) {
      WARNING("Null pointer");
      return ;
    }
    trVector grad(xstl.size()) ;
    f->computeFunctionAndDerivatives(&xstl,&grad,NULL,&success,err) ;
    if (err != NULL) {
      WARNING(err->describe()) ;
      return ;
    }
    grad = -grad ;
    if (!success) {
      return ;
    }
    for (unsigned short k = 0 ; k < n ; ++k) {
      gradgi[k+1] = grad[k] ;
    }
    return ;
  }
  else  {
    // This is a linear constraint
    if (globalDonlp2Problem->nLinearIneq() > 0) {
      j -= globalDonlp2Problem->nNonLinearIneq() ;
      pair<patVariables,patReal> f = globalDonlp2Problem->getLinInequality(j,err) ;
      if (err != NULL) {
	WARNING(err->describe()) ;
	return ;
      }
      patVariables grad = -f.first ;
      for (unsigned short k = 0 ; k < n ; ++k) {
	gradgi[k+1] = grad[k] ;
      }
      return ;
    }
    
  }
  // If it is a bound, no need to provide the gradient...

//    j -= globalDonlp2Problem->nLinearIneq() ;
//    if (j < n) {
//      // This is a lower bound
//        for (unsigned short k = 0 ; k < n ; ++k) {
//  	gradgi[k+1] = 0.0 ;
//        }
//        gradgi[j+i] = 1.0 ;
//    }
//    else {
//      // This is an upper bound
//        for (unsigned short k = 0 ; k < n ; ++k) {
//  	gradgi[k+1] = 0.0 ;
//        }
//        gradgi[j+i-n] = -1.0 ;
//    }
}

void eval_extern(INTEGER mode) {
#define  X extern
#include "o8comm.h"
#include "o8fint.h"
#undef   X
#include "o8cons.h"

    return;
}


patULong patDonlp2::nbrIter() {
  // NO YET IMPLEMENTED.
  return 0 ;
}
