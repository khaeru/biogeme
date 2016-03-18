//-*-c++-*------------------------------------------------------------
//
// File name : patTestNormal.h
// Author :    \URL[Michel Bierlaire]{http://roso.epfl.ch/mbi}
// Date :      Wed Jun 15 16:38:22 2005
//
//--------------------------------------------------------------------

#ifndef patTestNormal_h
#define patTestNormal_h

#include "trFunction.h"
#include "patNormalCdf.h"

/**
 */
class patTestNormal : public  trFunction {
public:
  /**
   */
  patReal getFunction(trVector* x,
		      patBoolean* success,
		      patError*& err)  ;
  /**
   */
  trVector* getGradient(trVector* x,
			trVector* grad,
			patBoolean* success,
			patError*& err) ;


  /**
  */
   trVector getHessianTimesVector(trVector* x,
				  trVector* v,
				  patBoolean* success,
				  patError*& err)
     const ;
    
  
  /**
   */
   trHessian* computeHessian(trVector* x,
			     trHessian& hessian,
			     patBoolean* success,
			     patError*& err) ;


  /**
  */
   trHessian* getCheapHessian(trVector* x,
			      trHessian& hessian,
			      patBoolean* success,
			      patError*& err) ;

  /**
   */
  virtual patBoolean isCheapHessianAvailable()  ;

  /**
   */
  virtual patBoolean isGradientAvailable(patError*& err) const ;
  /**
   */
  virtual patBoolean isHessianTimesVectorAvailable(patError*& err) const ;

  /**
   */
  virtual patBoolean isHessianAvailable(patError*& err) const ;


  /**
   */
  unsigned long getDimension(patError*& err) const  ;

};

#endif

