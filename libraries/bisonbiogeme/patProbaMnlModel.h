//-*-c++-*------------------------------------------------------------
//
// File name : patProbaMnlModel.h
// Author :    \URL[Michel Bierlaire]{http://rosowww.epfl.ch/mbi}
// Date :      Tue Aug 26 14:40:41 2003
//
//--------------------------------------------------------------------

#ifndef patProbaMnlModel_h
#define patProbaMnlModel_h

#include "patError.h"

class patGEV ;
class patUtility ;

#include "patProbaModel.h"

/**
   @doc Defines the probability model for multinomial logit model, based on the specification of the utility functions and a sample
   @author \URL[Michel Bierlaire]{http://rosowww.epfl.ch/mbi}, EPFL ( Tue Aug 26 14:40:41 2003)
 */

class patProbaMnlModel : public patProbaModel {

public:
  /**
     Sole constructor
   */
  patProbaMnlModel(patUtility* aUtility) ;

  /**
     Destructor
  */
  ~patProbaMnlModel() ;

  /**
     Evaluates the logarithm of the probability given by the model that the
     individiual chooses alternative 'index', knowing the utilities, for a given draw number. If
     requested, the derivatives are evaluated as well. Note that the derivatives
     are cumulated to previous values. In order to have the value of the
     derivatives, the corresponding storage area must be initialized to zero
     before calling this method.
     
     \textbf{WARNING: the documentation must ne rewritten for MNL}

     Noting $i$ the chosen alternative, $V$ the $J$ utilities, $p$ the $\ell$ model parameters, and $\alpha$ the scale parameter, we have
     \[
     P(i,V_1,\ldots,V_J,p_1,\ldots,p_\ell,\alpha,\mu) = \frac{e^{\alpha V_i+\ln G_i(e^{\alpha V_1},\ldots, e^{\alpha V_J},p_1,\ldots,p_\ell,\mu)}}{\sum_j e^{\alpha V_j+\ln G_j(e^{\alpha V_1},\ldots, e^{\alpha V_J},p_1,\ldots,p_\ell,\mu)}}
     \]
     and
     \[
     \begin{array}{rcl}
     \ln P(i,V_1,\ldots,V_J,p_1,\ldots,p_\ell,\alpha,\mu) &=& \alpha V_i\\ &+& \ln G_i(e^{\alpha V_1},\ldots, e^{\alpha V_J},p_1,\ldots,p_\ell,\mu) \\&-& \ln\left( \sum_j e^{\alpha V_j}G_j(e^{\alpha V_1},\ldots, e^{\alpha V_J},p_1,\ldots,p_\ell,\mu)\right)
     \end{array}
     \]
     
     The derivatives with respect to $\beta_k$ are given by
     
     \[
     \frac{\partial}{\partial \beta_k} \ln P = \alpha \frac{\partial V_i}{\partial \beta_k} + \frac{1}{G_i} \sum_{j=1}^J \frac{\partial G_i}{\partial x_j} e^{\alpha V_j} \alpha \frac{\partial V_j}{\partial \beta_k} - \frac{1}{\Delta}\sum_j e^{\alpha V_j}  \left( \alpha \frac{\partial V_j}{\partial \beta_k} G_j + \sum_{n=1}^J \frac{\partial G_j}{\partial x_n} e^{\alpha V_n} \alpha \frac{\partial V_n}{\partial \beta_k} \right) 
     \]
     where
\[
\Delta = \sum_j e^{\alpha V_j} G_j
\]
Note that we do not assume here that the $V_j$ are linear-in-parameters, so that $\partial V_j/\partial \beta_k$ may be not trivial. 

The derivatives with respect to the GEV model parameters $p_k$ are given by
\[
\frac{\partial}{\partial p_k} \ln P = \frac{1}{G_i} \frac{\partial G_i}{\partial p_k} - \frac{1}{\Delta} \sum_j e^{\alpha V_j} \frac{\partial G_j}{\partial p_k}
\]

The derivative with respect to the GEV model parameter $\mu$ is given by
\[
\frac{\partial}{\partial \mu} \ln P = \frac{1}{G_i} \frac{\partial G_i}{\partial \mu} - \frac{1}{\Delta} \sum_j e^{\alpha V_j} \frac{\partial G_j}{\partial \mu}
\]

The derivative with respect to the scale parameter $\alpha$ is given by
\[
\frac{\partial}{\partial \alpha} \ln P = V_i + \frac{1}{G_i} \frac{\partial G_i}{\partial \alpha} - \frac{1}{\Delta} \frac{\partial \Delta}{\partial \alpha}
\]
where
\[
\frac{\partial G_i}{\partial \alpha} = \sum_j V_j e^{\alpha V_j} \frac{\partial G_i}{\partial x_j}
\]
and
\[
 \frac{\partial \Delta}{\partial \alpha} = \sum_j \left(V_j e^{\alpha V_j} G_j + e^{\alpha V_j} \frac{\partial G_j}{\partial \alpha} \right)
\]
@param logOfProba if patTRUE, the log of the probability is returned
@param indivId identifier of the individual
    @param drawNumber number of the requested draw
    @param index index of the chosen alternative
    @param utilities value of the utility of each alternative
    @param beta value of the beta parameters
    @param x value of the characteristics
    @param parameters value of the structural parameters of the GEV function
    @param scale scale parameter
    @param noDerivative patTRUE: no derivative is computed. 
                        patFALSE: derivatives may be computed.
    @param compBetaDerivatives If entry k is patTRUE, the derivative with 
                                respect to $\beta_k$ is computed
    @param compParamDerivatives If entry k is patTRUE, the derivative with 
                                respect to parameter $k$ is computed
    @param compMuDerivative If patTRUE, the derivative with respect to $\mu$ 
                            is computed.
    @param betaDerivatives pointer to the vector where the derivatives with respoect to $\beta$ will be accumulated. The result will be meaningful only if noDerivatives is patFALSE and for entries corresponding to patTRUE in compBetaDerivatives. All other cells are undefined. 
    @param Derivatives pointer to the vector where the derivatives with respect to the GEV parameters will be accumulated. The result will be meaningful only if noDerivatives is patFALSE and for entries corresponding to patTRUE in compParamDerivatives. All other cells are undefined. 
    @param muDerivative pointer to the patReal where the derivative with respect to mu will be accumulated. The result will be meaningful only if noDerivatives is patFALSE and  compMuDerivative is patTRUE.
    @param available describes which alternative are available
    @param mu value of the homogeneity parameter of the GEV function
    @param err ref. of the pointer to the error object. 
    @return pointer to the variable where the result will be stored. 
    Same as outputvariable
    
    The derivatives computed in this routine have been successfully checked
    through a comparison with finite difference computation (Wed Jan 24
    22:06:22 2001)
    
    IMPORTANT: if SNP terms are added, each derivative is multiplied
    by factorForDerivOfSnpTerms before being accumulated

  */
    
  patReal evalProbaPerDraw(patBoolean logOfProba, 
			   patObservationData* individual,
			   unsigned long drawNumber,
			   patVariables* beta,
			   const patVariables* parameters,
			   patReal scale,
			   patBoolean noDerivative ,
			   const vector<patBoolean>& compBetaDerivatives,
			   const vector<patBoolean>& compParamDerivatives,
			   patBoolean compMuDerivative,
			   patBoolean compScaleDerivative,
			   patVariables* betaDerivatives,
			   patVariables* paramDerivatives,
			   patReal* muDerivative,
			   patReal* scaleDerivative,
			   const patReal* mu,
			   patSecondDerivatives* secondDeriv,
			   patBoolean snpTerms,
			   patReal factorForDerivOfSnpTerms,
			   patBoolean* success,
			   patError*& err) ;
  



  /**
   */
  patString getModelName(patError*& err) ;

  /**
     @return largest argument used in all calls to the exp() function
   */
  patReal getMaxExpArgument () ;

  /**
     @return smallest argument used in all calls to the exp() function
   */
  patReal getMinExpArgument () ;

  /**
   */
  unsigned long getUnderflow() ;
  /**
   */
  unsigned long getOverflow() ;

  /**
   */
  virtual patString getInfo() ;

  /**
   */
  virtual void generateCppCodePerDraw(ostream& str,
				      patBoolean logOfProba,
				      patBoolean derivatives, 
				      patBoolean secondDerivatives, 
				      patError*& err) ;
  
private :

  patReal minExpArgument ;
  patReal maxExpArgument ;

  unsigned long overflow ;
  unsigned long underflow ;

  vector<patObservationData::patAttributes>* x ;
  unsigned long index ;
  patReal weight ;
  unsigned long indivId ;
  unsigned long J ;
  patVariables V ;
  patVariables expV ;
  patVariables scaleMuV ;
  patReal maxUtility ;
  unsigned long  i,j,k,kk ;
  patReal sumExp ;
  patVariables probVec ;
  patReal lsumExp ;
  patReal result ;
  patReal tmp ;
  unsigned long altId ;
  unsigned long nBeta ;
  vector<vector<patReal> > DV ;
  vector<patVariables> utilDerivatives ;
  patReal theDeriv ;
  patReal theProba ;
  patVariables term ;
  patReal theHessianEntry ;
};

#endif
