//-*-c++-*------------------------------------------------------------
//
// File name : bioMain.cc
// Author :    Michel Bierlaire and Mamy Fetiarison
// Date :      Fri Jul 17 16:41:28 2009
//
//--------------------------------------------------------------------

#include <Python.h>
#include <iomanip>
#include "patFileExists.h"
#include "patMath.h"
#include "bioExpressionRepository.h"
#include "bioBayesianResults.h"
#include "bioMinimizationProblem.h"
#include "patErrNullPointer.h"
#include "patTimeInterval.h"
#include "bioMain.h"
#include "bioRandomDraws.h"
#include "bioParameters.h"
#include "bioSample.h"
#include "bioRowIterator.h"
#include "patDisplay.h"
#include "patFileNames.h"
#include "bioVersion.h"
#include "bioModelParser.h"
#include "trBounds.h"
#include "trSimpleBoundsAlgo.h"
#include "patCfsqp.h"
#include "patSolvOpt.h"
#include "patDonlp2.h"
#include "bioRawResults.h"
#include "bioReporting.h"
#include "bioStatistics.h"
#include "bioSimulation.h"
#include "bioConstraints.h"
#include "bioIterationBackup.h"
#include "bioPrecompiledFunction.h"
#include "bioArithPrint.h"
#include "bioArithBayes.h"

#include "patLap.h"

#ifdef __linux__
#include <sched.h>
#endif
#ifdef __APPLE__
#include <sys/types.h>
#include <sys/sysctl.h>
#endif


bioMain::bioMain() : theFunction(NULL),theSample(NULL), theModel(NULL), theStatistics(NULL), theConstraints(NULL), theRawResults(NULL), theReport(NULL), theRepository(NULL) {

}

bioMain::~bioMain() {
  DELETE_PTR(theFunction) ;
  DELETE_PTR(theSample) ;
  DELETE_PTR(theModel) ;
  DELETE_PTR(theStatistics) ;
  DELETE_PTR(theConstraints) ;
  DELETE_PTR(theRawResults) ;
  DELETE_PTR(theReport) ;
  DELETE_PTR(theRepository) 
}

void bioMain::setFunction(bioPythonWrapper* f) {
  theFunction = f ;
}

void bioMain::setStatistics(bioStatistics* s) {
  theStatistics = s ;
}

void bioMain::setConstraints(bioConstraints* c) {
  theConstraints= c ;
}

void bioMain::run(int argc, char *argv[]) {

  patError* err(NULL) ;


  if (patFileExists()(bioParameters::the()->getValueString("stopFileName"))) {  
    FATAL("Remove the file " << bioParameters::the()->getValueString("stopFileName")) ;
  }

  theReport = new bioReporting(err) ;
  if (err != NULL) {
    WARNING(err->describe()) ;
    exit(-1) ;
  }
  // Reads the sample file

  patString sampleFile("sample.dat") ;
  if (argc >= 3) {
    sampleFile = patString(patString(argv[2])) ;
  }
  //  DEBUG_MESSAGE("Sample file: " << sampleFile) ;
  
  theSample = new bioSample(sampleFile,err) ;
  if (err != NULL) {
    WARNING(err->describe()) ;
    exit(-1) ;
  }

  //  DEBUG_MESSAGE("Number of columns: " << theSample->getNumberOfColumns()) ;

  // Init Python
#ifdef STANDALONE_EXECUTABLE
  Py_NoSiteFlag = 1;
#endif 
  Py_Initialize() ;

  PyRun_SimpleString("import sys") ;
  PyRun_SimpleString("sys.path.append('./module/')") ;



  // Get model name


  patString fn(argv[1]) ;

  //DEBUG_MESSAGE("Model file: " << argv[1]) ;
  patString pythonFileExtension(".py") ;
  int extensionPos = fn.size()-pythonFileExtension.size() ;

  //Initializing lap couting
  patLap::next();

  // Remove .py extension if present
  if (fn.find(pythonFileExtension, extensionPos)!=string::npos) {
    fn.erase(extensionPos) ;
  }

  patFileNames::the()->setModelName(fn) ;

  //  GENERAL_MESSAGE("READ MODEL FROM " << fn) ;
  bioModelParser* parser = new bioModelParser(fn) ;
  parser->setSampleFile(sampleFile) ;
  theModel = parser->readModel(err) ;
  if (err != NULL) {
    WARNING(err->describe()) ;
    exit(-1) ;
  }

  int displayLevel = bioParameters::the()->getValueInt("biogemeDisplay",err) ; 
  if (displayLevel == 0) {
    patDisplay::the()->setScreenImportanceLevel(patImportance::patGENERAL) ;
  }
  else if (displayLevel == 1) {
    patDisplay::the()->setScreenImportanceLevel(patImportance::patDETAILED) ;
  }
  else {
    patDisplay::the()->setScreenImportanceLevel(patImportance::patDEBUG) ;

  }

  displayLevel = bioParameters::the()->getValueInt("biogemeLogDisplay",err) ; 
  if (displayLevel == 0) {
    patDisplay::the()->setLogImportanceLevel(patImportance::patGENERAL) ;
  }
  else if (displayLevel == 1) {
    patDisplay::the()->setLogImportanceLevel(patImportance::patDETAILED) ;
  }
  else {
    patDisplay::the()->setLogImportanceLevel(patImportance::patDEBUG) ;
  }
  
  theRepository = theModel->getRepository() ;
  if (theRepository == NULL) {
    err = new patErrNullPointer("bioExpressionRepository") ;
    WARNING(err->describe()) ;
    exit(-1) ;
  }
  //  GENERAL_MESSAGE("MODEL READ") ;


  bioExpression* excludeExpression = theRepository->getExpression(theModel->getExcludeExpr()) ;
  vector<pair<patString, patULong> >* userExpressions = theModel->getUserExpressions() ;
  vector<pair<patString,bioExpression*> > ptrUserExpressions ;
  for (vector<pair<patString, patULong> >::iterator i = userExpressions->begin() ;
       i != userExpressions->end() ;
       ++i) {
    bioExpression* aExpr = theRepository->getExpression(i->second) ;
    ptrUserExpressions.push_back(pair<patString,bioExpression*>(i->first,aExpr)) ;
  }
  DEBUG_MESSAGE("Prepare sample") ;
  theSample->prepare(excludeExpression,&ptrUserExpressions,err) ;
  if (err != NULL) {
    WARNING(err->describe()); 
    exit(-1) ;
  }

  // DEBUG_MESSAGE("PRINT SAMPLE") ;
  // theSample->printSample() ;
  // exit(-1) ;

  setParameters(err) ;
  if (err != NULL) {
    WARNING(err->describe()); 
    exit(-1) ;
  }

  bioLiteralRepository::the()->prepareMemory();


  if (theReport == NULL) {
    err = new patErrNullPointer("bioReporting") ;
    WARNING(err->describe());
    exit(-1) ;
  }
  theReport->computeFromSample(theSample,err) ;

  DEBUG_MESSAGE("Sample of size: " << theSample->size()) ;


  ofstream po("__parametersUsed.py") ;
  po << bioParameters::the()->printPythonCode() << endl;
  po.close() ;

  bioRandomDraws::the()->populateListOfIds(theSample,err) ;
  if (err != NULL) {
    WARNING(err->describe()); 
    exit(-1) ;
  }

  bioRandomDraws::the()->generateDraws(theSample,err) ;
  if (err != NULL) {
    WARNING(err->describe()); 
    exit(-1) ;
  }

  //  DEBUG_MESSAGE(*bioRandomDraws::the()) ;

  map<patString, patULong>* stats = theModel->getStatistics() ;
  map<patString, patULong>* formulas = theModel->getFormulas() ;

  bioIteratorSpan __theThreadSpan ;


  if (theModel->mustEstimate()) {

    
    if (theModel->involvesMonteCarlo()) {
      theReport->involvesMonteCarlo() ;
    }
    bioExpressionRepository* theRepository = theModel->getRepository();
    theRepository->setReport(theReport) ;
    patULong exprId = theModel->getFormula() ;
    theFunction = new bioPrecompiledFunction(theRepository,exprId,err) ;
    if (err != NULL) {
      WARNING(err->describe()); 
      exit(-1) ;
    }
    estimate(err) ;
    if (err != NULL) {
      WARNING(err->describe()); 
      exit(-1) ;
    }
  }
  else if (theModel->mustBayesEstimate()) {
    bayesian(err) ;
    if (err != NULL) {
      WARNING(err->describe()); 
      exit(-1) ;
    }
    
  }

  if (stats != NULL) {
    for (map<patString, patULong>::iterator i = stats->begin() ;
	 i != stats->end() ;
	 ++i) {
      bioExpression* expr = theRepository->getExpression(i->second) ;
      if (expr == NULL) {
	err = new patErrNullPointer("bioExpression") ;
	WARNING(err->describe()) ;
	return ;
      }
      expr->setThreadSpan(__theThreadSpan) ;
      expr->setSample(theSample) ;
      expr->setReport(theReport) ;
      patReal result = expr->getValue(patFALSE, patLapForceCompute, err) ;
      if (err != NULL) {
	WARNING(err->describe()); 
	return ;
      }
      theReport->addStatistic(i->first,result) ;
    }

  }


  if (formulas != NULL) {
    for (map<patString, patULong>::iterator i = formulas->begin() ;
	 i != formulas->end() ;
	 ++i) {
      bioExpression* expr = theRepository->getExpression(i->second) ;
      if (expr == NULL) {
	err = new patErrNullPointer("bioExpression") ;
	WARNING(err->describe()) ;
	return ;
      }
      patString theFormula = expr->getExpression(err) ;
      if (expr == NULL) {
	err = new patErrNullPointer("bioExpression") ;
	WARNING(err->describe()) ;
	return ;
      }
      theReport->addFormula(i->first,theFormula) ;
    }
  }






  if (theModel->mustSimulate()) {

    //For simulation, only one thread is used
    int nt = bioParameters::the()->getValueInt("numberOfThreads",err) ; 
    if (err != NULL) {
      WARNING(err->describe()); 
      exit(-1) ;
    }
    if (nt > 1) {
      WARNING("Only one thread is used for simulation."); 
      bioParameters::the()->setValueInt("numberOfThreads",1,err) ; 
    }
    simulate(err) ;
    if (err != NULL) {
      WARNING(err->describe()); 
      exit(-1) ;
    }

  }
  
  patString htmlFile = patFileNames::the()->getHtmlFile(err) ;
  if (err != NULL) {
    WARNING(err->describe()) ;
    return ;
  }
  
  theReport->writeHtml(htmlFile,err) ;
  if (err != NULL) {
    WARNING(err->describe()) ;
    return ;
  }

  if (theModel->mustEstimate()) {
    patString latexFile = patFileNames::the()->getLatexFile(err) ;
    if (err != NULL) {
      WARNING(err->describe()) ;
      return ;
    }
    
    theReport->writeLaTeX(latexFile,err) ;
    if (err != NULL) {
      WARNING(err->describe()) ;
      return ;
    }
  }

}

void bioMain::estimate(patError*& err) {

  if (theSample == NULL) {
    err = new patErrMiscError("No sample has been loaded. Biogeme cannot estimate the parameters") ;
    WARNING(err->describe()) ;
    return ;
  }

  theFunction->setSample(theSample,err) ;
  if (err != NULL) {
    WARNING(err->describe()) ;
    return ;
  }


  // Bounds

  patVariables l = bioLiteralRepository::the()->getLowerBounds() ;
  patVariables u = bioLiteralRepository::the()->getUpperBounds() ;
  trBounds theBounds(l,u) ;
  
  patVariables beta = bioLiteralRepository::the()->getBetaValues(patFALSE) ;
  DEBUG_MESSAGE("Number of betas = " << beta.size()) ;
  patBoolean success ;
  patReal initLikelihood(0.0) ;
  patBoolean computeInit = (bioParameters::the()->getValueInt("computeInitLoglikelihood",err) != 0) ;
  if (err != NULL) {
    WARNING(err->describe()) ;
    return ;
  }
  
  patAbsTime beg ;
  patAbsTime endc ;

  
  if (computeInit) {
    DEBUG_MESSAGE("computing init likelihood");
    beg.setTimeOfDay() ;
    initLikelihood = -theFunction->computeFunction(&beta,&success,err) ;
    DEBUG_MESSAGE("computation results obtained");
    endc.setTimeOfDay() ;
    patTimeInterval til0(beg,endc) ;
    GENERAL_MESSAGE("Init. log-likelihood: " << initLikelihood << " ["<<til0.getLength()<<"]") ;
    if (patAbs(initLikelihood) > 1.0e300) {
      patString theErr("There is a numerical problem with the initial log likelihood.\n It typically happens when one observation is associated with a very low probability,\n so that taking the log generates a very high number.\n Modify the starting values of the parameters.\n You may want to use the SIMULATE feature of pythonbiogeme to identify the cause of the problem.") ;
      err = new patErrMiscError(theErr) ;
      WARNING(err->describe()) ;
      return ;

    }
  }
  patBoolean checkDeriv = (bioParameters::the()->getValueInt("checkDerivatives",err) != 0) ;
  if (err != NULL) {
    WARNING(err->describe()) ;
    return ;
  }

  patBoolean computeHessian =  (bioParameters::the()->getValueInt("deriveAnalyticalHessian",err)) != 0 ;

  // WARNING("Derivatives are always checked") ;
  // checkDeriv = patTRUE ;

  if (checkDeriv) {

    DEBUG_MESSAGE("computing derivative");

    patVariables grad(beta.size(),0.0) ;
    trHessian hess(bioParameters::the()->getTrParameters(err),beta.size()) ;
    beg.setTimeOfDay() ;

    patReal theFct ;
    if (computeHessian) {
      theFct = theFunction->computeFunctionAndDerivatives(&beta,
							  &grad,
							  &hess,
							  &success,
							  err) ;
    }
    else {
      theFct = theFunction->computeFunctionAndDerivatives(&beta,
							  &grad,
							  NULL,
							  &success,
							  err) ;

    }
    endc.setTimeOfDay() ;
    patTimeInterval tig(beg,endc) ;

    patVariables finDiffGrad(beta.size(),0.0) ;
    trHessian finDiffHess(bioParameters::the()->getTrParameters(err),beta.size()) ;
    beg.setTimeOfDay() ;
    theFunction->computeFinDiffGradient(&beta,&finDiffGrad,&success,err) ;
    endc.setTimeOfDay() ;
    patTimeInterval tigfd(beg,endc) ;

    GENERAL_MESSAGE("Value: " << theFct) ;
    GENERAL_MESSAGE("First derivatives") ;
    GENERAL_MESSAGE("+++++++++++++++++") ;
    GENERAL_MESSAGE("Analytical\t\tFin. diff.\tError\tName") ;
    for (patULong i = 0 ; i < grad.size() ; ++i) {
      patString name = bioLiteralRepository::the()->getBetaName(i,patFALSE,err) ;
      if (err != NULL) {
	WARNING(err->describe()) ;
	return ;
      }
      GENERAL_MESSAGE( setprecision(7) << setiosflags(ios::scientific|ios::showpos) << grad[i] << '\t' << finDiffGrad[i] << '\t' << setprecision(2) << 100*(grad[i]-finDiffGrad[i])/grad[i]<<"% \t" << name) ;
    }

    if (computeHessian) {
      theFunction->computeFinDiffHessian(&beta,&finDiffHess,&success,err) ;
      if (err != NULL) {
	WARNING(err->describe()) ;
	return ;
      }
      GENERAL_MESSAGE("Second derivatives") ;
      GENERAL_MESSAGE("++++++++++++++++++") ;
      GENERAL_MESSAGE("Analytical\t\tFin. diff.\tError\tName i\tName j") ;

      for (patULong i = 0 ; i < grad.size() ; ++i) {
	patString namei = bioLiteralRepository::the()->getBetaName(i,patFALSE,err) ;
	for (patULong j = i ; j < grad.size() ; ++j) {
	  patString namej = bioLiteralRepository::the()->getBetaName(j,patFALSE,err) ;
	  if (err != NULL) {
	    WARNING(err->describe()) ;
	    return ;
	  }
	  patReal an = hess.getElement(i,j,err) ;
	  if (err != NULL) {
	    WARNING(err->describe()) ;
	    return ;
	  }
	  patReal fd = finDiffHess.getElement(i,j,err) ;
	  if (err != NULL) {
	    WARNING(err->describe()) ;
	    return ;
	  }
	  GENERAL_MESSAGE( setprecision(7) << setiosflags(ios::scientific|ios::showpos) << an << '\t' << fd << '\t' << setprecision(2) << 100*(an-fd)/an<<"% \t" << namei << '\t' << namej) ;
	}
      }
    }
  }


   // WARNING("Exit for debug") ;
   // exit(0) ;
  //Optimization problem

  bioMinimizationProblem theProblem(theFunction, theBounds,theTrParameters) ;
  

  // Constraints
  if (theConstraints != NULL) {
    theConstraints->addConstraints() ;
    vector<bioConstraintWrapper*>* c = theConstraints->getConstraints() ;
    for (vector<bioConstraintWrapper*>::iterator i = c->begin() ;
	 i != c->end() ;
	 ++i) {
      theProblem.addEqualityConstraint(*i) ;
    }
    
  }
  
  // Optimization algorithm

  beg.setTimeOfDay() ;

  patString selectedAlgo = bioParameters::the()->getValueString("optimizationAlgorithm",err) ;
  if (err != NULL) {
    WARNING(err->describe()) ;
    return ;
  }

  if (theProblem.nNonTrivialConstraints() > 0 && selectedAlgo == "BIO") {
    if (patCfsqp::isAvailable) {
      selectedAlgo = "CFSQP" ;
      GENERAL_MESSAGE("The problem has " << theProblem.nNonTrivialConstraints() << " constraints. Algorithm BIO cannot be used. CFSQP is used instead.") ; 
    }
    else {
      selectedAlgo = "DONLP2" ;
      GENERAL_MESSAGE("The problem has " << theProblem.nNonTrivialConstraints() << " constraints. Algorithm BIO cannot be used. DONLP2 is used instead.") ; 
    }

    patULong nle = theProblem.nNonLinearEq() ;
    for (patULong i = 0 ; i < nle ; ++i) {
      trFunction* f = theProblem.getNonLinEquality(i,err) ;
      if (err != NULL) {
	WARNING(err->describe()) ;
	return ;
      }
      patBoolean success ;
      f->computeFunction(&beta,&success,err) ;
      if (err != NULL) {
	WARNING(err->describe()) ;
	return ;
      }
    }
  }


  if (!theProblem.isFeasible(beta,err)) {
    WARNING("Starting point not feasible") ;
  }

  trNonLinearAlgo* algo ;
  if (selectedAlgo == "BIO") {
    patIterationBackup* aBackup = new bioIterationBackup() ;
    algo = new trSimpleBoundsAlgo(&theProblem,
				  beta,
				  theTrParameters,
				  aBackup,
				  err) ;
    
    if (err != NULL) {
      WARNING(err->describe()) ;
      return ;
    }
  }

  else if (selectedAlgo == "CFSQP") {
      
    DETAILED_MESSAGE("Use CFSQP") ;
      
    patIterationBackup* aBackup = new bioIterationBackup() ;
    patCfsqp* algoCfsqp = new patCfsqp(aBackup,&theProblem) ;

    int mode = bioParameters::the()->getValueInt("cfsqpMode",err) ;
    if (err != NULL) {
      WARNING(err->describe()) ;
      return ;
    }
    int iprint = bioParameters::the()->getValueInt("cfsqpIprint",err) ;
    if (err != NULL) {
      WARNING(err->describe()) ;
      return ;
    }
    int miter = bioParameters::the()->getValueInt("cfsqpMaxIter",err) ;
    if (err != NULL) {
      WARNING(err->describe()) ;
      return ;
    }
    patReal eps = bioParameters::the()->getValueReal("cfsqpEps",err) ;
    if (err != NULL) {
      WARNING(err->describe()) ;
      return ;
    }
    patReal epseqn = bioParameters::the()->getValueReal("cfsqpEpsEqn",err) ;
    if (err != NULL) {
      WARNING(err->describe()) ;
      return ;
    }
    patReal udelta = bioParameters::the()->getValueReal("cfsqpUdelta",err) ;
    if (err != NULL) {
      WARNING(err->describe()) ;
      return ;
    }
    patString stopFile = bioParameters::the()->getValueString("stopFileName",err) ;
    if (err != NULL) {
      WARNING(err->describe()) ;
      return ;
    }
    
    algoCfsqp->setParameters(mode,
			     iprint,
			     miter,
			     eps,
			     epseqn,
			     udelta,
			     stopFile) ;
    algo = algoCfsqp ;
    if (algo == NULL) {
      err = new patErrNullPointer("patCfsqp") ;
      return ;
    }
    algo->defineStartingPoint(beta) ;
  }
  else if (selectedAlgo == "SOLVOPT") {
    patIterationBackup* aBackup = new bioIterationBackup() ;
    algo = new patSolvOpt(theSolvoptParameters,&theProblem) ;
    if (algo == NULL) {
      err = new patErrNullPointer("patSolvOpt") ;
      WARNING(err->describe()) ;
      return ;
    }
    algo->setBackup(aBackup) ;
    algo->defineStartingPoint(beta) ;
  }
  else if (selectedAlgo == "DONLP2") {
    patIterationBackup* aBackup = new bioIterationBackup() ;
    algo = new patDonlp2(theDonlp2Parameters,&theProblem) ;
    if (algo == NULL) {
      err = new patErrNullPointer("patDonlp2") ;
      WARNING(err->describe()) ;
      return ;     
    }
    algo->setBackup(aBackup) ;
    algo->defineStartingPoint(beta) ;
  }
  else {
    stringstream str ;
    str << "Unknown algorithm: " << selectedAlgo ;
    err = new patErrMiscError(str.str()) ;
    WARNING(err->describe()) ;
    return ;
  }

  patString algoName(algo->getName()) ;
  
  DEBUG_MESSAGE("About to run algorithm: " << algoName) ;
  
  patString diag = algo->run(err) ;
  if (err != NULL) {
    WARNING(err->describe()) ;
    return ;
  }
  endc.setTimeOfDay() ;
  patTimeInterval ti0(beg,endc) ;
  
  DETAILED_MESSAGE("Run time interval: " << ti0 ) ;
  DETAILED_MESSAGE("Run time: " << ti0.getLength()) ;

  DEBUG_MESSAGE("get solutions");
  patVariables solution = algo->getSolution(err) ;
  if (err != NULL) {
    WARNING(err->describe()) ;
    return ;
  }

  DEBUG_MESSAGE("set beta values");
  bioLiteralRepository::the()->setBetaValues(solution,err) ;
  if (err != NULL) {
    WARNING(err->describe()) ;
    return ;
  }

  DEBUG_MESSAGE("Raw values");
  theRawResults = new bioRawResults(&theProblem,solution,err) ;
  if (err != NULL) {
    WARNING(err->describe()) ;
    return ;
  }

  DEBUG_MESSAGE("compute estimation results");
  theReport->computeEstimationResults(theRawResults,err) ;
  if (err != NULL) {
    WARNING(err->describe()) ;
    return ;
  }


  theReport->setDiagnostic(diag) ;
  theReport->setIterations(algo->nbrIter()) ;
  theReport->setRunTime(ti0.getLength()) ;
  theReport->setInitLL(initLikelihood) ;
  patString pythonEstimated = patFileNames::the()->getPythonEstimated(err) ;
  if (err != NULL) {
    WARNING(err->describe()) ;
    return ;
  }
  theReport->printEstimatedParameters(pythonEstimated,err) ;

  DEBUG_MESSAGE("estimate done");
}

void bioMain::bayesian(patError*& err) {
  bioArithBayes* expr = theModel->getBayesian() ;
  if (theSample == NULL) {
    err = new patErrMiscError("No sample has been loaded. Biogeme cannot estimate the parameters") ;
    WARNING(err->describe()) ;
    return ;
  }

  expr->setSample(theSample) ;
  expr->setReport(theReport) ;

  patAbsTime beg ;
  patAbsTime endc ;
  beg.setTimeOfDay() ;
  bayesianResults = expr->generateDraws(err) ;
  if (err != NULL) {
    WARNING(err->describe()) ;
    return ;
  }
  endc.setTimeOfDay() ;
  patTimeInterval runTime(beg,endc) ;
  GENERAL_MESSAGE("Run time: " << runTime.getLength()) ;
  theReport->setRunTime(runTime.getLength()) ;

  theReport->computeBayesianResults(&bayesianResults,err) ;
  if (err != NULL) {
    WARNING(err->describe()) ;
    return ;
  }

}

void bioMain::simulate(patError*& err) {



  bioArithPrint* expr = theModel->getSimulation() ;
  expr->checkMonteCarlo(patFALSE,err) ;
  if (err != NULL) {
    WARNING(err->describe()) ;
  }

  if (theSample == NULL) {
    err = new patErrMiscError("No sample has been loaded. Biogeme cannot estimate the parameters") ;
    WARNING(err->describe()) ;
    return ;
  }

  expr->setSample(theSample) ;
  expr->setReport(theReport) ;

  
  //  patVariables beta = bioLiteralRepository::the()->getBetaValues(patFALSE) ;
  //patBoolean success ;
  //patReal result = theFunction->getFunction(&beta,&success,err) ;

  patULong weightId = theModel->getWeightExpr() ;
  bioExpression* weight = theRepository->getExpression(weightId) ;
  expr->simulate(theModel->getVarCovarForSensitivity(),
		 theModel->getNamesForSensitivity(),
		 weight,
		 err) ;
  if (err != NULL) {
    WARNING(err->describe()) ;
    return ;
  }
  
  theReport->setSimulationResults(expr->getSimulationResults(),theModel->mustPerformSensitivityAnalysis()) ;

  if (theModel->mustPerformSensitivityAnalysis()) {
    patVariables beta = bioLiteralRepository::the()->getBetaValues(patFALSE) ;
    
  }
}


void bioMain::statistics(patError*& err) {

  if (theSample == NULL) {
    err = new patErrMiscError("No sample has been loaded. Biogeme cannot estimate the parameters") ;
    WARNING(err->describe()) ;
    return ;
  }

  if (theStatistics == NULL) {
    return ;
  }

  theStatistics->bio__generateStatistics(theSample, err) ;
  if (err != NULL) {
    WARNING(err->describe()) ;
    return ;
  }
  vector<pair<patString,patReal> >* theStats = theStatistics->getStatistics() ;
  for (vector<pair<patString,patReal> >::iterator i = theStats->begin() ;
       i != theStats->end() ;
       ++i) {
    theReport->addStatistic(i->first,i->second) ;
  }
}

void bioMain::setParameters(patError*& err) {
  theTrParameters =bioParameters::the()->getTrParameters(err) ;
  if (err != NULL) {
    WARNING(err->describe()) ;
    return ;
  }

  // Update the number of threads

  patULong nThreads = bioParameters::the()->getValueInt("numberOfThreads",err) ;
  if (err != NULL) {
    WARNING(err->describe()) ;
    return ;
  }
  
  // Actual number of physical cores
  if (nThreads == 0) {
    short numCPU = 4;
#ifdef __APPLE__
    int procnum;
    size_t length = sizeof( procnum );
    if(!sysctlbyname("hw.physicalcpu", &procnum, &length, NULL, 0)) {
      numCPU=procnum;}
    else if (!sysctlbyname("hw.ncpu", &procnum, &length, NULL, 0)) {
      numCPU=procnum;}
#elif defined(__linux__)
    cpu_set_t procset;
    if(!sched_getaffinity(0, sizeof(procset), &procset)) {
      numCPU = CPU_COUNT(&procset);}
#elif defined HAVE_UNISTD_H 
    numCPU = sysconf( _SC_NPROCESSORS_CONF ) ;
#endif
    long int shareProc = bioParameters::the()->getValueInt("shareOfProcessors",err) ;
    if (numCPU <= 0) {
      nThreads = 4;}
    else if (shareProc > 0 ) {
      nThreads = ceil(patReal(numCPU) * patReal(shareProc) / 100.0)  ;
    }
    else {
      nThreads = numCPU ;}
    bioParameters::the()->setValueInt("numberOfThreads",nThreads,err) ;
    if (err != NULL) {
      WARNING(err->describe()) ;
      return ;
    }
    GENERAL_MESSAGE("Nbr of cores reported by the system: " << numCPU) ;
  }

  GENERAL_MESSAGE("Nbr of cores used by biogeme: "        << nThreads) ;

  theSolvoptParameters.errorArgument = bioParameters::the()->getValueReal("solvoptErrorArgument",err) ;
 if (err != NULL) {
   WARNING(err->describe()) ;
   return ;
 }
 theSolvoptParameters.errorFunction = bioParameters::the()->getValueReal("solvoptErrorFunction",err) ;
 if (err != NULL) {
   WARNING(err->describe()) ;
   return ;
 }
 theSolvoptParameters.maxIter = bioParameters::the()->getValueInt("solvoptMaxIter",err) ;
 if (err != NULL) {
   WARNING(err->describe()) ;
   return ;
 }
 theSolvoptParameters.display = bioParameters::the()->getValueInt("solvoptDisplay",err) ;
 if (err != NULL) {
   WARNING(err->describe()) ;
   return ;
 }

 theDonlp2Parameters.epsx = bioParameters::the()->getValueReal("donlp2epsx",err) ;
 if (err != NULL) {
   WARNING(err->describe()) ;
   return ;
 }
 theDonlp2Parameters.delmin = bioParameters::the()->getValueReal("donlp2delmin",err) ;
 if (err != NULL) {
   WARNING(err->describe()) ;
   return ;
 }
 theDonlp2Parameters.smallw = bioParameters::the()->getValueReal("donlp2smallw",err) ;
 if (err != NULL) {
   WARNING(err->describe()) ;
   return ;
 }
 theDonlp2Parameters.epsdif = bioParameters::the()->getValueReal("donlp2epsdif",err) ;
 if (err != NULL) {
   WARNING(err->describe()) ;
   return ;
 }
 theDonlp2Parameters.nreset = bioParameters::the()->getValueInt("donlp2nreset",err) ;
 if (err != NULL) {
   WARNING(err->describe()) ;
   return ;
 }
 theDonlp2Parameters.stopFileName = bioParameters::the()->getValueString("stopFileName",err) ;
 if (err != NULL) {
   WARNING(err->describe()) ;
   return ;
 }

}
