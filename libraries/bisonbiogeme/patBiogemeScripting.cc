//-*-c++-*------------------------------------------------------------
//
// File name : patBiogemeScripting.cc
// Author :    \URL[Michel Bierlaire]{http://roso.epfl.ch/mbi}
// Date :      Tue Dec 13 15:57:07 2005
//
//--------------------------------------------------------------------

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <stdlib.h>
#include "patBiogemeScripting.h"
#include "patFileNames.h"
#include "patVersion.h"
#include "patCfsqp.h"
#include "patDisplay.h"
#include "patParameters.h"
#include "patTimeInterval.h"

patBiogemeScripting::patBiogemeScripting() : fastBio(patFALSE), pythonRes(NULL) {
  cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl ;
  cout << patVersion::the()->getVersionInfoDate() << endl ;
  cout << patVersion::the()->getVersionInfoAuthor() << endl ;
  cout << patVersion::the()->getVersionInfoCompiled() << endl ;
  cout << "See http://biogeme.epfl.ch" << endl  ;
  if (patCfsqp::isAvailable()) {
    cout << "                    !! CFSQP is available !!" << endl ;
  }
  else {
    cout << "                    !! CFSQP is NOT available !!" << endl ;
  }
  cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl  ;
  cout << "\t\"In every non-trivial program there is at least one bug.\"" << endl  ;
  cout << "" << endl  ;

}

void patBiogemeScripting::estimate(int argc, char *argv[]) {
  //%%%%%%%%%%%%%%%%%%
  // 1 Test the arguments
  //%%%%%%%%%%%%%%%%%%

  if ((argc == 2) && patString(argv[1]) == patString("-h")) {
    cout << "Usage: " << argv[0] << " model_name sampleFile1 sampleFile2 sampleFile3 ... " << endl ;
    exit(0 );
  }

  //%%%%%%%%%%%%%%%%%%
  // 2 Define the model name
  //%%%%%%%%%%%%%%%%%%

  
  if (argc > 1) {
    patString modelName(argv[1]) ;
    patFileNames::the()->setModelName(modelName) ;
  }

  //%%%%%%%%%%%%%%%%%%
  // 3 Define the name of the sample files
  //%%%%%%%%%%%%%%%%%%

  if (argc > 2) {
    for (unsigned short i=2 ; i < argc ; ++i) {
      patString sampleFileName(argv[i]) ;
      patFileNames::the()->addSamFile(sampleFileName) ;
    }
  }  
  
  cout << endl ;

  //%%%%%%%%%%%%%%%%%%
  // 5 Invoke biogeme
  //%%%%%%%%%%%%%%%%%%

  invokeBiogeme() ;
}

void patBiogemeScripting::simulation(int argc, char *argv[]) {


  //%%%%%%%%%%%%%%%%%%
  // S1 Test the arguments
  //%%%%%%%%%%%%%%%%%%


  patError* err = NULL ;
  if ((argc == 2) && patString(argv[1]) == patString("-h")) {
    cout << "Usage: " << argv[0] << " model_name sampleFile1 sampleFile2 sampleFile3 ... " << endl ;
    exit(0 );
  }

  if (argc > 1) {
    patString modelName(argv[1]) ;
    patFileNames::the()->setModelName(modelName) ;
  }
  if (argc > 2) {
    for (unsigned short i=2 ; i < argc ; ++i) {
      patString sampleFileName(argv[i]) ;
      patFileNames::the()->addSamFile(sampleFileName) ;
    }
  }  
  
  cout << endl ;

  //%%%%%%%%%%%%%%%%%%
  // S2 Messages
  //%%%%%%%%%%%%%%%%%%


  //%%%%%%%%%%%%%%%%%%
  // S3 Read parameter file
  //%%%%%%%%%%%%%%%%%%

  biogeme.loadModelAndSample(err) ;
  if (err != NULL) {
    WARNING(err->describe()) ;
    exit(-1) ;
  }

#ifndef GIANLUCA
  biogeme.sampleEnumeration(NULL,0,0,err) ;
  if (err != NULL) {
    WARNING(err->describe()) ;
    exit(-1) ;
  }
#endif

}

void patBiogemeScripting::estimate(std::string modelName, 
				   patPythonReal** data, int row, int columns,
				   vector<std::string> head,
				   patPythonResults* theResults) {

  pythonRes = theResults ;
  GENERAL_MESSAGE("Estimate model " << modelName) ;

  patError* err = NULL ;

  cout << "Estimate model " << modelName << endl ;

  patFileNames::the()->setModelName(patString(modelName)) ;

  //%%%%%%%%%%%%%%%%%%
  // 3 Define the  sample 
  //%%%%%%%%%%%%%%%%%%

  biogeme.externalData(data,row,columns,head) ;


  //%%%%%%%%%%%%%%%%%%
  // 4 Messages
  //%%%%%%%%%%%%%%%%%%




  //%%%%%%%%%%%%%%%%%%
  // 5 Invoke biogeme
  //%%%%%%%%%%%%%%%%%%

  invokeBiogeme() ;

  
}

void patBiogemeScripting::invokeBiogeme() {

  patError* err = NULL ;
  patAbsTime t1 ;
  patAbsTime t2 ;
  patAbsTime t3 ;
  t1.setTimeOfDay() ;
  biogeme.loadModelAndSample(err) ;
  if (err != NULL) {
    WARNING(err->describe()) ;
    exit(-1) ;
  }
  t2.setTimeOfDay() ;
  GENERAL_MESSAGE("Run time for data processing: " << patTimeInterval(t1,t2).getLength()) ;

  biogeme.estimate(pythonRes,err) ;
  if (err != NULL) {
    WARNING(err->describe()) ;
    exit(-1) ;
  }

  t3.setTimeOfDay() ;
  GENERAL_MESSAGE("Run time for estimation:      " << patTimeInterval(t2,t3).getLength()) ;
  GENERAL_MESSAGE("Total run time:               " << patTimeInterval(t1,t3).getLength()) ;
  if (biogeme.typeOfRun == patParallelRun) {
    GENERAL_MESSAGE("---- Run complete with " << patParameters::the()->getgevNumberOfThreads() << " processors ----") ;
  }
}


void patBiogemeScripting::simulate(std::string modelName, 
				   patPythonReal** data, int row, int columns,
				   vector<std::string> head,
				   patPythonReal** results, int resRow, int resColumns) {


  GENERAL_MESSAGE("Simulate model " << modelName) ;

  patError* err = NULL ;

  patFileNames::the()->setModelName(patString(modelName)) ;

  //%%%%%%%%%%%%%%%%%%
  // 3 Define the  sample 
  //%%%%%%%%%%%%%%%%%%

  biogeme.externalData(data,row,columns,head) ;


  //%%%%%%%%%%%%%%%%%%
  // 4 Messages
  //%%%%%%%%%%%%%%%%%%




  //%%%%%%%%%%%%%%%%%%
  // 5 Invoke biogeme
  //%%%%%%%%%%%%%%%%%%

  biogeme.loadModelAndSample(err) ;
  if (err != NULL) {
    WARNING(err->describe()) ;
    exit(-1) ;
  }

#ifndef GIANLUCA
  biogeme.sampleEnumeration(results,resRow,resColumns,err) ;
  if (err != NULL) {
    WARNING(err->describe()) ;
    exit(-1) ;
  }
#endif



}

void patBiogemeScripting::fastBiogeme(patFastBiogeme* aFct) {
  if (aFct == NULL) {
    biogeme.typeOfRun = patGeneratingCode ;
  }
  else {
    biogeme.typeOfRun = patParallelRun ;
    biogeme.theFastFunction = aFct ;
  }
}
