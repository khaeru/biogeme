//-*-c++-*------------------------------------------------------------
//
// File name : patDisplay.cc
// Author :    \URL[Michel Bierlaire]{http://transp-or2.epfl.ch}
// Date :      Sat Jun 23 11:33:06 2007
//
//--------------------------------------------------------------------

#include <cassert>
#include "patDisplay.h"
#include "patFileNames.h"
#include "patAbsTime.h"
#include "patVersion.h"

#include "patWinBiogeme.h" 


patDisplay::patDisplay() : screenImportance(patImportance::patDEBUG), 
			   logImportance(patImportance::patDEBUG),
			   logFile(patFileNames::the()->getLogFile().c_str()),
			   userPointer(NULL) {
  
  patAbsTime now ;
  now.setTimeOfDay() ;
  logFile << "This file has automatically been generated." << endl ;
  logFile << now.getTimeString(patTsfFULL) << endl ;
  logFile << patVersion::the()->getCopyright() << endl ;
  logFile << endl ;
  logFile << patVersion::the()->getVersionInfoDate() << endl ;
  logFile << patVersion::the()->getVersionInfoAuthor() << endl ;
  logFile << endl ;
  logFile << flush ;

}

patDisplay::~patDisplay() {
  logFile.close() ;
}

patDisplay* patDisplay::the() {
  static patDisplay* singleInstance = NULL ;
  if (singleInstance == NULL) {
    singleInstance = new patDisplay ;
    assert(singleInstance != NULL) ;
  }
  return singleInstance ;
}

void patDisplay::addMessage(const patImportance& aType,
			    const patString& text,
			    const patString& fileName,
			    const patString& lineNumber) {

  patWinBiogeme* thePtr = (patWinBiogeme*) userPointer ;
  patMessage theMessage ;
  theMessage.theImportance = aType ;
  theMessage.text = text ;
  theMessage.fileName = fileName ;
  theMessage.lineNumber = lineNumber ;

  char colorCode ;
  switch (aType()) {
  case patImportance::patFATAL :
    colorCode = 'A' ;
    break ;
  case patImportance::patERROR :
    colorCode = 'A' ;
    break ;
  case patImportance::patWARNING :
    colorCode = 'A' ;
    break ;
  case patImportance::patGENERAL :
    colorCode = 'B' ;
    break ;
  case patImportance::patDETAILED:
    colorCode = 'B' ;
    break ;
  case patImportance::patDEBUG:
    colorCode = 'C' ;
    break ;
    
  }
  
  patAbsTime now ;
  now.setTimeOfDay() ;
  theMessage.theTime = now ;
  messages.push_back(theMessage) ;
  if (thePtr != NULL) {
    if (aType <= screenImportance) {
      if (screenImportance < patImportance::patDEBUG) {
	thePtr->addLogMessage(theMessage.shortVersion()) ;
      }
      else {
	thePtr->addLogMessage(theMessage.fullVersion()) ;
      }
    }
  }
  if (aType <= logImportance) {
    if (logImportance < patImportance::patDEBUG) {
      logFile << theMessage.shortVersion() << endl << flush  ;
    }
    else {
      logFile << theMessage.fullVersion() << endl << flush  ;
    }
  }
}

void patDisplay::setScreenImportanceLevel(const patImportance& aType) {
  screenImportance = aType ;
}
  
void patDisplay::setLogImportanceLevel(const patImportance& aType) {
  logImportance = aType ;
}

void patDisplay::initProgressReport(const patString message,
				    unsigned long upperBound) {
  
  // TO BE IMPLEMENTED
}

patBoolean patDisplay::updateProgressReport(unsigned long currentValue) {
  return patTRUE ;
}

void patDisplay::terminateProgressReport() {
  // TO BE IMPLEMENTED
}

void patDisplay::setUserPointer(void* up) {
  userPointer = up ;
}
