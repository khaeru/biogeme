//-*-c++-*------------------------------------------------------------
//
// File name : patDisplay.cc
// Author :    \URL[Michel Bierlaire]{http://rosowww.epfl.ch/mbi}
// Date :      Mon Jun  2 22:59:15 2003
//
//--------------------------------------------------------------------

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <cassert>
#include <iostream>
#include "patDisplay.h"
#include "patFileNames.h"
#include "patAbsTime.h"
#include "patVersion.h"

patDisplay::patDisplay() : screenImportance(patImportance::patDEBUG), 
			   logImportance(patImportance::patDEBUG),
			   logFile(patFileNames::the()->getLogFile().c_str()) {

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
  messages.erase(messages.begin(),messages.end()) ;
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

  patMessage theMessage ;

  theMessage.theImportance = aType ;
  theMessage.text = text ;
  theMessage.fileName = fileName ;
  theMessage.lineNumber = lineNumber ;

  patAbsTime now ;
  now.setTimeOfDay() ;
  theMessage.theTime = now ;

  if (aType <= screenImportance) {
    if (screenImportance < patImportance::patDEBUG) {
      cout << theMessage.shortVersion() << endl << flush  ;
    }
    else {
      cout << theMessage.fullVersion() << endl << flush  ;
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
  messages.push_back(theMessage) ;
}

void patDisplay::setScreenImportanceLevel(const patImportance& aType) {
  screenImportance = aType ;
}
  
void patDisplay::setLogImportanceLevel(const patImportance& aType) {
  logImportance = aType ;
}

void patDisplay::initProgressReport(const patString message,
			unsigned long upperBound) {
  

}

patBoolean patDisplay::updateProgressReport(unsigned long currentValue) {
  return patTRUE ;
}

void patDisplay::terminateProgressReport() {

}
