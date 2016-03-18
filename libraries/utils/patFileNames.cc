//-*-c++-*------------------------------------------------------------
//
// File name : patFileNames.cc
// Author :    \URL[Michel Bierlaire]{http://rosowww.epfl.ch/mbi}
// Date :      Thu Aug 15 21:50:46 2002
//
//--------------------------------------------------------------------

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <sstream>
#include <cassert>
#include "patDisplay.h"
#include "patConst.h"
#include "patFileNames.h"
#include "patFileExists.h"
#include "patErrMiscError.h"

patFileNames::~patFileNames() {

}
  
patFileNames::patFileNames() : modelName("default"), 
			       defaultName("default") {
  
}

patFileNames* patFileNames::the() {
  static patFileNames* singleInstance = NULL ;
  if (singleInstance == NULL) {
    singleInstance = new patFileNames ;
    assert(singleInstance != NULL) ;
  }
  return singleInstance ;
}

  
patString patFileNames::getModelName() {
  return modelName ;
}

void patFileNames::setModelName(const patString& name) {
  modelName = name ;
}

patString patFileNames::getInputFile(patString extension) {
  stringstream modstr ;
  modstr << modelName << extension ;
  patString inputFile = modstr.str() ;
  if (!patFileExists()(inputFile) && (modelName != defaultName)) {
    GENERAL_MESSAGE(inputFile << " does not exist") ;
    stringstream defstr ;
    defstr << defaultName << extension ;
    inputFile = defstr.str() ;
    GENERAL_MESSAGE("Trying " << inputFile << " instead") ;
  }
  return inputFile ;
}

patString patFileNames::getOutputFile(patString extension, patError*& err) {
  stringstream str ;
  short count = 0 ;
  str << modelName << extension ;
  patString outputFile = str.str() ;
  while (patFileExists()(outputFile)) {
    if (count == 9999) {
      stringstream errstr ;
      errstr << outputFile << " exists. Too many files for model " 
	     << modelName 
	     << ". Please purge your files." ;
      err = new patErrMiscError(errstr.str()) ;
      WARNING(err->describe()) ;
      return patString() ;
    }
    else {
      ++count ;
    }
    stringstream defstr ;
    defstr << modelName << '~' << count << extension ;
    outputFile = defstr.str() ;
  }
  return outputFile ;
}

patString patFileNames::getOutputFileNoBackup(patString extension) {
  stringstream str ;
  str << modelName << extension ;
  patString outputFile = str.str() ;
  return outputFile ;
}


patString patFileNames::getParFile() {
  static patBoolean first = patTRUE ;
  if (first) {
    parFile = getInputFile(".par") ;
    first = patFALSE ;
  }
  return parFile ;
}
  
patString patFileNames::getModFile() {
  static patBoolean first = patTRUE ;
  if (first) {
    modFile = getInputFile(".mod") ;
    first = patFALSE ;
  }
  return modFile ;
}

void patFileNames::addSamFile(patString name) {
  samFile.push_back(name) ;
}
  
patString patFileNames::getSamFile(unsigned short nbr, patError*& err) {
  if (samFile.empty()) {
    samFile.push_back(patString("sample.dat")) ;
  }
  if (nbr >= samFile.size()) {
    stringstream str ;
    str << "Sample file #" << nbr << " does not exist. Only " <<  samFile.size() << " files, numbered from 0 to " << samFile.size()-1 ;
    err = new patErrMiscError(str.str()) ;
    WARNING(err->describe()) ;
    return patString();
  }
  return samFile[nbr] ;
}
  
patString patFileNames::getStaFile(patError*& err) {
  if (err != NULL) {
    WARNING(err->describe()) ;
    return patString() ;
  }
  static patBoolean first = patTRUE ;
  if (first) {
    staFile = getOutputFile(".sta",err) ;
    if (err != NULL) {
      WARNING(err->describe()) ;
      return patString() ;
    }
    first = patFALSE ;
  }
  return staFile ;
}

patString patFileNames::getRepFile(patError*& err) {
  if (err != NULL) {
    WARNING(err->describe()) ;
    return patString() ;
  }
  static patBoolean first = patTRUE ;
  if (first) {
    repFile = getOutputFile(".rep",err) ;
    if (err != NULL) {
      WARNING(err->describe()) ;
      return patString() ;
    }
    first = patFALSE ;
  }
  return repFile ;
}
  
patString patFileNames::getALogitFile(patError*& err) {
  if (err != NULL) {
    WARNING(err->describe()) ;
    return patString() ;
  }
  static patBoolean first = patTRUE ;
  if (first) {
    f12File = getOutputFile(".F12",err) ;
    if (err != NULL) {
      WARNING(err->describe()) ;
      return patString() ;
    }
    first = patFALSE ;
  }
  return f12File ;
}
  
patString patFileNames::getHtmlFile(patError*& err) {
  if (err != NULL) {
    WARNING(err->describe()) ;
    return patString() ;
  }
  static patBoolean first = patTRUE ;
  if (first) {
    htmlFile = getOutputFile(".html",err) ;
    if (err != NULL) {
      WARNING(err->describe()) ;
      return patString() ;
    }
    first = patFALSE ;
  }
  return htmlFile ;
}
  
patString patFileNames::getPythonSpecFile(patError*& err) {
  if (err != NULL) {
    WARNING(err->describe()) ;
    return patString() ;
  }
  static patBoolean first = patTRUE ;
  if (first) {
    pythonSpecFile = getOutputFile(".py",err) ;
    if (err != NULL) {
      WARNING(err->describe()) ;
      return patString() ;
    }
    first = patFALSE ;
  }
  return pythonSpecFile ;
}
  
patString patFileNames::getLatexFile(patError*& err) {
  if (err != NULL) {
    WARNING(err->describe()) ;
    return patString() ;
  }
  static patBoolean first = patTRUE ;
  if (first) {
    latexFile = getOutputFile(".tex",err) ;
    if (err != NULL) {
      WARNING(err->describe()) ;
      return patString() ;
    }
    first = patFALSE ;
  }
  return latexFile ;
}
  
patString patFileNames::getLogFile() {
  static patBoolean first = patTRUE ;
  if (first) {
    logFile = getOutputFileNoBackup(".log") ;
    first = patFALSE ;
  }
  return logFile ;
}
  
patString patFileNames::getEnuFile(patError*& err) {
  if (err != NULL) {
    WARNING(err->describe()) ;
    return patString() ;
  }
  static patBoolean first = patTRUE ;
  if (first) {
    enuFile = getOutputFile(".enu",err) ;
    if (err != NULL) {
      WARNING(err->describe()) ;
      return patString() ;
    }
    first = patFALSE ;
  }
  return enuFile ;
}
 
patString patFileNames::getZhengFosgerauLatex(patError*& err) {
  if (err != NULL) {
    WARNING(err->describe()) ;
    return patString() ;
  }
  static patBoolean first = patTRUE ;
  if (first) {
    stringstream str ;
    short count = 0 ;
    str << modelName << "_zheng.tex" ;
    zfLatex = str.str() ;
    while (patFileExists()(zfLatex)) {
      if (count == 99) {
	stringstream errstr ;
	errstr << zfLatex << " exists. Too many files for model " 
	       << modelName 
	       << ". Please purge your files." ;
	err = new patErrMiscError(errstr.str()) ;
	WARNING(err->describe()) ;
	return patString() ;
      }
      else {
	++count ;
      }
      stringstream defstr ;
      defstr << modelName << "_zheng" << '~' << count << ".tex" ;
      zfLatex = defstr.str() ;
    }

    first = patFALSE ;
  }
  return zfLatex ;

}

patString patFileNames::getZhengFosgerau(patError*& err) {
  if (err != NULL) {
    WARNING(err->describe()) ;
    return patString() ;
  }
  static patBoolean first = patTRUE ;
  if (first) {
    stringstream str ;
    short count = 0 ;
    str << modelName << "_zheng.enu" ;
    zfDetails = str.str() ;
    while (patFileExists()(zfDetails)) {
      if (count == 99) {
	stringstream errstr ;
	errstr << zfDetails << " exists. Too many files for model " 
	       << modelName 
	       << ". Please purge your files." ;
	err = new patErrMiscError(errstr.str()) ;
	WARNING(err->describe()) ;
	return patString() ;
      }
      else {
	++count ;
      }
      stringstream defstr ;
      defstr << modelName << "_zheng" << '~' << count << ".enu" ;
      zfDetails = defstr.str() ;
    }

    first = patFALSE ;
  }
  return zfDetails ;

}

patString patFileNames::getResFile(patError*& err) {
  if (err != NULL) {
    WARNING(err->describe()) ;
    return patString() ;
  }
  static patBoolean first = patTRUE ;
  if (first) {
    resFile = getOutputFile(".res",err) ;
    if (err != NULL) {
      WARNING(err->describe()) ;
      return patString() ;
    }
    first = patFALSE ;
  }
  return resFile ;
}
  
patString patFileNames::getBckFile(patError*& err) {
  if (err != NULL) {
    WARNING(err->describe()) ;
    return patString() ;
  }
  static patBoolean first = patTRUE ;
  if (first) {
    bckFile = getOutputFile(".bck",err) ;
    if (err != NULL) {
      WARNING(err->describe()) ;
      return patString() ;
    }
    first = patFALSE ;
  }
  return bckFile ;
}
  


patString patFileNames::getModelDebug() {
  return patString("model.debug") ;
}

patString patFileNames::getParamDebug() {
  return patString("parameters.out") ;
}

patString patFileNames::getSpecDebug() {
  return patString("__specFile.debug");
}

patString patFileNames::getGnuplotFile(patError*& err) {
  if (err != NULL) {
    WARNING(err->describe()) ;
    return patString() ;
  }
  static patBoolean first = patTRUE ;
  if (first) {
    gnuplotFile = getOutputFile(".gp",err) ;
    if (err != NULL) {
      WARNING(err->describe()) ;
      return patString() ;
    }
    first = patFALSE ;
  }
  return gnuplotFile ;
}

void patFileNames::setNormalDrawsFile(patString fileName) {
  normalDrawFile = fileName ;
}

void patFileNames::setUnifDrawsFile(patString fileName) {
  unifDrawFile = fileName ;
}

patString patFileNames::getNormalDrawsFile(patError*& err) {
  return normalDrawFile ;
}

patString patFileNames::getUnifDrawsFile(patError*& err) {
  return unifDrawFile ;
}

patString patFileNames::getNormalDrawLogFile(patError*& err) {
  if (err != NULL) {
    WARNING(err->describe()) ;
    return patString() ;
  }
  static patBoolean first = patTRUE ;
  if (first) {
    normalDrawLogFile = getOutputFile(".ndraws",err) ;
    if (err != NULL) {
      WARNING(err->describe()) ;
      return patString() ;
    }
    first = patFALSE ;
  }
  return normalDrawLogFile ;
}

patString patFileNames::getUnifDrawLogFile(patError*& err) {
  if (err != NULL) {
    WARNING(err->describe()) ;
    return patString() ;
  }
  static patBoolean first = patTRUE ;
  if (first) {
    unifDrawLogFile = getOutputFile(".udraws",err) ;
    if (err != NULL) {
      WARNING(err->describe()) ;
      return patString() ;
    }
    first = patFALSE ;
  }
  return unifDrawLogFile ;
}


unsigned short patFileNames::getNbrSampleFiles() {
  if (samFile.empty()) {
    samFile.push_back(patString("sample.dat")) ;
  }
  return samFile.size() ;
}

patString patFileNames::getCcCode(patError*& err) {
  if (err != NULL) {
    WARNING(err->describe()) ;
    return patString() ;
  }
  static patBoolean first = patTRUE ;
  if (first) {
    ccCodeName = getOutputFileNoBackup(".cc") ;
    if (err != NULL) {
      WARNING(err->describe()) ;
      return patString() ;
    }
    first = patFALSE ;
  }
  return ccCodeName ;
}

patString patFileNames::getCcCode_f(patError*& err) {
  if (err != NULL) {
    WARNING(err->describe()) ;
    return patString() ;
  }
  static patBoolean first = patTRUE ;
  if (first) {
    stringstream str ;
    str << modelName << "_f.cc" ;
    ccfCodeName = patString(str.str()) ;
    first = patFALSE ;
  }
  return ccfCodeName ;
}

patString patFileNames::getCcCode_s(patError*& err) {
  if (err != NULL) {
    WARNING(err->describe()) ;
    return patString() ;
  }
  static patBoolean first = patTRUE ;
  if (first) {
    stringstream str ;
    str << modelName << "_s.cc" ;
    ccfCodeName = patString(str.str()) ;
    first = patFALSE ;
  }
  return ccfCodeName ;
}


patString patFileNames::getCcCode_g(patError*& err) {
  if (err != NULL) {
    WARNING(err->describe()) ;
    return patString() ;
  }
  static patBoolean first = patTRUE ;
  if (first) {
    stringstream str ;
    str << modelName << "_g.cc" ;
    ccfCodeName = patString(str.str()) ;
    first = patFALSE ;
  }
  return ccfCodeName ;
}

patString patFileNames::getPythonEstimated(patError*& err) {
  if (err != NULL) {
    WARNING(err->describe()) ;
    return patString() ;
  }
  static patBoolean first = patTRUE ;
  if (first) {
    pythonEstimated = getOutputFile("_param.py",err) ;
    if (err != NULL) {
      WARNING(err->describe()) ;
      return patString() ;
    }
    first = patFALSE ;
  }
  return pythonEstimated ;
}

patString patFileNames::getMonteCarloLog(patError*& err) {
  if (err != NULL) {
    WARNING(err->describe()) ;
    return patString() ;
  }
  static patBoolean first = patTRUE ;
  if (first) {
    mcFile = getOutputFile(".mclg",err) ;
    if (err != NULL) {
      WARNING(err->describe()) ;
      return patString() ;
    }
    first = patFALSE ;
  }
  return mcFile ;
}
