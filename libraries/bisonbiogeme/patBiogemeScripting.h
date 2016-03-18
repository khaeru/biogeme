//-*-c++-*------------------------------------------------------------
//
// File name : patBiogemeScripting.h
// Author :    \URL[Michel Bierlaire]{http://roso.epfl.ch/mbi}
// Date :      Tue Dec 13 15:48:16 2005
//
//--------------------------------------------------------------------

#ifndef patBiogemeScripting_h
#define patBiogemeScripting_h

class patFastBiogeme ;

#ifdef GIANLUCA
#include "patGianlucaBiogeme.h"
#else
#include "patBiogeme.h"
#endif
#include "patPythonResults.h"

/**
   @doc This class is designed to be called from a scripting language
   like Python, Perl or Tcl, through a SWIG interface
   @author \URL[Michel Bierlaire]{http://roso.epfl.ch/mbi}, EPFL (Tue Dec 13 15:48:16 2005)
*/

class patBiogemeScripting {
 public:
  patBiogemeScripting() ;

  // Interface called from main()
  void estimate(int argc, char *argv[]) ;
  void simulation(int argc, char *argv[]) ;

  // Interface called from python
  void estimate(std::string modelName, 
		patPythonReal** data, int row, int columns,
		vector<std::string> headers,
		struct patPythonResults* theResults) ;
  void simulate(std::string modelName, 
		patPythonReal** data, int row, int columns,
		vector<std::string> headers,
		patPythonReal** results, int resRow, int resColumns) ;


  void fastBiogeme(patFastBiogeme* aFct = NULL) ;
 private:
  void invokeBiogeme() ;

 private:

  patBoolean fastBio ;

#ifdef GIANLUCA
  xxxpatGianlucaBiogeme biogeme ;
#else
  patBiogeme biogeme ;
#endif
  patPythonResults* pythonRes ;
};

#endif
