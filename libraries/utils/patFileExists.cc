//-*-c++-*------------------------------------------------------------
//
// File name : patFileExists.cc
// Author :    Michel Bierlaire
// Date :      Mon Dec 21 16:27:31 1998
//
//--------------------------------------------------------------------
//

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <fstream>
#include "patDisplay.h"
#include "patFileExists.h"

/// Checks if a file exists using the stat() function
patBoolean patFileExists::operator()(const patString& fileName) {
  ifstream f(fileName.c_str()) ;
  if (!f) {
    return patFALSE ;
  }
  else {
    f.close() ;
    return patTRUE ;
  }
  
    //struct stat buf ;
  //  int staterr = stat (fileName.c_str(),&buf) ;
  //  stat (fileName.c_str(),&buf) ;
  // return  (S_ISREG(buf.st_mode) || S_ISLNK(buf.st_mode)) ;
}
