//-*-c++-*------------------------------------------------------------
//
// File name : patVersion.cc
// Author:     Michel Bierlaire, EPFL
// Date  :     Mon Jul  9 15:18:09 2001
//
//--------------------------------------------------------------------

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <sstream>
#include <cassert>
#include "patVersion.h"
#include "patCompilationDate.h"
#include "patUser.h"

patVersion::patVersion() :
  version(PACKAGE_VERSION),
  date(__PATCOMPILATIONDATE) ,
  user(__PATUSER) ,
  copyright("Michel Bierlaire, EPFL") {
  stringstream str1 ;
  str1 << PACKAGE_STRING ;
  versionInfo = patString(str1.str()) ;
  stringstream str4 ;
  str4 << date  ;
  versionDate = patString(str4.str()) ;
  versionInfoDate = patString(str1.str() + " [" + str4.str() + "]") ;
  stringstream str2 ;
  str2 << "Michel Bierlaire, EPFL" ;
  versionInfoAuthor = patString(str2.str()) ;
  stringstream str3 ;
  str3 << "-- Compiled by " << user  ;
  versionInfoCompiled = patString(str3.str()) ;
}
  

patVersion* patVersion::the() {
  static patVersion* singleInstance = NULL ;
  if (singleInstance == NULL) {
    singleInstance = new patVersion() ;
    assert(singleInstance != NULL) ;
  }
  return singleInstance ;
}

patString patVersion::getVersion() const {
  return version ;
}


patString patVersion::getDate() const {
  return date ;
}

patString patVersion::getVersionInfoDate() const {
  
  return patString(versionInfoDate) ;
}

patString patVersion::getVersionInfo() const {
  
  return patString(versionInfo) ;
}

patString patVersion::getVersionDate() const {
  
  return patString(versionDate) ;
}

patString patVersion::getVersionInfoAuthor() const {
  return versionInfoAuthor ;
}

patString patVersion::getVersionInfoCompiled() const {
  return versionInfoCompiled ;
}

patString patVersion::getCopyright() const {
  return copyright ;
}

patString patVersion::getUser() const {
  return user ;
}
