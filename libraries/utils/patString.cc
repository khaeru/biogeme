//-*-c++-*------------------------------------------------------------
//
// File name : patString.cc
// Author :    \URL[Michel Bierlaire]{http://rosowww.epfl.ch/mbi}
// Date :      Wed Nov 21 12:13:05 2001
//
//--------------------------------------------------------------------

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "patDisplay.h"
#include "patString.h"
patString noSpace(const patString& aString) {
  patString result = aString ;
  for (unsigned long i = 0 ; i < result.size() ; ++i) {
    if (result[i] == ' ') {
      result[i] = '_' ;
    }
  }
  return result ;
}

patString* replaceAll(patString* theString, patString chain, patString with) {

  if (theString == NULL) {
    return NULL ;
  }

  patString::size_type ip = theString->find(chain) ;
 
  while (ip != string::npos ) {
    theString->replace(ip,chain.size(),with) ;
    ip = theString->find(chain) ;
  }
  return theString ;
}

patString fillWithBlanks(const patString& theStr, 
			 unsigned long n, short justifyLeft) {
  unsigned long currentLength = theStr.size() ;
  if (currentLength >= n) {
    return theStr ;
  }

  patString result ;
  patString blanks(n-currentLength,' ') ;
  if (justifyLeft) {
    result = theStr + blanks ; 
  }
  else {
    result = blanks + theStr ;
  }
  return result ;
}
