//-*-c++-*------------------------------------------------------------
//
// File name : testing.cc
// Author :    Michel Bierlaire
// Date :      Sun Jul 12 16:54:24 2015
//
//--------------------------------------------------------------------

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "patDisplay.h"
#include "patStatistics.h"
int main(int argc, char *argv[]) {

  patStatistics theStat ;
  for (patULong i=0; i < 100 ; ++i) {
    vector<patReal> d ;
    d.push_back(patReal(i)) ;
    d.push_back(patReal(i/2)) ;
    //    theStat.addData(d) ;
  }

}
