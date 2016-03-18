//-*-c++-*------------------------------------------------------------
//
// File name : donlp2Parameters.h
// Date :      Fri Oct 22 14:12:31 2010
//
//--------------------------------------------------------------------

#ifndef donlp2Parameters_h
#define donlp2Parameters_h

class donlp2Parameters {
 public:
  patReal epsx ; // patParameters::the()->getdonlp2Epsx() ;
  patReal delmin ; // patParameters::the()->getdonlp2Delmin() ;
  patReal smallw ; // patParameters::the()->getdonlp2Smallw() ;
  patReal epsdif ; // patParameters::the()->getdonlp2Epsdif() ;
  patULong nreset ; // patParameters::the()->getdonlp2NReset() ;
  patString stopFileName ;

};

#endif 
