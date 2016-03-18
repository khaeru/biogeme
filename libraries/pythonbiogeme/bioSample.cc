
//-*-c++-*------------------------------------------------------------
//
// File name : bioSample.cc
// Author :    Michel Bierlaire and Mamy Fetiarison
// Date :      Fri Jun 19 10:35:31  2009
//
//--------------------------------------------------------------------


#include <set>
#include <iostream>
#include <fstream>
#include "bioSample.h"
#include "bioParameters.h"
#include "patDisplay.h"
#include "patErrMiscError.h"
#include "patErrOutOfRange.h"
#include "bioIteratorInfoRepository.h"
#include "bioMetaIterator.h"
#include "bioRowIterator.h"
#include "bioIdIterator.h"
#include "bioLiteralRepository.h"
#include "bioLiteralValues.h"
#include "bioExpression.h"

bioSample::bioSample(patString n, patError*& err) : fileName(n) {
  // Read the headers
  // Open the file
  patString line ;
  patString str ;
  ifstream in(fileName.c_str());
  if (!in || in.fail())
  {
    stringstream str ;
    str << "Cannot open the file " << fileName ;
    err = new patErrMiscError(str.str()) ;
    return ;
  }

  // The first header is the row id

  patString rowIdHeader =  bioParameters::the()->getValueString(patString("HeaderOfRowId"),err) ; 
  if (err != NULL) {
    WARNING(err->describe()) ;
    return ;
  }
  headers.push_back(rowIdHeader) ;
  // The column ID is set to patBadId. An virtualcolumn ID will be
  // assigned by the literal repository.
  bioLiteralRepository::the()->addVariable(rowIdHeader,patBadId,err) ;


  // Read headers
  patULong count = 1 ;
  getline(in, line) ;
  istringstream iss(line);
  while (iss >> str) {
    this->headers.push_back(str) ;
    bioLiteralRepository::the()->addVariable(str,count,err) ;
    if (err != NULL) {
      WARNING(err->describe()) ;
      return ;
    }
    ++count ;
  }

  numberOfRealColumns = count ;
  //  DETAILED_MESSAGE("There are " << numberOfRealColumns << " columns in the file, including a column containing the row number, with header " << rowIdHeader) ;
  

  // DETAILED_MESSAGE("Headers from the data file:") ;
  // printHeaders() ;
  in.close() ;

  generatePythonHeaders(err) ;
   if (err != NULL) {
     WARNING(err->describe()) ;
     return ;
   }

}

bioSample::~bioSample() {
}

void bioSample::readFile(bioExpression* exclude, 
			 vector<pair<patString, bioExpression*> >* userExpressions, 
			 patError*& err) {

  DEBUG_MESSAGE("***** " << fileName << " *****") ;
  // only to displa ytotal  line number while reading
  ifstream f(fileName.c_str());
  patString l;
  int nbLines=0;
  while(std::getline(f, l)){
    nbLines++;
  }

  patString rowIdHeader =  bioParameters::the()->getValueString(patString("HeaderOfRowId"),err) ; 
  if (err != NULL) {
    WARNING(err->describe()) ;
    return ;
  }
  headers.push_back(rowIdHeader) ;

  for (vector<pair<patString, bioExpression*> >::iterator i = userExpressions->begin() ;
       i != userExpressions->end() ;
       ++i) {
    headers.push_back(i->first) ;
  }

  ifstream in(fileName.c_str());
  if (!in || in.fail()) {
    stringstream str ;
    str << "Cannot open the file " << fileName ;
    err = new patErrMiscError(str.str()) ;
    return ;
  }
  // Skip headers
  patString line ;
  getline(in, line) ;


  GENERAL_MESSAGE("Read sample file: " << fileName) ;

  patReal rowNumber(0.0) ;
  patULong physicalRowNumber(0) ;
  bioLiteralValues::the()->eraseValues() ;

  
  // patULong nbrUserExpressions((userExpressions == NULL)
  // 			      ?0
  // 			      :userExpressions->size()) ;

  vector<patReal> listOfV(headers.size()) ;

  // Read and store values
  //fprintf(stdout, "\rReading\n");

  char chars[] = {'-', '\\', '|', '/'};
  
  nbrOfExcludedRows = 0 ;

  while(in && !in.eof()) {
    ++physicalRowNumber ;

    //fprintf(stdout, "%c [%5lu/%i]", chars[physicalRowNumber % sizeof(chars)], physicalRowNumber, nbLines);
    //fflush(stdout);

    listOfV[0] = rowNumber ;
    in >> listOfV[1] ;
    if (!in.eof()) {
      for (unsigned long i = 2 ;  i < numberOfRealColumns ; ++i) {
	       if (!in.eof()) {
	         in >> listOfV[i] ;
	       }
      }
      // Include user expressions
      if (userExpressions != NULL) {
      	for (vector<pair<patString, bioExpression*> >::iterator i = 
      	       userExpressions->begin() ;
      	     i != userExpressions->end() ;
      	     ++i) {

      	  i->second->setVariables(&listOfV) ;
      	  patReal value = i->second->getValue(patFALSE, patLapForceCompute, err) ;   
      	  patULong col = 
      	    bioLiteralRepository::the()->getColumnIdOfVariable(i->first,err) ;
      	  if (err != NULL) {
      	    WARNING(err->describe()) ;
      	    return ;
      	  }
      	  listOfV[col] = value ;
      	}
      }
      patReal ex ;
      if (exclude != NULL) {
      	exclude->setVariables(&listOfV) ;
      	ex = exclude->getValue(patFALSE, patLapForceCompute, err) ;
      	if (err != NULL) {
      	  WARNING(err->describe()) ;
      	  return ;
      	}
      }
      else {
	     ex = 0 ;
      }

      if (ex == 0) {
      	rowNumber += 1.0 ;
      	samples.push_back(listOfV) ;
      }
      else {
	++nbrOfExcludedRows ;
      }
    }
    //fprintf(stdout, "\r");
  }

  // Close file
  in.close() ;


}


patULong bioSample::size() const {
  return samples.size() ;
}


vector<patString> bioSample::getHeaders() const {
  return headers ;
}


patVariables* bioSample::at(patULong index) {
  if (index >= samples.size()) {
    return NULL ;
  }
  return &(samples[index]) ;
}


patULong bioSample::getIndexFromHeader(patString h,patError*& err) const {
  patULong result =  bioLiteralRepository::the()->getColumnIdOfVariable(h,err) ;
  if (err != NULL) {
    WARNING(err->describe()) ;
    return patBadId ;
  }
  return result ;
}


patString bioSample::getFileName() {
  return fileName ;
}


patBoolean bioSample::isEqual(bioSample s) {
  return s.getFileName() == this->fileName ;
}


void bioSample::printSample() {
  DEBUG_MESSAGE("*****************") ;
  printHeaders() ;
  patULong cpt = 0;
  for (vector<patVariables>::const_iterator it=samples.begin() ;
       it != samples.end() ;
       it++ ) {
    std::cout << "[" << cpt++ << "] " << *it << endl ;
  }
  std::cout << std::endl ;
}


void bioSample::addObservationData(patVariables data) {
  samples.push_back(data) ;
}


void bioSample::printHeaders() {
  for (vector<patString>::const_iterator it=headers.begin() ;
       it != headers.end();
       it++ ) {
     std::cout << *it << "    ";
  }
  std::cout << std::endl ;
}


void bioSample::computeMapOfDatabase(patError*& err) {

  vector<patString> listOfIterators = 
    bioIteratorInfoRepository::the()->getListOfIterators() ;

  map<patULong, patReal> currentValue ;
  
  for (vector<patString>::iterator i = listOfIterators.begin() ;
       i != listOfIterators.end() ;
       ++i) {
    bioIteratorType theType = bioIteratorInfoRepository::the()->getType(*i,err) ;  if (err != NULL) {
      WARNING(err->describe()) ;
      return ;
    }

    if (theType == ROW || theType == META) {
      patString indexName =  bioIteratorInfoRepository::the()->getIndexName(*i,err) ;
      if (err != NULL) {
	WARNING(err->describe()) ;
	return ;
      }
      patULong colId = getIndexFromHeader(indexName,err) ;
      if (err != NULL) {
	WARNING(err->describe()) ;
	return ;
      }
      if (colId == patBadId) {
	stringstream str ;
	str << "Column " << indexName << " characterizing iterator " << *i << " is unknown in file " << fileName ;
	err = new patErrMiscError(str.str()) ;
	WARNING(err->describe()) ;
	return ;
      }
      currentValue[colId] = patReal() ;
    }
  }

  for (patULong row = 0 ;
       row < samples.size() ;
       ++row) {
    if (row == 0) {
      for (map<patULong, patReal>::iterator col = currentValue.begin() ;
	   col != currentValue.end() ;
	   ++col) {
	bioIteratorInfoRepository::the()->addNewRowId(headers[col->first],0) ;
        col->second = samples[0][col->first] ;
      }
    }
    else {
      for (map<patULong, patReal>::iterator col = currentValue.begin() ;
	   col != currentValue.end() ;
	   ++col) {
	if (col->second != samples[row][col->first]) {
	  bioIteratorInfoRepository::the()->addNewRowId(headers[col->first],row) ;
	  col->second = samples[row][col->first] ;
	}
      }
    }
  }

  bioIteratorInfoRepository::the()->computeRowPointers(samples.size(),err) ;
  if (err != NULL) {
    WARNING(err->describe()) ;
    return ;
  }

}

bioMetaIterator* bioSample::createMetaIterator(bioIteratorSpan theSpan, bioIteratorSpan theThreadSpan, patError*& err) const {
  
  //  DEBUG_MESSAGE("CREATE META ITERATOR FOR THREAD " << theThreadSpan) ;

  if (theSpan.lastRow == patBadId) {
    // theSpan.lastRow = size()-1 ;
    theSpan.lastRow = size() ;
  }
  bioMetaIterator* ptr = new bioMetaIterator(&samples,
					     theSpan, 
					     theThreadSpan,
					     err) ;
  if (err != NULL) {
    WARNING(err->describe()) ;
  }
  return ptr ;
}

bioRowIterator* bioSample::createRowIterator(bioIteratorSpan theSpan, bioIteratorSpan theThreadSpan, patBoolean printMessages, patError*& err) const {
  
  if (theSpan.lastRow == patBadId) {
    // theSpan.lastRow = size()-1 ;
    theSpan.lastRow = size() ;
  }
  bioRowIterator* ptr = new bioRowIterator(&samples,theSpan,theThreadSpan,printMessages,err) ;
  if (err != NULL) {
    WARNING(err->describe()) ;
    return NULL ;
  }

  return ptr ;
}


bioRowIterator* bioSample::createIdIterator(bioIteratorSpan theSpan, bioIteratorSpan theThreadSpan, patString header, patError*& err) const {

  if (theSpan.lastRow == patBadId) {
    // theSpan.lastRow = size()-1 ;
    theSpan.lastRow = size() ;
  }

  patULong indexId = getIndexFromHeader(header,err) ;
  if (err != NULL) {
    WARNING(err->describe()) ;
    return NULL ;
  }
  if (indexId == patBadId) {
    stringstream str ;
    str << "Header " << header << " is unknown" ;
    err = new patErrMiscError(str.str()) ;
    WARNING(err->describe()) ;
    return NULL ;
  }

  DEBUG_MESSAGE("Create an iterator on column " << indexId << ": " << header) ;
  bioIdIterator* ptr = new bioIdIterator(&samples,theSpan,theThreadSpan,indexId,err) ;
  if (err != NULL) {
    WARNING(err->describe()) ;
    return NULL ;
  }

  return ptr ;
  
}


void bioSample::prepare(bioExpression* exclude, 
			vector<pair<patString, bioExpression*> >* userExpressions, 
			patError*& err) {
  DEBUG_MESSAGE("Prepare data...") ;

  readFile(exclude,userExpressions,err) ;
  if (err != NULL) {
    WARNING(err->describe()) ;
    return ;
  }
  computeMapOfDatabase(err) ;
  if (err != NULL) {
    WARNING(err->describe()) ;
    return ;
  }
  
}

void bioSample::generatePythonHeaders(patError*& err) {
  ofstream f("headers.py") ;
  f << "from biogeme import *" << endl ;
  for (vector<patString>::iterator i = headers.begin() ;
       i != headers.end() ;
       ++i) {
    f << *i << "=Variable('" << *i << "')" << endl ;
  }
  f.close() ;
  // ofstream g("headersSimulation.py") ;
  // for (unsigned long i = 0 ; i < headers.size() ; ++i) {
  //   g << headers[i] << " = currentRow[" << i << "]" << endl ;
  // }
  // g.close() ;
}

vector<patULongLong> bioSample::getListOfIds(patString header, patError*& err) {
  patULong indexId = getIndexFromHeader(header,err) ;
  if (err != NULL) {
    WARNING(err->describe()) ;
    return vector<patULongLong>() ;
  }
  if (indexId == patBadId) {
    stringstream str ;
    str << "Header " << header << " is unknown" ;
    err = new patErrMiscError(str.str()) ;
    WARNING(err->describe()) ;
    return vector<patULongLong>() ;
  }
  set<patULongLong> listOfUniqueIds ;
  for (vector<patVariables>::iterator row = samples.begin();
       row != samples.end() ;
       ++row) {
    listOfUniqueIds.insert(patULongLong((*row)[indexId])) ;
  }
  DETAILED_MESSAGE("Variable " << header << " has " << listOfUniqueIds.size() << " different values in the sample") ;
  vector<patULongLong> result ;
  for (set<patULongLong>::iterator i = listOfUniqueIds.begin() ;
       i != listOfUniqueIds.end() ;
       ++i) {
    result.push_back(*i) ;
  }
  return result ;
}

patULong bioSample::getNumberOfColumns() const {
  return headers.size() ;
}


void bioSample::setColumn(patString column, 
			  vector<patReal> values, 
			  patBoolean duplicate, 
			  patString individual,
			  patError*& err) {

  patULong theColumn = 
    bioLiteralRepository::the()->getColumnIdOfVariable(column,err) ;
  if (err != NULL) {
    WARNING(err->describe()) ;
    return ;
  }

  patULong theId = 
    bioLiteralRepository::the()->getColumnIdOfVariable(individual,err) ;
  if (err != NULL) {
    WARNING(err->describe()) ;
    return ;
  }

  patReal currentIndividual = patReal(patBadId) ;
  patULong currentIndex = 0 ;
  for ( patULong row = 0 ;
       row < samples.size() ;
       ++row) {
    if (samples[row][theId] != currentIndividual) {
      if (currentIndex >= values.size()) {
	err = new patErrOutOfRange<patULong>(currentIndex,0,values.size()-1) ;
	WARNING(err->describe()) ;
	return ;
      }
      samples[row][theColumn] = values[currentIndex] ;
      ++currentIndex ;
      currentIndividual = samples[row][theId] ;
    }
    else {
      if (duplicate) {
	samples[row][theColumn] = values[currentIndex] ;
      }
    }
  }
}

vector<patReal> bioSample::getColumn(patString column, patString individual, patError*& err) {

  vector<patReal> result ;
  patULong theColumn = 
    bioLiteralRepository::the()->getColumnIdOfVariable(column,err) ;
  if (err != NULL) {
    WARNING(err->describe()) ;
    return vector<patReal>();
  }

  patULong theId = 
    bioLiteralRepository::the()->getColumnIdOfVariable(individual,err) ;
  if (err != NULL) {
    WARNING(err->describe()) ;
    return vector<patReal>();
  }

  patReal currentIndividual = patReal(patBadId) ;
  //  patULong currentIndex = 0 ;
  for ( patULong row = 0 ;
       row < samples.size() ;
       ++row) {
    if (samples[row][theId] != currentIndividual) {
      result.push_back(samples[row][theColumn]) ;
      currentIndividual = samples[row][theId] ;
    }
  }
  return result ;
}


patULong bioSample::getNbrExcludedRow() const {
  return nbrOfExcludedRows ;
}

const patVariables* bioSample::getFirstRow() {
  if (samples.empty()) {
    return NULL ;
  }
  return &(samples[0]) ;
}
