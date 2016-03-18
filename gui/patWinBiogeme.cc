//-*-c++-*------------------------------------------------------------
//
// File name : patWinBiogeme.cc
// Author :    \URL[Michel Bierlaire]{http://transp-or2.epfl.ch}
// Date :      Fri Jun 22 17:54:39 2007
//
//--------------------------------------------------------------------

#include <iostream> 

#include <FL/Enumerations.H>
#include <FL/Fl_Help_View.H>
#include "config.h"
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#else
#include <direct.h>
#endif
#include "patWinBiogeme.h"
#include "patVersion.h"
#include "patFileNames.h"
#include "patBiogeme.h"

patWinBiogeme::patWinBiogeme() :   
  Fl_Window(800,600),
  modelFile(120,20,660-120,30,"Model spec. file:"),
  modelFileButton(680,20,780-680,30,"Select file"),
  dataFile(120,60,660-120,30,"Data file:"),
  dataFileButton(680,60,780-680,30,"Select file"),
  workingDir(120,100,660-120,30,"Working directory:"),
  estimateButton(120,140,250,30,"Estimate"),
  simulateButton(410,140,250,30,"Simulate"),
  reportFile(120,560,660-120,30,"Report file:"),
  reportFileButton(680,560,780-680,30,"Display file"),
  directory(".") {

  output = new Fl_Text_Display(20, 190, 760, 550-190) ;
  outputBuffer = new Fl_Text_Buffer() ;
  output->textfont(FL_COURIER) ;
  patDisplay::the()->setUserPointer((void*)this) ;
  modelFileButton.callback(selectModelFileName,(void*)this) ;
  dataFileButton.callback(selectDataFileName,(void*)this) ;
  estimateButton.callback(estimateBiogeme,(void*)this) ;
  simulateButton.callback(simulateBiogeme,(void*)this) ;
  reportFileButton.callback(displayReportFile,(void*)&reportFileName) ;
  reportFileButton.deactivate() ;
  modelFileButton.box(FL_PLASTIC_DOWN_BOX);
  dataFileButton.box(FL_PLASTIC_DOWN_BOX);
  estimateButton.box(FL_PLASTIC_DOWN_BOX);
  simulateButton.box(FL_PLASTIC_DOWN_BOX);
  reportFileButton.box(FL_PLASTIC_DOWN_BOX);

  color((Fl_Color)26) ;
  //  color((Fl_Color)17) ;
  output->buffer(*outputBuffer) ;
  label(patVersion::the()->getVersionInfoDate().c_str()) ;
  resizable(output) ;
  end() ;
  show() ;
}


void patWinBiogeme::addLogMessage(patString m) {
  static int row = 0 ;
  outputBuffer->append(m.c_str()) ;
  outputBuffer->append("\n") ;
  output->scroll(row,0) ;
  ++row ;
  Fl::check() ;

}


void selectModelFileName(Fl_Widget *w, void* userdata) {
  patWinBiogeme* thePtr = (patWinBiogeme*) userdata ;

  // Fl_File_Chooser modelFileChooser(thePtr->directory.c_str(),
  // 				  "*.mod",
  // 				  Fl_File_Chooser::SINGLE,
  // 				  "Model specification file") ;

  Fl_Native_File_Chooser modelFileChooser ;
  modelFileChooser.title("Model specification file") ;
  modelFileChooser.type(Fl_Native_File_Chooser::BROWSE_FILE) ;

  modelFileChooser.filter("*.mod") ;
  
  switch (modelFileChooser.show()) {
  case -1:
    break ;
  case 1:
    break ;
  default:
    patString theFile = thePtr->getFileNameFromFullPath(patString(modelFileChooser.filename()),patFALSE) ;
    patString theDirectory = thePtr->getDirectoryFromFullPath(patString(modelFileChooser.filename())) ;
    thePtr->modelSpecFileName = theFile ;
    thePtr->directory = theDirectory;
#ifdef HAVE_CHDIR
    chdir(theDirectory.c_str()) ;
#endif
    thePtr->workingDir.value(theDirectory.c_str()) ;
    thePtr->modelFile.value(theFile.c_str()) ;
    break;    
  }
}

void selectDataFileName(Fl_Widget *w, void* userdata) {
  patWinBiogeme* thePtr = (patWinBiogeme*) userdata ;
  Fl_Native_File_Chooser dataFileChooser;
  dataFileChooser.title("Data file");
  dataFileChooser.type(Fl_Native_File_Chooser::BROWSE_FILE);
  dataFileChooser.filter("*.dat");
  switch(dataFileChooser.show()) {case -1: break; case 1: break;
  default:
    patString theFile = patString(dataFileChooser.filename()) ;
    thePtr->dataFileName = theFile ;
    thePtr->dataFile.value(dataFileChooser.filename()) ;
    break;
  };
}

void estimateBiogeme(Fl_Widget *w, void* userdata) {
  patWinBiogeme* thePtr = (patWinBiogeme*) userdata ;
  patFileNames::the()->setModelName(thePtr->modelSpecFileName) ;  
  patFileNames::the()->addSamFile(thePtr->dataFileName) ;
  DEBUG_MESSAGE("*** Start estimation ***") ;

  patBiogeme biogeme ;
  patError* err = NULL ;
  DEBUG_MESSAGE("*** Load model and sample ***") ;
  biogeme.loadModelAndSample(err) ;
  if (err != NULL) {
    //    Fl_Color currentCol = thePtr->output->textcolor() ;
    //    thePtr->output->textcolor(FL_RED) ;
    WARNING(err->describe()) ;
//     thePtr->output->textcolor(currentCol) ;
    thePtr->deactivateButtons() ;
    return ;
  }
  biogeme.estimate(NULL,err) ;
  if (err != NULL) {
//     Fl_Color currentCol = thePtr->output->textcolor() ;
//     thePtr->output->textcolor(FL_RED) ;
    thePtr->outputBuffer->append(err->describe().c_str()) ;
//     thePtr->output->textcolor(currentCol) ;
    thePtr->deactivateButtons() ;
    return ;
  }

  thePtr->reportFileName = patFileNames::the()->getRepFile(err) ;
  if (err != NULL) {
    WARNING(err->describe()) ;
    thePtr->deactivateButtons() ;
    return ;
  }

  thePtr->reportFile.value(thePtr->reportFileName.c_str()) ;
  thePtr->reportFileButton.activate() ;
  thePtr->deactivateButtons() ;
  return ;
    
}

void simulateBiogeme(Fl_Widget *w, void* userdata) {
  patWinBiogeme* thePtr = (patWinBiogeme*) userdata ;
  patFileNames::the()->setModelName(thePtr->modelSpecFileName) ;  
  patFileNames::the()->addSamFile(thePtr->dataFileName) ;
  DEBUG_MESSAGE("*** Start estimation ***") ;
  patBiogeme biogeme ;
  patError* err = NULL ;
  DEBUG_MESSAGE("*** Load model and sample ***") ;
  biogeme.loadModelAndSample(err) ;
  
  if (err != NULL) {
    WARNING(err->describe()) ;
    return ;
  }
  
  biogeme.sampleEnumeration(NULL,0,0,err) ;
  if (err != NULL) {
    WARNING(err->describe()) ;
    return ;
  }
  
  thePtr->reportFileName = patFileNames::the()->getEnuFile(err) ;
  if (err != NULL) {
    WARNING(err->describe()) ;
    return ;
  }
  
  thePtr->reportFile.value(thePtr->reportFileName.c_str()) ;
  thePtr->reportFileButton.activate() ;
  thePtr->deactivateButtons() ;
  return ;
  
}

patString patWinBiogeme::getFileNameFromFullPath(patString fullPath,
						 patBoolean keepExtension) {
  if (keepExtension) {
    patString result(fullPath.begin() + fullPath.find_last_of("/") + 1,fullPath.end()) ;
    return result ;
  }
  else {
    patString result(fullPath.begin() + fullPath.find_last_of("/") + 1,
		     fullPath.begin() + fullPath.find_last_of(".")) ;
    return result ;

  }
}

patString patWinBiogeme::getDirectoryFromFullPath(patString fullPath) {
  patString result(fullPath.begin() , 
		   fullPath.begin() + fullPath.find_last_of("/")) ;
  return result ;
}

void displayReportFile(Fl_Widget *w, void* userdata) {
  patString* theFileName = (patString*) userdata ;
  if (theFileName == NULL) {
    return ;
  }

  Fl_Window* theReportWin = new Fl_Window(800,600,theFileName->c_str()) ;
  theReportWin->begin() ;
  Fl_Text_Display* theView = new Fl_Text_Display(0,0,800,600) ;
  Fl_Text_Buffer* fileBuffer = new Fl_Text_Buffer() ;
  theView->buffer(*fileBuffer) ;
  theView->textfont(FL_COURIER) ;
  fileBuffer->insertfile(theFileName->c_str(),0) ;
  theReportWin->resizable(theView) ;
  theReportWin->end() ;
  theReportWin->show() ;
  while(theReportWin->shown()){ 
    Fl::wait(); 
  }
  return ;
}

void patWinBiogeme::deactivateButtons() {
  modelFileButton.deactivate() ;
  dataFileButton.deactivate() ;
  estimateButton.deactivate() ;
  simulateButton.deactivate() ;

}

void patWinBiogeme::setOutputColor(Fl_Color theCol) {
  output->textcolor(theCol) ;
}
