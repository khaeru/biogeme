//-*-c++-*------------------------------------------------------------
//
// File name : patWinBiogeme.h
// Author :    \URL[Michel Bierlaire]{http://transp-or2.epfl.ch}
// Date :      Fri Jun 22 17:43:34 2007
//
//--------------------------------------------------------------------

#ifndef patWinBiogeme_h
#define patWinBiogeme_h

#include <FL/Fl.H>
#include <FL/Fl_Window.H>
//#include <FL/Fl_File_Chooser.H>
#include <FL/Fl_Native_File_Chooser.H>
#include <FL/Fl_Text_Display.H>
#include <FL/Fl_Output.H>
#include <FL/Fl_Button.H>
#include "patString.h"
#include "patConst.h"

/**
   userdata is the pointer to patWinBiogeme
*/
void selectModelFileName(Fl_Widget *w, void* userdata) ;
/**
   userdata is the pointer to patWinBiogeme
*/
void selectDataFileName(Fl_Widget *w, void* userdata) ;
/**
   userdata is the pointer to patWinBiogeme
*/
void estimateBiogeme(Fl_Widget *w, void* userdata) ;

/**
   userdata is the pointer to patWinBiogeme
*/
void simulateBiogeme(Fl_Widget *w, void* userdata) ;

/**
 */
void displayReportFile(Fl_Widget *w, void* userdata) ;


class patWinBiogeme : public Fl_Window {

  friend void selectModelFileName(Fl_Widget *w, void* userdata) ;
  friend void selectDataFileName(Fl_Widget *w, void* userdata) ;
  friend void estimateBiogeme(Fl_Widget *w, void* userdata) ;
  friend void simulateBiogeme(Fl_Widget *w, void* userdata) ;
  friend void displayReportFile(Fl_Widget *w, void* userdata) ;


  
 public:
  /**
   */
  patWinBiogeme() ;
  /**
   */
  void addLogMessage(patString m) ;
  /**
   */
  patString getFileNameFromFullPath(patString fullPath, patBoolean keepExtension = patTRUE) ;
  patString getDirectoryFromFullPath(patString fullPath) ;

  /**
   */
  void deactivateButtons() ;

  /**
   */
  void setOutputColor(Fl_Color theCol) ;

 protected:

  Fl_Text_Display* output ;
  Fl_Text_Buffer* outputBuffer ;

  Fl_Button modelFileButton ;
  Fl_Output modelFile ;
  Fl_Button dataFileButton ;
  Fl_Output dataFile ;
  Fl_Output workingDir ;

  Fl_Button reportFileButton ;
  Fl_Output reportFile ;


  Fl_Button estimateButton ;
  Fl_Button simulateButton ;

  patString directory ;
  patString modelSpecFileName ;
  patString dataFileName ;
  patString reportFileName ;  


};
#endif
