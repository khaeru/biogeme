ACLOCAL_AMFLAGS = -I m4

noinst_LTLIBRARIES =

if BUILD_BISONBIOGEME
noinst_LTLIBRARIES += libutils.la
endif

if BUILD_GUI
noinst_LIBRARIES = libutilsfltk.a
endif

if BUILD_PYTHON
noinst_LTLIBRARIES += libutilspython.la
endif

if BUILD_STATIC
libutils_la_LDFLAGS = -static
libutilspython_la_LDFLAGS = -static
endif

libutilsfltk_a_SOURCES = $(libutils_la_SOURCES)
libutilsfltk_a_CXXFLAGS = $(FLTK_CXXFLAGS)

libutilspython_la_SOURCES = $(libutils_la_SOURCES)
libutilspython_la_CXXFLAGS = $(PYTHON_CXXFLAGS)


libutils_la_SOURCES = patAbsTime.cc \
    patAngle.cc \
    patCompareCorrelation.cc \
    patErrFileNotFound.cc \
    patErrMiscError.cc \
    patErrNullPointer.cc \
    patError.cc \
    patFileExists.cc \
    patFileNames.cc \
    patFormatRealNumbers.cc \
    patFunction.cc \
    patCorrelation.cc \
    patGeoCoordinates.cc \
    patImportance.cc \
    patLargeNumbers.cc \
    patFileSize.cc \
    patStatistics.cc \
    patLoopTime.cc \
    patMessage.cc \
    patNonParamRegression.cc \
    patNormalCdf.cc \
    patNormalPdf.cc \
    patOutputTable.cc \
    patPValue.cc \
    patPower.cc \
    patPrimeNumbers.cc \
    patPsTricks.cc \
    patString.cc \
    patTimeInterval.cc \
    patTimer.cc \
    patVariables.cc \
    patVersion.cc \
    patQuantiles.cc \
    patLap.cc \
FlexLexer.h \
patAbsTime.h \
patAngle.h \
patBiogemeRun.h \
patClass.h \
patCompare.h \
patCompareCorrelation.h \
patConst.h \
patCorrelation.h \
patDisplay.h \
patErrFileNotFound.h \
patErrMiscError.h \
patErrNullPointer.h \
patErrOutOfRange.h \
patError.h \
patFileExists.h \
patFileNames.h \
patFormatRealNumbers.h \
patFunction.h \
patGenerateCombinations.h \
patGeoCoordinates.h \
patImportance.h \
patIterationBackup.h \
patIterator.h \
patLargeNumbers.h \
patFileSize.h \
patStatistics.h \
patLoop.h \
patLoopTime.h \
patMath.h \
patMessage.h \
patMtl.h \
patNonParamRegression.h \
patNormalCdf.h \
patNormalPdf.h \
patOutputTable.h \
patPValue.h \
patPower.h \
patPrimeNumbers.h \
patPsTricks.h \
patSequenceIterator.h \
patShuffleVector.h \
patStlSetIterator.h \
patStlVectorIterator.h \
patString.h \
patSubvectorIterator.h \
patSubvectorPtrIterator.h \
patThreeStrings.h \
patTimeInterval.h \
patTimer.h \
patType.h \
patUtilFunction.h \
patVariables.h \
patVectorIterator.h \
patVectorPtrIterator.h \
patVersion.h \
trNonLinearAlgo.h \
patQuantiles.h \
patLap.h



nodist_libutils_la_SOURCES = \
patUser.h \
patCompilationDate.h

BUILT_SOURCES = \
patUser.h \
patCompilationDate.h

patUser.h:
	echo "#define __PATUSER \"`id -un` on "`uname`"\""> patUser.h

patCompilationDate.h:
	echo "#define __PATCOMPILATIONDATE \""`date`"\""> patCompilationDate.h

clean-local:
	rm -f patCompilationDate.h patUser.h