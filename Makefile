CPPFLAGS = -std=c++0x

all: SetupConfFile CustomChan ConstructCalibrationFile ChanCheck Calibrate Calibrate_LaBr ResCheck Alphacal EffCheck Efficiency

SetupConfFile: Conf_Setup.cxx functions.h CalibrationParameters.h
	g++ Conf_Setup.cxx $(CPPFLAGS) -o SetupConfFile
CustomChan: CustomChannelSetup.cxx
	g++ CustomChannelSetup.cxx $(CPPFLAGS) -o CustomChan
ConstructCalibrationFile: CalFileConstructor.C
	g++ CalFileConstructor.C $(CPPFLAGS) -o ConstructCalibrationFile 
ChanCheck: chan_check.cxx
	g++ chan_check.cxx $(CPPFLAGS) -I$(GRSISYS)/include -L$(GRSISYS)/libraries `grsi-config --cflags --all-libs` `root-config --cflags --libs` -lTreePlayer -lSpectrum -o ChanCheck
Calibrate: calib-ge.cxx
	g++ calib-ge.cxx $(CPPFLAGS) -o Calibrate
Calibrate_LaBr: LaBr/calib-labr.cxx
	g++ LaBr/calib-labr.cxx $(CPPFLAGS) -o Calibrate_LaBr
ResCheck: res_check.cxx
	g++ res_check.cxx $(CPPFLAGS) -I$(GRSISYS)/include -L$(GRSISYS)/lib -I$(GRSISYS)/GRSIData/include -L$(GRSISYS)/GRSIData/lib `grsi-config --cflags --all-libs --GRSIData-libs` `root-config --cflags --libs` -lTreePlayer -lMathMore -lSpectrum -lMinuit -lPyROOT -o ResCheck 
Alphacal: AlphaCalibration.c
	g++ AlphaCalibration.c $(CPPFLAGS) -I$(GRSISYS)/include -L$(GRSISYS)/lib `grsi-config --cflags --all-libs --GRSIData-libs` -I$(GRSISYS)/GRSIData/include -L$(GRSISYS)/GRSIData/lib `root-config --cflags --libs` -lTreePlayer -lMathMore -lSpectrum -lMinuit -lPyROOT -o Alphacal
EffCheck: eff_check.cxx
	g++ eff_check.cxx -std=c++0x -I$(GRSISYS)/include -L$(GRSISYS)/lib `grsi-config --cflags --all-libs --GRSIData-libs` -I$(GRSISYS)/GRSIData/include -L$(GRSISYS)/GRSIData/lib `root-config --cflags --libs` -lTreePlayer -lMathMore -lSpectrum -lMinuit -lPyROOT -o EffCheck 
Efficiency: eff-ge.cxx
	g++ eff-ge.cxx $(CPPFLAGS) -o Efficiency

clean:
	rm -rf *~ SetupConfFile CustomChan ConstructCalibrationFile ChanCheck Calibrate Calibrate_LaBr ResCheck Alphacal EffCheck Efficiency *tmpdatafile*
