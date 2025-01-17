LINKFLAGS = -Wl,--no-as-needed $(shell root-config --cflags --libs --glibs) -lSpectrum -lMinuit -lGuiHtml -lTreePlayer -lTMVA -L/opt/local/lib -lX11 -lXpm
CXXFLAGS = -O2 -Wl,--copy-dt-needed-entries
GRSISORT = -L/opt/local/lib -lX11 -lXpm $(shell grsi-config --cflags --all-libs --GRSIData-libs) -I$(GRSISYS)/GRSIData/include

all: SetupConfFile CustomChan ConstructCalibrationFile ChanCheck Calibrate Calibrate_LaBr ResCheck Alphacal EffCheck Efficiency

SetupConfFile: Conf_Setup.cxx functions.h CalibrationParameters.h
	g++ Conf_Setup.cxx $(CXXFLAGS) -o SetupConfFile
CustomChan: CustomChannelSetup.cxx
	g++ CustomChannelSetup.cxx $(CXXFLAGS) -o CustomChan
ConstructCalibrationFile: CalFileConstructor.C
	g++ CalFileConstructor.C $(CXXFLAGS) -o ConstructCalibrationFile 
ChanCheck: chan_check.cxx
	g++ chan_check.cxx $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o ChanCheck
Calibrate: calib-ge.cxx
	g++ calib-ge.cxx $(CXXFLAGS) -o Calibrate
Calibrate_LaBr: LaBr/calib-labr.cxx
	g++ LaBr/calib-labr.cxx $(CXXFLAGS) -o Calibrate_LaBr
ResCheck: res_check.cxx
	g++ res_check.cxx $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o ResCheck 
Alphacal: AlphaCalibration.c
	g++ AlphaCalibration.c $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o Alphacal
EffCheck: eff_check.cxx
	g++ eff_check.cxx $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o EffCheck 
Efficiency: eff-ge.cxx
	g++ eff-ge.cxx $(CXXFLAGS) -o Efficiency

clean:
	rm -rf *~ SetupConfFile CustomChan ConstructCalibrationFile ChanCheck Calibrate Calibrate_LaBr ResCheck Alphacal EffCheck Efficiency *tmpdatafile*
