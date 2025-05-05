# GRSIConfCodes

Codes for setting up the TIGRESS/GRIFFIN DAQ systems

Original author: S. Gillespie

Revisions (Oct 2020 - present): J. Williams, D. Yates, Y. Zhu

# Running

Requires GRSISort to be installed, with the `$GRSISYS` environment variable configured.

Instructions on how to use these codes are on the [GRSI wiki](https://grsi.wiki.triumf.ca/index.php/Stephens_page_of_things) or [here](doc/stephens_page_of_things.pdf).

# Documentation for specific codes

## AlphaCalibration.c
Auto calibration code for charged partile detector with triple alpha source;<span style="color:red"> GRSISort Required.</span>

0. **Compiling Command (1st line in the code txt):** 

```
g++ AlphaCalibration.c -Wl,--no-as-needed `root-config --cflags --libs --glibs` -lSpectrum -lMinuit -lGuiHtml -lTreePlayer -lTMVA -L/opt/local/lib -lX11 -lXpm -O2 -Wl,--copy-dt-needed-entries -L/opt/local/lib -lX11 -lXpm `grsi-config --cflags --all-libs --GRSIData-libs` -I$GRSISYS/GRSIData/include -o bin/Alphacal
```

1. **Input: FragmentTree + CalibrationFile + starting CH + ending CH;**

2. **Calibration Math Formula: Pu+Am+Cm.**
    - Edit lin 328 to `TF1 *fc = tasf(hist, Form("fc_CH%i",ich), min,max,"cl");` for **Gd+Th+Cm**;
3. **Output:**</br>
    - **Calibration.txt**: includes two array, gain and offset;</br>
    - **Hist.root:** includes calibrated summary TH2 for calibration quick check;</br>
    - **Values Printed on Screen:** FWHM of three alpha peaks after calibration.</br>

<span style="color:red">Return FWHM = -1, GAIN = 1, OFFSET = 0, **if the channel is empty!**</span>

# Dev history

These codes were originally developed by S. Gillespie as part of his postdoctoral work.  This repository collects the various DAQ setup and calibration codes in a single location, and may host future updates to the codes.  For now, it is recommended to use the original versions (specified here: https://grsi.wiki.triumf.ca/index.php/Stephens_page_of_things) unless you *really* know what you're doing.