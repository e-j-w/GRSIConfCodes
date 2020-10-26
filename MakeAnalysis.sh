#!/bin/bash
DATADIR=MidasFiles # File Path of Midas Files
FRAGDIR=FragmentTrees # Directory to store Fragment Trees
ANALDIR=AnalysisTrees # Directory to store Analysis Trees

DLEN=${#DATADIR}
 
m=51644 # Only Sort after this run

if [ $# -eq 0 ]
	then
for f in $DATADIR/*.mid; do # Sorts all runs
 g=${f:DLEN+4}
 h=${g:0:${#g}-4} 
 i=${g:0:${#g}-8} 
 FFILE=$FRAGDIR/fragment$h.root
 AFILE=$ANALDIR/analysis$h.root
 CFILE=CalibrationFile.cal # Calibration File
 #echo "Processing run$g $h $i"

if [ ! -f $AFILE ] && [ $i -gt $m ]; # Checks if file exists if not it sorts it  
		then
			echo "File $AFILE does not exist."
			grsisort --recommended $f $CFILE
			mv -f analysis$h.root $ANALDIR
			mv -f fragment$h.root $FRAGDIR

	fi
     
if [ -f $AFILE ]
then 

if [ "$AFILE" -ot "$f" ] && [ $i -gt $m ]; # Checks if the midas file is newer than analysis tree
		then
			echo "File $AFILE exists but is older than $f"
			grsisort --recommended $f $CFILE
			mv -f analysis$h.root $ANALDIR
			mv -f fragment$h.root $FRAGDIR
fi
fi
done 

else 
for f in $DATADIR/run"$@"_*.mid; # Sorts sub runs from a single file

do
 g=${f:DLEN+4}
 h=${g:0:${#g}-4} 
 i=${g:0:${#g}-8} 
 FFILE=$FRAGDIR/fragment$h.root
 AFILE=$ANALDIR/analysis$h.root
 CFILE=CalibrationFile.cal
 echo "Processing run$g "

if [ ! -f $AFILE ];  
		then
			echo "File $AFILE does not exist."
			grsisort --recommended $f $CFILE
			mv -f analysis$h.root $ANALDIR
			mv -f fragment$h.root $FRAGDIR

	fi

if [ -f $AFILE ]
then 
    
if [ "$AFILE" -ot "$f" ] && [ $i -gt $m ];
		then
			echo "File $AFILE exists but is older than $f"
			grsisort --recommended $f $CFILE
			mv -f analysis$h.root $ANALDIR
			mv -f fragment$h.root $FRAGDIR
fi
fi
done
fi
