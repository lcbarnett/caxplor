#!/bin/bash

scriptname=$(basename $0 .sh)
echo -e "\n*** Running batch script '"$scriptname"' ***\n"

codedir=$HOME/git/caxplor
currdir=$(pwd -P)

flamres=10
nthreads=12
nfpert=400
odir=$currdir

for i in $(seq 1 $flamres); do

	# the log file
	logfile=$currdir/caddr\_$i.log

	# run .caxplor
	$codedir/caxplor ddr -jobidx $i -flamres $flamres -nthreads $nthreads -nfpert $nfpert -odir $odir > $logfile < /dev/null 2>&1

done
