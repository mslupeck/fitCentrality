#!/bin/bash

START=$1
STOP=$2
if [ -z "$STOP" ]; then STOP=$((START+1)); fi #if stop run is unset assume only single run is to be handled

REMOTEBASEFOLDER="/alice/cern.ch/user/m/mslupeck/2018-07-fitCentrRes-B500_FullImpact"
LOCALBASEFOLDER="/home/mss/workspace/fitCentrality/inputs/B500_FullImpact"

for (( i=$START; i<$STOP; i++ ))
do
   printf -v NJOB "%0.3i" $i
   echo ${NJOB}
   NFOLDER=${LOCALBASEFOLDER}/${NJOB}
   if [ ! -d ${NFOLDER} ]; then mkdir ${NFOLDER}; fi
   alien_cp -v alien:${REMOTEBASEFOLDER}/outputs/246392/${NJOB}/log_archive.zip ${NFOLDER}
   alien_cp -v alien:${REMOTEBASEFOLDER}/outputs/246392/${NJOB}/root_archive.zip ${NFOLDER}
   unzip ${NFOLDER}/log_archive.zip -d ${NFOLDER}
   unzip ${NFOLDER}/root_archive.zip -d ${NFOLDER}
done

