#!/bin/bash
RAWDATA=false
#RAWDATA=true
#SFR=true
SFR=false
Zstar=true
#Zstar=false
TS=$(date +%s)
FNAME="v20_q050.d15a.h5"
BASEFNAME="${FNAME%%.h5}"

[ ! -d "log" ] && mkdir log
[ ! -d "hdf5" ] && mkdir hdf5

$RAWDATA && (
    echo "python -m cProfile -s cumulative rawdata.py @centaurus_rawdata.args &> log/rawdata-${TS}.log"
    python -m cProfile -s cumulative rawdata.py @centaurus_rawdata.args &> log/rawdata-${TS}.log
    echo "cp -v $FNAME hdf5/${BASEFNAME}_rawdata.h5"
    cp -v $FNAME hdf5/${BASEFNAME}_rawdata.h5
)

$SFR && (
    echo "python -m cProfile -s cumulative SFR.py -H $FNAME @centaurus_SFR.args &> log/SFR-${TS}.log"
    python -m cProfile -s cumulative SFR.py -H $FNAME @centaurus_SFR.args &> log/SFR-${TS}.log
    echo "cp -v $FNAME hdf5/${BASEFNAME}_rawdata_SFR.h5"
    cp -v $FNAME hdf5/${BASEFNAME}_rawdata_SFR.h5
)

$Zstar && (
    echo "python -m cProfile -s cumulative Zstar.py -H $FNAME @centaurus_Zstar.args &> log/Zstar-${TS}.log"
    python -m cProfile -s cumulative Zstar.py -H $FNAME @centaurus_Zstar.args &> log/Zstar-${TS}.log
    echo "cp -v $FNAME hdf5/${BASEFNAME}_rawdata_SFR_Zstar.h5"
    cp -v $FNAME hdf5/${BASEFNAME}_rawdata_SFR_Zstar.h5
)

#$Zstar && (
#    echo " python Zstar.py -H $FNAME @centaurus_Zstar.args &> log/Zstar-${TS}.log"
#    python Zstar.py -H $FNAME @centaurus_Zstar.args &> log/Zstar-${TS}.log
#    echo "cp -v $FNAME hdf5/${BASEFNAME}_rawdata_SFR_Zstar.h5"
#    cp -v $FNAME hdf5/${BASEFNAME}_rawdata_SFR_Zstar.h5
#)
