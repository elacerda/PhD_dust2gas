#!/bin/bash
benchmark=false
rawdata=false
SFR=false
Zstar=false
radial_profiles=true
TS=$(date +%s)
FNAME="v20_q050.d15a.h5"
BASEFNAME="${FNAME%%.h5}"
BASEDIR="${HOME}/dev/astro/PhD_dust2gas"
SRCDIR="${HOME}/dev/astro/PhD_dust2gas/src/PhD_dust2gas"
RUNDIR="${BASEDIR}/runs"
LOGDIR="${RUNDIR}/log"
H5DIR="${RUNDIR}/hdf5"
ARGSDIR="${RUNDIR}/args"
H5FNAME="${RUNDIR}/${FNAME}"
prog_rawdata="${SRCDIR}/rawdata.py"
prog_SFR="${SRCDIR}/SFR.py"
prog_Zstar="${SRCDIR}/Zstar.py"
prog_radial_profiles="${SRCDIR}/radial_profiles.py"

[ ! -d $LOGDIR ] && mkdir -p ${RUNDIR}/log
[ ! -d $H5DIR ] && mkdir -p ${LOGDIR}/hdf5

if ( $benchmark )
then
    prog_rawdata="-m cProfile -s cumulative ${prog_rawdata}"
    prog_SFR="-m cProfile -s cumulative ${prog_SFR}"
    prog_Zstar="-m cProfile -s cumulative ${prog_Zstar}"
    prog_radial_profiles="-m cProfile -s cumulative ${prog_radial_profiles}"
fi

$rawdata && (
    echo "python $prog_rawdata -H ${H5FNAME} @${ARGSDIR}/centaurus_rawdata.args &> ${LOGDIR}/rawdata-${TS}.log"
    python $prog_rawdata -H ${H5FNAME} @${ARGSDIR}/centaurus_rawdata.args &> ${LOGDIR}/rawdata-${TS}.log
    [ -f ${H5FNAME} ] && (
        echo "cp -v ${H5FNAME} ${H5DIR}/${BASEFNAME}_rawdata.h5"
        cp -v ${H5FNAME} ${H5DIR}/${BASEFNAME}_rawdata.h5
    )
)

$SFR && (
    echo "python $prog_SFR -H ${H5FNAME} @${ARGSDIR}/centaurus_SFR.args &> ${LOGDIR}/SFR-${TS}.log"
    python $prog_SFR -H ${H5FNAME} @${ARGSDIR}/centaurus_SFR.args &> ${LOGDIR}/SFR-${TS}.log
    [ -f ${H5FNAME} ] && (
        echo "cp -v ${H5FNAME} ${H5DIR}/${BASEFNAME}_rawdata_SFR.h5"
        cp -v ${H5FNAME} ${H5DIR}/${BASEFNAME}_rawdata_SFR.h5
    )
)

$Zstar && (
    echo "python $prog_Zstar -H ${H5FNAME} @${ARGSDIR}/centaurus_Zstar.args &> ${LOGDIR}/Zstar-${TS}.log"
    python $prog_Zstar -H ${H5FNAME} @${ARGSDIR}/centaurus_Zstar.args &> ${LOGDIR}/Zstar-${TS}.log
    [ -f ${H5FNAME} ] && (
        echo "cp -v ${H5FNAME} ${H5DIR}/${BASEFNAME}_rawdata_SFR_Zstar.h5"
        cp -v ${H5FNAME} ${H5DIR}/${BASEFNAME}_rawdata_SFR_Zstar.h5
    )
)

$radial_profiles && (
    echo "python $prog_radial_profiles -H ${H5FNAME} @${ARGSDIR}/centaurus_radial_profiles.args &> ${LOGDIR}/radial_profiles-${TS}.log"
    python $prog_radial_profiles -H ${H5FNAME} @${ARGSDIR}/centaurus_radial_profiles.args &> ${LOGDIR}/radial_profiles-${TS}.log
)
