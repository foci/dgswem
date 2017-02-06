#!/bin/bash

# run directory
RUNDIR=$HOME/performance_testing/square_1/16_hpx_2

cd $RUNDIR || {
    echo "Couldn't cd to RUNDIR = $RUNDIR , exiting"
    exit 
}

TRIAL=1

#output directory prefix for vtune results
OUTPREFIX="results"

# output directory for vtune results
OUTDIR="${OUTPREFIX}_${TRIAL}"

# if directory already exists, 
while [ -z $OUTDIR ]; do 
    echo -n "OUTDIR $OUTDIR already exists, trying "
    TRIAL=$(( $TRIAL + 1 ))
    OUTDIR="$OUTPREFIX_$TRIAL"
    echo "$OUTDIR"
done

# mpirun (or ibrun?)
MPIRUN=mpirun

#path to the executable
APP=$HOME/dgswem-hpx/builds/build_RelWithDebInfo_intel/src/dgswem_lgd_dataflow

# path to vtune amplxe-cl
VTUNE="amplxe-cl"

VTUNE_FLAGS="-data-limit=4000 -r $OUTDIR -collect hotspots -knob enable-user-tasks=true"

APP_FLAGS="--chunksize=100 --n_timesteps=100"

HPX_FLAGS="--hpx:localities=1 --hpx:threads=4 --hpx:ini=hpx.stacks.small_size=0x20000 --hpx:ini=hpx.use_itt_notify=1"
HPX_PERF_COUNTER_FLAGS="--hpx:print-counter=/threads/time/cumulative  --hpx:print-counter=/threads/idle-rate  --hpx:print-counter=/threads/time/cumulative-overhead  --hpx:print-counter=/threadqueue/length --hpx:print-counter-interval=10"



CMD="$MPIRUN $VTUNE ${VTUNE_FLAGS} $APP ${APP_FLAGS} ${HPX_FLAGS} ${HPX_PERF_COUNTER_FLAGS}"

CMD_NO_VTUNE="$MPIRUN $APP ${APP_FLAGS} ${HPX_FLAGS} ${HPX_PERF_COUNTER_FLAGS}"

echo "CMD = $CMD"
echo "CMD_NO_VTUNE = $CMD_NO_VTUNE"


# GDB command:
# run --chunksize=100 --n_timesteps=100 --hpx:localities=1 --hpx:threads=1 --hpx:ini=hpx.stacks.small_size=0x20000 --hpx:ini=hpx.use_itt_notify=1 
