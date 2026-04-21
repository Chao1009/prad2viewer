#!/bin/bash

# Usage: ./run_gain_monitor.sh <run_number> <num_cpus> [subfile_min] [subfile_max]
# Example: ./run_gain_monitor.sh 023545 10
#          ./run_gain_monitor.sh 023545 10 5 20   (only sub-files 5 through 20)

RUN=$1
NCPU=$2
SUBFILE_MIN=$3
SUBFILE_MAX=$4

INPUTDIR="${INPUTDIR:-/data/evio/data}"
OUTPUTDIR="${OUTPUTDIR:-/home/clasrun/prad2_daq/gain_monitoring/gain_monitor_output}"

RUNDIR="${INPUTDIR}/prad_${RUN}"

if [ ! -d "${RUNDIR}" ]; then
    echo "ERROR: run directory not found: ${RUNDIR}"
    exit 1
fi

# Collect all file numbers from evio filenames, sorted numerically
ALL_FILENUMS=($(ls "${RUNDIR}/prad_${RUN}.evio."* 2>/dev/null \
    | sed 's/.*\.evio\.//' \
    | sort -n))

# Filter by sub-file range if third and fourth parameters are given
if [ -n "${SUBFILE_MIN}" ] && [ -n "${SUBFILE_MAX}" ]; then
    echo "Filtering sub-files from ${SUBFILE_MIN} to ${SUBFILE_MAX}"
    FILENUMS=()
    for NUM in "${ALL_FILENUMS[@]}"; do
        N=$(( 10#${NUM} ))
        if [ "${N}" -ge "${SUBFILE_MIN}" ] && [ "${N}" -le "${SUBFILE_MAX}" ]; then
            FILENUMS+=("${NUM}")
        fi
    done
else
    FILENUMS=("${ALL_FILENUMS[@]}")
fi

NFILES=${#FILENUMS[@]}
if [ "${NFILES}" -eq 0 ]; then
    echo "ERROR: no evio files found in ${RUNDIR}"
    exit 1
fi
echo "Found ${NFILES} evio files for run ${RUN}"

mkdir -p "${OUTPUTDIR}"

# Divide files among CPUs
BASE=$(( NFILES / NCPU ))
REM=$(( NFILES % NCPU ))

PIDS=()
PARTFILES=()
IDX=0
for (( i=0; i<NCPU; i++ )); do
    if [ $i -lt $REM ]; then
        COUNT=$(( BASE + 1 ))
    else
        COUNT=$BASE
    fi

    START_NUM=$(( 10#${FILENUMS[$IDX]} ))
    END_NUM=$(( 10#${FILENUMS[$(( IDX + COUNT - 1 ))]} ))
    IDX=$(( IDX + COUNT ))

    PARTFILE="${OUTPUTDIR}/prad_${RUN}_LMS_file_${START_NUM}_${END_NUM}.root"
    PARTFILES+=("${PARTFILE}")

    echo "Job ${i}: file numbers ${START_NUM} to ${END_NUM} (${COUNT} files)"
    ./bin/gain_monitor -r "${RUN}" -s "${START_NUM}" -e "${END_NUM}" \
        -i "${INPUTDIR}" -o "${OUTPUTDIR}" &
    PIDS+=($!)
done

# Wait for all background jobs
echo "Waiting for ${#PIDS[@]} jobs..."
for PID in "${PIDS[@]}"; do
    wait "${PID}"
done
echo "All jobs done."

# Merge partial root files
FINALFILE="${OUTPUTDIR}/prad_${RUN}_LMS.root"
echo "Merging into ${FINALFILE}"
hadd -f "${FINALFILE}" "${PARTFILES[@]}"
rm -f "${PARTFILES[@]}"

# Fit the merged file
echo "Fitting ${FINALFILE}"
./bin/gain_fitter -r "${RUN}" -d "${OUTPUTDIR}"

echo "Done."
