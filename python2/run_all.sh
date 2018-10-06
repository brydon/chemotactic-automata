#!/bin/bash

ARGS="cl ch norm cl-nn ch-nn norm-nn"

for RUN in ${ARGS} ; do
    echo Starting ${RUN}
    STR=$( date +%b%d%y_%H%M%S )
    screen -d -m ./keepup.py ./cellularautomata.py run_${RUN}_${STR} ${RUN}
    echo Started
done
