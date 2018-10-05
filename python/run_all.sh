#!/bin/bash

ARGS="cl ch norm cl-nn ch-nn norm-nn"

for RUN in ${ARGS} ; do
    echo Starting ${RUN}
    screen -d -m ./keepup.py ./cellularautomata.py ${RUN}
    echo Started
done
