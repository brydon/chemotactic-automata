#!/bin/bash

ARGS="cl ch norm"

for RUN in ${ARGS} ; do
    echo Starting ${RUN}
    screen -d -m ./keepup.py ./cellularautomata.py ${RUN}
    echo Started
done
