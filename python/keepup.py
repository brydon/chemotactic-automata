#!/software/.admin/bins/bin/python2.7
import sys
import time
import subprocess
from numpy.random import shuffle

print "Starting keepup"

"""
Keep a process up and running

Use it like this:
./keepup.py ./script.py args
"""

MAX = 100

if __name__ == "__main__":
    cmd = ' '.join(sys.argv[1:])

    def start_subprocess():
        return subprocess.Popen(cmd, shell=True)

    p = start_subprocess()
    i = 0

    while True:
        res = p.poll()
        if res is not None:
            if i == MAX - 1: # Counting starts at 0
                print "Should be done"
                break
            else:
                i += 1
            p = start_subprocess()

            print p.pid, 'was killed, restarting it with cmd ' + cmd

        time.sleep(0.5)
