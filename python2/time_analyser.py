#!/software/.admin/bins/bin/python2.7

import sys
import numpy as np

def analyse(fn):
    with open(fn,'r') as f:
        _text = f.readlines()

    dispersed, lcd, dumped = [], [], []

    for line in _text:
        if line.strip().startswith("dispersed in"):
            dispersed.append(float(line[line.find("in")+2:]))
        elif line.strip().startswith("life cycled in"):
            lcd.append(float(line[line.find("in")+2:]))
        elif line.strip().startswith("dumped in"):
            dumped.append(float(line[line.find("in")+2:]))

    return np.array(dispersed), np.array(lcd), np.array(dumped)


if __name__ == "__main__":
    disp, life, dump = analyse(sys.argv[-1])

    print "Life cycle avg time", np.mean(life)
    print "Dispersion avg time", np.mean(disp)
    print "Dumping avg time", np.mean(dump)

    print "Overhead avg time", np.mean(dump+life)
    print [np.round(x,2) for x in dump+life]
    print "Frame avg time", np.mean(dump+life+disp)
