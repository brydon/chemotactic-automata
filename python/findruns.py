import os
import pickle as pik

dat_files = [x for x in os.listdir(os.getcwd()) if x.endswith(".dat")]
dat_files.sort()

print dat_files

fn = dat_files[-1]

print fn

with open(fn,'rb') as f:
    [tumor, phi] = pik.load(f)
