#!/usr/bin/env python
from subprocess import check_call
import numpy as np

for i, f in enumerate(np.linspace(0.75, 5, 200)):
    print "Starting run", i, f
    check_call(["./cupgame-nowait", "512", "%f" % f, "0.9", "/media/scratch/cupgame/loop/temp_%04d.npy" % i])
