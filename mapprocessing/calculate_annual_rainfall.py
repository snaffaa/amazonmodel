'''
compute the annual rainfall for the selected year.
'''

from pcraster import *
from pcraster.framework import *
import os, sys
import datetime

workdir = os.getcwd()
inputdir = os.path.join(workdir, 'amazon_fine')
mapdir= os.path.join(workdir,'mapinput')
txtdir = os.path.join(workdir, 'txtdata')


clone_filename= os.path.join(mapdir, 'amazon_5min_mask.map')
setclone(clone_filename)
Pannual_test= scalar(0)
paths = [p for p in os.listdir(inputdir) if p.startswith('Pre90')]
for path in paths:
    Pmonth = readmap(os.path.join(inputdir,path))
    Pannual_test+=Pmonth
mv = -999.9
Pannual_test= cover(Pannual_test , mv)
report(Pannual_test, os.path.join(inputdir, "Pan_90.map"))
