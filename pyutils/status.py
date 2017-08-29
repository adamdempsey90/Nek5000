import numpy as np
import matplotlib.pyplot as plt
execfile('/projects/p20850/amd616/Nek5000/pyutils/nek.py')

dd = np.loadtxt('zavg_down.dat')
du = np.loadtxt('zavg_up.dat')
dat = dd + du
dat[:,0] /= 2
dat[:,1] /= 2

fld = Zavg(dat=dat)
plt.ion()
fld.tss()
if fld.nt > 50:
    fld.summary(average=True,window=30)
elif fld.nt > 100:
    fld.summary(average=True,window=70)
elif fld.nt > 130:
    fld.summary(average=True,window=100)
else:
    fld.summary(average=False)
plt.show(block=True)


