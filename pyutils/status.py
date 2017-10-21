import numpy as np
import matplotlib.pyplot as plt
try:
    execfile('/projects/p20850/amd616/Nek5000/pyutils/nek.py')
except NameError:
    with open('/projects/p20850/amd616/Nek5000/pyutils/nek.py') as f:
        exec(f.read())


fld,_,_,_,_ = loadit()
print(fld.nt)
plt.ion()
fld.tss(savefig='tss.png')
if fld.nt > 50:
    fld.summary(average=True,window=30,savefig='summary.png')
elif fld.nt > 100:
    fld.summary(average=True,window=70,savefig='summary.png')
elif fld.nt > 130:
    fld.summary(average=True,window=100,savefig='summary.png')
else:
    fld.summary(average=False,savefig='summary.png')
plt.show(block=True)


