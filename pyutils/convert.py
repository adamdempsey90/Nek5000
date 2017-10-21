import numpy as np

with open('RB.rea') as f:
    for line in f.readlines():
        if 'p005' in line:
            dad = float(line.split()[0])

fd = 'zd_{:d}.bin'.format(int(dad))
fu = 'zu_{:d}.bin'.format(int(dad))


dd = np.loadtxt('zavg_down.dat')
s = dd.shape
dd = np.hstack((float(s[0]),float(s[1]),dad,dd.ravel()))
dd.tofile(fd)
du = np.loadtxt('zavg_up.dat')
s = du.shape
du = np.hstack((float(s[0]),float(s[1]),dad,du.ravel()))
du.tofile(fu)




