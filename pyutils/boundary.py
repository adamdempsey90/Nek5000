import numpy as np
from scipy.optimize import fsolve
import sys



d = -abs(float(sys.argv[1]))
a = float(sys.argv[2])
n2 = -abs(float(sys.argv[3]))
st = float(sys.argv[4])



if st > 5:
    # Assume this is entropy and we want overshoot
    func = lambda z: -2*d + n2 + 2*(d+n2)*z + 2*d/a * np.log((1+a)/(1+a*z)) - st
    print(fsolve(func,.05,factor=.01)[0])
else:
    # Assume this is overshoot and we want entropy
    z = st
    print(-2*d + n2 + 2*(d+n2)*z + 2*d/a * np.log((1+a)/(1+a*z)))
