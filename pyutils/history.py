class History(Parameters):
    def __init__(self,n=6,directory='',fname = 'RB.his'):
        if len(directory) > 1:
            if directory[-1] != '/':
                directory += '/'

        Parameters.__init__(self,directory)
        dat = np.loadtxt(directory + fname)


        with open(directory + 'hpts.in') as f:
            lines = f.readlines()

        self.nhis = int(lines[0].strip())

        self.pts = np.zeros((self.nhis,3))

        for i in range(self.nhis):
            self.pts[i,:] = [float(_) for _ in lines[i+1].strip().split()]



        self.n = n
        self.nt = np.size(dat)//(self.n*self.nhis)

        dat = dat.ravel().reshape(self.nt,self.nhis,self.n)
        self.t = dat[:,0,0]
        self.vx = dat[:,:,1].copy()
        self.vy = dat[:,:,2].copy()
        self.vz = dat[:,:,3].copy()
        self.p = dat[:,:,4].copy()
        self.T = dat[:,:,5].copy()

        self.s = self.T.copy()

        for i in range(self.nhis):
            self.s[:,i] -= self.delad*self.pts[i,2]

    def plot(self,q,i='all',ax=None,fig=None,savefig=None,logx=False,logy=False,**kargs):

        if ax is None:
            fig,ax = plt.subplots(figsize=(8,6))

        try:
            if i == 'all':
                ax.plot(self.t,getattr(self,q),**kargs)
            else:
                ax.plot(self.t,getattr(self,q)[:,i],**kargs)
        except TypeError:
            if i == 'all':
                ax.plot(self.t,q,**kargs)
            else:
                ax.plot(self.t,q[:,i],**kargs)

        if logx:
            ax.set_xscale('log')
        if logy:
            ax.set_yscale('log')

        ax.set_xlabel('$t$')

        ax.minorticks_on()
        fig.tight_layout()
        return fig,ax




