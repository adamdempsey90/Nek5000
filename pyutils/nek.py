import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad,cumtrapz
from scipy.interpolate import interp1d
import copy


class Parameters():
    def __init__(self,directory=''):
        """ 118 parameters in nek """
        if len(directory) > 1:
            if directory[-1] != '/':
                directory += '/'

        with open(directory + 'RB.rea') as f:
            lines = f.readlines()


        names = ['p{:03d}'.format(i) for i in range(1,119)]
        self.params = {}
        for line in lines:
            for name in names:
                if name in line:
                    self.params[name] = float(line.split()[0])

        with open(directory + 'RB.box') as f:
            lines = f.readlines()

        self.nelx, self.nely, self.nelz = (-int(x) for x in lines[9].split()[:3])

        with open(directory + 'SIZE') as f:
            lines = f.readlines()

        for line in lines:
            if 'lx1=' in line:
                line = line.split()[1].split(',')
                break

        lx1 = line[0].split('=')[-1]
        ly1 = line[1].split('=')[-1]
        lz1 = line[2].split('=')[-1]


        try:
            lx1 = int(lx1)
            lx1_flag = True
        except ValueError:
            lx1_flag=  False
        try:
            ly1 = int(ly1)
            ly1_flag = True
        except ValueError:
            ly1_flag=  False
        try:
            lz1 = int(lz1)
            lz1_flag = True
        except ValueError:
            lz1_flag=  False

        if not lx1_flag:
            if lx1 == 'ly1':
                lx1 = ly1
            elif lx1 == 'lz1':
                lx1 = lz1
            else:
                print("Can't set lx1!")
        if not ly1_flag:
            if ly1 == 'lx1':
                ly1 = lx1
            elif ly1 == 'lz1':
                ly1 = lz1
            else:
                print("Can't set ly1!")
        if not lz1_flag:
            if lz1 == 'lx1':
                lz1 = lx1
            elif lz1 == 'ly1':
                lz1 = ly1
            else:
                print("Can't set lz1!")

        self.lx1 = lx1
        self.ly1 = ly1
        self.lz1 = lz1
        self.nx = self.lx1 * self.nelx
        self.ny = self.ly1 * self.nely
        self.nz = self.lz1 * self.nelz

        self.grid_params = ['lx1','ly1','lz1','nx','ny','nz']

        for key,val in self.params.items():
            setattr(self,key, val)


        # Special parameters

       # self.named_params = { 'p002': 'nu',
       #            'p008': 'k0',
       #            'p075': 'omega',
       #            'p076': 'theta',
       #            'p077': 'nk',
       #            'p078': 'ks',
       #            'p079': 'delta',
       #            'p080': 'delad'}
        self.named_params = { 'p002': 'nu',
                   'p008': 'k0',
                   'p003': 'omega',
                   'p004': 'theta',
                   'p009': 'ks',
                   'p006': 'delta',
                   'p005': 'delad'}

        for key,val in self.named_params.items():
            setattr(self,val,getattr(self,key))



        self.kapfunc = lambda x: np.exp( np.log(self.ks) + .5*np.log(self.ks/self.k0)*(np.tanh((x-1)/self.delta -2.) - np.tanh(x/self.delta+2.)))
        self.dkapfunc = lambda x: 0
        self.intkapfunc = lambda x: 0

        self.delad *= -1

#        if int(self.kprof) == 0:
#            self.kapfunc = lambda x: self.k0 * ( (1+x)**(-self.nk) + (2-x)**(-self.nk))/(1 + 2.**(-self.nk))
#            self.dkapfunc = lambda x: self.k0*self.nk * ( (2+x)**(-self.nk-1) - (1 + x)**(-self.nk-1))/(1 + 2.**(-self.nk))
#            self.intkapfunc = lambda x: np.array([ quad(lambda x1: 1./self.kapfunc(x1),0,i)[0] / (-self.k0) for i in x])
#        elif int(self.kprof) == 1:
#            self.kapfunc = lambda x: 10**( np.log10(self.ks) + .5*np.log10(self.ks/self.k0)*(np.tanh((x-1)/self.delta -2.) - np.tanh(x/self.delta+2.)))
#            self.dkapfunc = lambda x: 0
#            self.intkapfunc = lambda x: 0
#        else:
#            self.kapfunc = lambda x: self.k0
#            self.dkapfunc = lambda x: 0
#            self.intkapfunc = lambda x: self.k0*x
#





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









class Tavg(Parameters):
    cols = ['t','vz2','vp2','KE','enstrophy',
            'visc_diss_1','visc_diss_2','visc_diss_3',
            'therm_diss','dTdz','Fk','Fc','Fe',
            'vrms','Nu','Ftot','dsdz','Re']
    def __init__(self,directory='',fname = 'fort.50'):
        if len(directory) > 1:
            if directory[-1] != '/':
                directory += '/'
        try:
            dat = np.loadtxt(directory + fname)
            self.read_nek(directory, dat)
        except FileNotFoundError as e:
            raise e
        except Exception:
            self.readh5(fname)

    def read_nek(self,directory, dat):
       Parameters.__init__(self,directory)
       if len(dat) == 0:
           print('File has nothing in it!')
           return None

       self.t = dat[:,0]
       self.vz2 = dat[:,1]
       self.vp2 = dat[:,2]
       self.KE = dat[:,3]
       self.enstrophy = dat[:,4]
       self.visc_diss_1 = dat[:,5]
       self.visc_diss_2 = dat[:,6]
       self.visc_diss_3 = dat[:,7]
       self.therm_diss = dat[:,8]
       self.dTdz = dat[:,9]
       self.Fk = dat[:,10]
       self.Fc = dat[:,11]
       self.Fe = dat[:,12]
       self.Nu = self.Fc/self.Fk
       self.Ftot = self.Fk + self.Fc
       self.dsdz = self.dTdz - self.delad
       self.vrms = np.sqrt(self.vz2)
       self.Re = self.vrms/self.nu

    def readh5(self, f):
        for key,val in f.attrs.items():
            setattr(self, key,val)
        names = []
        f.visit(names.append)
        for name in names:
            setattr(self,name,f[name][...])

    def dump(self, f):
        for key,val in self.params.items():
            try:
                f.attrs[key] = val
            except RuntimeError:
                print('{} key already created! Skipping...'.format(key))
                pass
        for key in self.grid_params:
            try:
                f.attrs[key] = getattr(self,key)
            except RuntimeError:
                print('{} key already created! Skipping...'.format(key))
                pass
        for key,val in self.named_params.items():
            try:
                f.attrs[val] = getattr(self,val)
            except RuntimeError:
                print('{} key already created! Skipping...'.format(key))
                pass
        for key in self.cols:
            try:
                f.create_dataset(key,data = getattr(self,key))
            except RuntimeError:
                print('{} key already created! Skipping...'.format(key))
                pass

    def turb_plot(self,fig=None,ax=None,savefig=None,logx=False,ylims=None,xlims=None,norm=1,**kargs):
        if ax is None:
            fig,ax = plt.subplots(figsize=(8,6))

        if xlims is not None:
            ind = (self.t>=xlims[0])&(self.t<=xlims[1])
        else:
            ind= Ellipsis

        k = self.KE[ind]/norm
        eps = self.visc_diss_1[ind]/norm
        t = self.t[ind]

        ax.plot(t,np.sqrt(k),label='Turbulent Velocity, $k^{1/2}$',**kargs)
        ax.plot(t,k**1.5/eps,label='Large Eddy Size, $k^{3/2} / \\epsilon$',**kargs)
        ax.plot(t,k/eps,label='Eddy Turnover Time, $k / \\epsilon$',**kargs)
        ax.plot(t,k**2/eps,label='Eddy Viscosity, $k^2 / \\epsilon$',**kargs)

        ax.set_xlabel('$t$',fontsize=20)

        ax.tick_params(labelsize=15)
        ax.minorticks_on()

        if ylims is not None:
            ax.set_ylim(ylims)

        if logx:
            ax.set_xscale('log')

        ax.legend(loc='best')

        if savefig is not None:
            fig.savefig(savefig, bbox_inches='tight')

        return fig,ax

    def semilogx(self,q,**kargs):
        return self.plot(q, logx=True, logy=False, **kargs)
    def semilogy(self,q,**kargs):
        return self.plot(q, logx=False, logy=True, **kargs)
    def loglog(self,q,**kargs):
        return self.plot(q, logx=True, logy=True, **kargs)



    def plot(self,q,ax=None,fig=None,shift=0,norm=1,window=100,logx=False,logy=False,savefig=None,**kargs):
        if ax is None:
            fig,ax = plt.subplots(figsize=(8,6))
        try:
            avg,st = window_avg2(q,window=window)
            line, = ax.plot(self.t,avg,**kargs)
            ax.fill_between(self.t, avg-st, avg+st, facecolor=line._color,alpha=.1)

        except TypeError:
            try:
                dat = getattr(self,q)
            except AttributeError:
                print('{} not a valid choice!'.format(q))

            avg,st = window_avg((dat-shift)/norm,window)
            line, = ax.plot(self.t,avg,**kargs)
            ax.fill_between(self.t, avg-st, avg+st, facecolor=line._color,alpha=.1)
            ax.set_ylabel(q,fontsize=15)
        ax.set_xlabel('$t$',fontsize=20)
        ax.minorticks_on()
        if logx:
            ax.set_xscale('log')
        if logy:
            ax.set_yscale('log')
        if savefig is not None:
            fig.savefig(savefig)
        return fig,ax
    def __add__(self,fld):
        newfld = copy.deepcopy(self)
        newfld.vz2 += fld.vz2
        newfld.vp2 += fld.vp2
        newfld.KE += fld.KE
        newfld.enstrophy += fld.enstrophy
        newfld.visc_diss_1 += fld.visc_diss_1
        newfld.visc_diss_2 += fld.visc_diss_2
        newfld.visc_diss_3 += fld.visc_diss_3
        newfld.therm_diss += fld.therm_diss
        newfld.dTdz += fld.dTdz
        newfld.Fk += fld.Fk
        newfld.Fc += fld.Fc
        newfld.Fe += fld.Fe

        newfld.Nu = newfld.Fc/newfld.Fk
        newfld.Ftot = newfld.Fk + newfld.Fc
        newfld.dsdz = newfld.dTdz - newfld.delad
        newfld.vrms = np.sqrt(newfld.vz2)
        newfld.Re = newfld.vrms/newfld.nu
        return newfld


class Zavg(Parameters):
    cols = ['t',
        'z',
        'vx',
        'vy',
        'vz',
        'vx2',
        'vy2',
        'vz2',
        'KE',
        'dTdz',
        'T',
        'T2',
        'Fk',
        'Fc',
        'Fe',
        'kap',
        'vort',
        'area',
        'Fb',
        'Nu',
        'Ftot',
        'vrms',
        's',
        'dsdz',
        'sk',
        'vortx',
        'vorty',
        'visc3',
        'therm',
        'P',
        'dzP',
        'Pflux',
        'TdzP',
        'trip1',
        'trip2',
        'rey12',
        'rey13',
        'rey23',
        'vflux',
        'keflux',
        'tflux',
        'thermflux',
        'Fcflux1',
        'Fcflux2',
        'visc1',
        'visc2']
    def __init__(self,directory='',fname='zavg.dat',dat=None):
        if len(directory) > 1:
            if directory[-1] != '/':
                directory += '/'
        try:
            if dat is None:
                dat = np.loadtxt(directory + fname)
            self.read_nek(directory, dat)
        except FileNotFoundError as e:
            raise e
        except Exception:
            self.readh5(fname)
    def read_nek(self, directory, dat):
        Parameters.__init__(self,directory)
        self.t  = np.unique(dat[:,0])
        nt  = len(self.t)
        nz = dat.shape[0]//nt
        nf = dat.shape[1]
        self.nt = nt
        self.nz = nz
        self.nf = nf
        dat = dat.reshape(nt,nz,nf)
        dat = dat[:,:,1:]
        print(dat.shape)
        self.z = dat[0,:,0]
        self.vx = dat[:,:,1].T
        self.vy = dat[:,:,2].T
        self.vz = dat[:,:,3].T
        self.vx2 = dat[:,:,4].T
        self.vy2 = dat[:,:,5].T
        self.vz2 = dat[:,:,6].T
        self.KE = dat[:,:,7].T
        self.dTdz = dat[:,:,8].T
        if nf < 17:
            self.T = dat[:,:,9].T
            self.Fk = dat[:,:,10].T
            self.Fc = dat[:,:,11].T
            self.Fe = dat[:,:,12].T
            self.kap = dat[:,:,13].T
            self.vort = dat[:,:,14].T
            self.T2 = np.zeros(dat[:,:,9].T.shape)
        else:
            self.T = dat[:,:,9].T
            self.T2 = dat[:,:,10].T
            self.Fk = dat[:,:,11].T
            self.Fc = dat[:,:,12].T
            self.Fe = dat[:,:,13].T
            self.kap = dat[:,:,14].T
            self.vort = dat[:,:,15].T
            self.area = np.zeros(dat[:,:,15].T.shape)
            if nf > 17:
                self.area = dat[:,:,16].T
            if nf > 18:
                self.vortx = dat[:,:,17].T
                self.vorty = dat[:,:,18].T
                self.visc3 = dat[:,:,19].T
                self.therm = dat[:,:,20].T
                self.P = dat[:,:,21].T
                self.dzP = dat[:,:,22].T
                self.Pflux = dat[:,:,23].T
                self.TdzP = dat[:,:,24].T
                self.trip1 = dat[:,:,25].T
                self.trip2 = dat[:,:,26].T
                self.rey12 = dat[:,:,27].T
                self.rey13 = dat[:,:,28].T
                self.rey23 = dat[:,:,29].T
                self.vflux = dat[:,:,30].T
                self.keflux = dat[:,:,31].T
                self.tflux = dat[:,:,32].T
                self.thermflux = dat[:,:,33].T
                self.Fcflux1 = dat[:,:,34].T
                self.Fcflux2 = dat[:,:,35].T
                self.visc1 = dat[:,:,36].T
                self.visc2 = dat[:,:,37].T
            else:
                self.vortx = np.zeros(self.Fc.shape)
                self.vorty = np.zeros(self.Fc.shape)
                self.visc3 = np.zeros(self.Fc.shape)
                self.therm = np.zeros(self.Fc.shape)
                self.P = np.zeros(self.Fc.shape)
                self.dzP = np.zeros(self.Fc.shape)
                self.Pflux = np.zeros(self.Fc.shape)
                self.TdzP = np.zeros(self.Fc.shape)
                self.trip1 = np.zeros(self.Fc.shape)
                self.trip2 = np.zeros(self.Fc.shape)
                self.rey12 = np.zeros(self.Fc.shape)
                self.rey13 = np.zeros(self.Fc.shape)
                self.rey23 = np.zeros(self.Fc.shape)
                self.vflux = np.zeros(self.Fc.shape)
                self.keflux = np.zeros(self.Fc.shape)
                self.tflux = np.zeros(self.Fc.shape)
                self.thermflux = np.zeros(self.Fc.shape)
                self.Fcflux1 = np.zeros(self.Fc.shape)
                self.Fcflux2 = np.zeros(self.Fc.shape)
                self.visc1 = np.zeros(self.Fc.shape)
                self.visc2 = np.zeros(self.Fc.shape)

        if np.all(self.thermflux == 0):
            self.thermflux = np.zeros(self.T2.shape)
            dz = np.zeros(self.z.shape)
            dz[0] = self.z[1]-self.z[0]
            dz[-1] = self.z[-1]-self.z[-2]
            dz[1:-1] = self.z[2:]-self.z[:-2]
            self.thermflux[0,:] = -.5*(self.T2[1,:]-self.T2[0,:])/dz[0]
            self.thermflux[-1,:] = -.5*(self.T2[-1,:]-self.T2[-2,:])/dz[-1]
            self.thermflux[1:-1,:] = -.5*(self.T2[2:,:]-self.T2[:-2,:])/dz[1:-1,np.newaxis]

        self.Fb = -self.vz*self.delad*self.z[:,np.newaxis]
        self.Nu = self.Fc/self.Fk

        self.Ftot = self.Fk + self.Fc
        self.vrms = np.sqrt(self.vz2)

        self.s = self.T - self.delad*self.z[:,np.newaxis]
        self.dsdz = self.dTdz - self.delad
        self.sk = cumtrapz(-1./self.kap[:,-1]-self.delad,x=self.z,initial=0)
        tmp = cumtrapz((-1./self.kap[:,-1]-self.delad)[::-1],x=self.z[::-1],initial=0)[::-1]
        self.sk_t = self.s[-1,:] + tmp[:,np.newaxis]
        self.N2k = -1./self.kap[:,0] - self.delad
        #kfunc = interp1d(self.z, self.kap[:0])
        #self.Fconduc = np.array([ quad(kfunc, self.z[0],i)[0] for i in self.z])
    def readh5(self, f):
        for key,val in f.attrs.items():
            setattr(self, key,val)
        names = []
        f.visit(names.append)
        for name in names:
            setattr(self,name,f[name][...])
    def dump(self, f):
        for key,val in self.params.items():
            try:
                f.attrs[key] = val
            except RuntimeError:
                print('{} key already created! Skipping...'.format(key))
                pass
        for key in self.grid_params:
            try:
                f.attrs[key] = getattr(self,key)
            except RuntimeError:
                print('{} key already created! Skipping...'.format(key))
                pass
        for key,val in self.named_params.items():
            try:
                f.attrs[val] = getattr(self,val)
            except RuntimeError:
                print('{} key already created! Skipping...'.format(key))
                pass
        for key in self.cols:
            try:
                f.create_dataset(key,data = getattr(self,key))
            except RuntimeError:
                print('{} key already created! Skipping...'.format(key))
                pass
    def tss(self,ax=None,fig=None,**kargs):
        if ax is None:
            fig,ax= plt.subplots(figsize=(8,6))

        ax.plot(self.t,self.sk.max() - self.sk_t.min(axis=0),'.-')
        ax.set_xlabel('$t$',fontsize=20)
        ax.set_ylabel('$\\Delta s_{min}$',fontsize=20)
        ax.minorticks_on()
        return fig,ax
    def semilogx(self,q,**kargs):
        return self.plot(q, logx=True, logy=False, **kargs)
    def semilogy(self,q,**kargs):
        return self.plot(q, logx=False, logy=True, **kargs)
    def loglog(self,q,**kargs):
        return self.plot(q, logx=True, logy=True, **kargs)

    def plot(self,q,i=-1,ax=None,fig=None,shift=0,norm=1,logx=False, logy=False,average=False, savefig=None,window=100,**kargs):
        first = ax is None
        if first:
            fig,ax = plt.subplots(figsize=(8,6))

        try:
            ylbl = q
            q = getattr(self,q)
        except AttributeError:
            print('{} not an attribute!'.format(q))
            return
        except TypeError:
            ylbl=''
            pass

        if average:
            avg,st = window_avg((q-shift)/norm,window)
            line, = ax.plot(self.z,avg,**kargs)
            ax.fill_between(self.z, avg-st, avg+st, facecolor=line._color,alpha=.1)
        else:
            print(i,type(i))
            if i=='all':
                i  = range(0,self.nt)
                print(i,'all')
            if type(i) == int:
                i = [i]
                print('List', i)

            for ii in i:
                print(ii)
                ax.plot(self.z,(q[:,ii]-shift)/norm,**kargs)

        if first:
            ax.set_ylabel(ylbl,fontsize=20)
            ax.set_xlabel('$z$',fontsize=20)
            ax.minorticks_on()
            if logx:
                ax.set_xscale('log')
            if logy:
                ax.set_yscale('log')
        if savefig is not None:
            fig.savefig(savefig)
        return fig,ax
    def summary2(self,i=-1,axes=None,fig=None,savefig=None,**kargs):
        if axes is None:
            fig,axes = plt.subplots(2,2,sharex=True,figsize=(12,8))

        if i == 'all':
            i = range(self.nt)

        try:
            for ii in i:
                axes[0,0].plot(self.z, self.s[:,ii],**kargs)
                axes[1,0].plot(self.z, self.dsdz[:,ii],**kargs)
                axes[0,1].plot(self.z, self.vz2[:,ii],**kargs)
                axes[1,1].plot(self.z, self.Fc[:,ii],**kargs)
            ls=kargs.pop('ls','--')
            axes[1,1].set_prop_cycle(None)
            axes[0,1].set_prop_cycle(None)
            for ii in i:
                axes[0,1].plot(self.z, self.vx2[:,ii],ls='--',**kargs)
                axes[1,1].plot(self.z, self.Fk[:,ii],ls='--', **kargs)
            axes[1,1].set_prop_cycle(None)
            axes[0,1].set_prop_cycle(None)
            for ii in i:
                axes[0,1].plot(self.z, self.vy2[:,ii],ls=':',**kargs)
                axes[1,1].plot(self.z, self.Ftot[:,ii],ls=':',**kargs)
            axes[1,1].set_prop_cycle(None)
            for ii in i:
                axes[1,1].plot(self.z, self.Fe[:,ii],ls='-.',**kargs)
        except TypeError:
            axes[0,0].plot(self.z, self.s[:,i],**kargs)
            axes[1,0].plot(self.z, self.dsdz[:,i],**kargs)
            axes[0,1].plot(self.z, self.vz2[:,i],**kargs)
            axes[1,1].plot(self.z, self.Fc[:,i],**kargs)
            _=kargs.pop('ls','--')
            axes[1,1].set_prop_cycle(None)
            axes[0,1].set_prop_cycle(None)
            axes[1,1].plot(self.z, self.Fk[:,i],ls='--', **kargs)
            axes[0,1].plot(self.z, self.vx2[:,i],ls='--',**kargs)
            axes[1,1].set_prop_cycle(None)
            axes[0,1].set_prop_cycle(None)
            axes[1,1].plot(self.z, self.Ftot[:,i],ls=':',**kargs)
            axes[0,1].plot(self.z, self.vy2[:,i],ls=':',**kargs)
            axes[1,1].set_prop_cycle(None)
            axes[1,1].plot(self.z, self.Fe[:,i],ls='-.',**kargs)

        axes[0,0].set_ylabel('$\\langle s \\rangle$',fontsize=25)
        axes[1,0].set_ylabel('$\\langle ds/dz \\rangle$',fontsize=25)
        axes[0,1].set_ylabel('$\\langle v^2 \\rangle$',fontsize=25)
        axes[1,1].set_ylabel('$\\langle v_z T \\rangle$',fontsize=25)
        for ax in axes.flatten():
            ax.minorticks_on()
            ax.tick_params(labelsize=15)
            #ax.ticklabel_format(style='sci',axis='y',scilimits=(-3,3))
        for ax in axes[-1,:]:
            ax.set_xlabel('$z$',fontsize=25)
        fig.tight_layout()
        if savefig is not None:
            fig.savefig(savefig)
        return fig,axes
    def summary(self,axes=None,fig=None,savefig=None,**kargs):
        if axes is None:
            fig,axes = plt.subplots(2,2,sharex=True,figsize=(12,8))

        self.plot(self.s,fig=fig,ax=axes[0,0],**kargs)
        line = axes[0,0].lines[-1]
        ylim = axes[0,0].get_ylim()
        axes[0,0].plot(self.z,self.sk,ls='--',c=line._color)
        axes[0,0].set_ylim(ylim)
        self.plot(self.dsdz,fig=fig,ax=axes[1,0],**kargs)
        self.plot(self.vz2,fig=fig,ax=axes[0,1],**kargs)
        self.plot(self.Fc,fig=fig,ax=axes[1,1],**kargs)
        lv = axes[0,1].lines[-1]
        lf = axes[1,1].lines[-1]

        _=kargs.pop('ls','--')
        self.plot(self.vx2,fig=fig,ax=axes[0,1],ls='--',c=lv._color,**kargs)
        self.plot(self.Fk,fig=fig,ax=axes[1,1],ls='--',c=lf._color,**kargs)
        self.plot(self.vy2,fig=fig,ax=axes[0,1],ls=':',c=lv._color,**kargs)
        self.plot(self.Ftot,fig=fig,ax=axes[1,1],ls=':',c=lf._color,**kargs)
        self.plot(self.Fe,fig=fig,ax=axes[1,1],ls='-.',c=lf._color,**kargs)

        axes[0,0].set_ylabel('$\\langle s \\rangle$',fontsize=25)
        axes[1,0].set_ylabel('$\\langle ds/dz \\rangle$',fontsize=25)
        axes[0,1].set_ylabel('$\\langle v^2 \\rangle$',fontsize=25)
        axes[1,1].set_ylabel('$\\langle v_z T \\rangle$',fontsize=25)
        for ax in axes.flatten():
            ax.minorticks_on()
            ax.tick_params(labelsize=15)
            #ax.ticklabel_format(style='sci',axis='y',scilimits=(-3,3))
        for ax in axes[-1,:]:
            ax.set_xlabel('$z$',fontsize=25)
        fig.tight_layout()
        if savefig is not None:
            fig.savefig(savefig)
        return fig,axes
    def __add__(self,fld):
        newfld = copy.deepcopy(self)
        newfld.vx += fld.vx
        newfld.vy += fld.vy
        newfld.vz += fld.vz
        newfld.vx2 += fld.vx2
        newfld.vy2 += fld.vy2
        newfld.vz2 += fld.vz2
        newfld.KE += fld.KE
        newfld.dTdz += fld.dTdz
        newfld.T += fld.T
        newfld.T2 += fld.T2
        newfld.Fk += fld.Fk
        newfld.Fc += fld.Fc
        newfld.Fe += fld.Fe
        newfld.kap += fld.kap
        newfld.vort += fld.vort
        newfld.area += fld.area
        newfld.vortx += fld.vortx
        newfld.vorty += fld.vorty
        newfld.visc3 += fld.visc3
        newfld.therm += fld.therm
        newfld.P += fld.P
        newfld.dzP += fld.dzP
        newfld.Pflux += fld.Pflux
        newfld.TdzP += fld.TdzP
        newfld.trip1 += fld.trip1
        newfld.trip2 += fld.trip2
        newfld.rey12 += fld.rey12
        newfld.rey13 += fld.rey13
        newfld.rey23 += fld.rey23
        newfld.vflux += fld.vflux
        newfld.keflux += fld.keflux
        newfld.tflux += fld.tflux
        newfld.thermflux += fld.thermflux
        newfld.Fcflux1 += fld.Fcflux1
        newfld.Fcflux2 += fld.Fcflux2
        newfld.visc1 += fld.visc1
        newfld.visc2 += fld.visc2



        newfld.thermflux *= newfld.kap[:,-1][:,np.newaxis]
        newfld.Fb = -newfld.vz*newfld.delad*newfld.z[:,np.newaxis]
        newfld.Nu = newfld.Fc/newfld.Fk

        newfld.Ftot = newfld.Fk + newfld.Fc
        newfld.vrms = np.sqrt(newfld.vz2)

        newfld.s = newfld.T - newfld.delad*newfld.z[:,np.newaxis]
        newfld.dsdz = newfld.dTdz - newfld.delad
        newfld.sk = cumtrapz(-1./newfld.kap[:,-1]-newfld.delad,x=newfld.z,initial=0)
        tmp = cumtrapz((-1./newfld.kap[:,-1]-newfld.delad)[::-1],x=newfld.z[::-1],initial=0)[::-1]
        newfld.sk_t = newfld.s[-1,:] + tmp[:,np.newaxis]
        return newfld

def multiple_avg_plots(dirs, q, ylbl, lbls, ax=None,fig=None,legend=True,savefig=None, **kargs):
    if ax is None:
        fig,ax = plt.subplots(figsize=(8,6))
    for d,lbl in zip(dirs, lbls):
        Zavg(directory=d).plot(q,ax=ax,fig=fig,label=lbl, **kargs)

    ax.minorticks_on()
    ax.tick_params(labelsize=20)
    ax.set_xlabel('$z$',fontsize=25)
    ax.set_ylabel(ylbl,fontsize=25)
    fig.tight_layout()
    if legend:
        ax.legend(loc='best')

    if savefig is not None:
        fig.savefig(savefig)
    return fig,ax

def multiple_time_plots(dirs, q, ylbl, lbls, ax=None,fig=None, region=False, logx=False, logy=False, legend=True, savefig=None,**kargs):
    if ax is None:
        fig,ax = plt.subplots(figsize=(8,6))
    for d,lbl in zip(dirs, lbls):
        if region:
            Tavg(directory=d,fname='fort.52').plot(q,ax=ax,fig=fig,label=lbl, **kargs)
        else:
            Tavg(directory=d).plot(q,ax=ax,fig=fig,label=lbl,**kargs)

    if logy:
        ax.set_yscale('log')
    if logx:
        ax.set_xscale('log')
    ax.minorticks_on()
    ax.tick_params(labelsize=20)
    ax.set_xlabel('$t$',fontsize=25)
    ax.set_ylabel(ylbl,fontsize=25)
    fig.tight_layout()
    if legend:
        ax.legend(loc='best')
    if savefig is not None:
        fig.savefig(savefig)

    return fig,ax


def window_avg(y, weights=None,navg=100):
    try:
        return y[:,-navg:].mean(axis=1), y[:,-navg:].std(axis=1)
    except IndexError:
        return y[-navg:].mean(), y[-navg:].std()


def window_avg2(y, window=100):
    n = len(y)

    avg  = np.zeros(y.shape)
    st = np.zeros(y.shape)

    h = window//2

    avg[:h] = y[:window].mean()
    st[:h] = y[:window].std()
    for i in range(h,n-h):
        avg[i] = y[(i-h):(i+h)].mean()
        st[i] = y[(i-h):(i+h)].std()
    for i in range(n-h,n):
        avg[i] = y[-window:].mean()
        st[i] = y[-window:].std()
    return avg,st

def get_dir_values(d):
    """ ds30/t2_22/k003_nu003/"""
    res = d.split('/')
    ds = res[0]
    t = res[1]
    k = res[2]
    dsdz = ds.split('ds')[1]
    k1, delad = t.split('t')[-1].split('_')
    k = k.split('_')
    if k[0][:3] == 'kap':
        kap = '.{}'.format(k[0].split('kap')[-1])
    else:
        kap = '.{}'.format(k[0].split('k')[-1])
    nu = '.{}'.format(k[1].split('nu')[-1])

    return  float(dsdz), 1./float(k1), float(delad), float(kap), float(nu)

def get_dir_names(ds):
    k1 = [2,40]

    if ds in [20,40]:
        kstr = 'k'
    else:
        kstr = 'kap'

    dirs = []
    vals = []
    for i  in k1:
        vals.append((ds,i,i+ds,.003,.003))
        vals.append((ds,i,i+ds,.003,.01))
        vals.append((ds,i,i+ds,.01,.003))
        dirs.append('ds{:d}/t{:d}_{:d}/{}003_nu003/'.format(ds,i,i+ds,kstr))
        dirs.append('ds{:d}/t{:d}_{:d}/{}003_nu01/'.format(ds,i,i+ds,kstr))
        dirs.append('ds{:d}/t{:d}_{:d}/{}01_nu003/'.format(ds,i,i+ds,kstr))

    lbls = ['$\\nabla = {:d},\\kappa_1 = {:.2f}, \\nabla_{{ad}}={:d},\\kappa_2={:.3f},\\nu ={:.3f}$'.format(v[0],1./v[1],v[2],v[3],v[4]) for v in vals]

    return dirs,vals,lbls

def create_h5_dir(ds, f):
    for d in ds:
        print('Creating {}'.format(d))
        vavg = f.create_group(d + 'Volavg')
        vavg1 = f.create_group(d + 'Volavg1')
        zavg1 = f.create_group(d + 'Zavg')

        Tavg(directory=d).dump(vavg)
        Tavg(fname='fort.52',directory=d).dump(vavg1)
        Zavg(directory=d).dump(zavg1)
def create_h5_file(fname,dirs):
    import h5py
    with h5py.File(fname,'w') as f:
        create_h5_dir2(dirs,f)

def create_h5_dir2(ds, f):
    for d in ds:
        if d[-1] != '/':
            d += '/'
        print('Creating {}'.format(d))
        vavg_down = f.create_group(d + 'Volavg_down')
        vavg_up = f.create_group(d + 'Volavg_Up')
        vavg1_down = f.create_group(d + 'Volavg1_down')
        vavg1_up = f.create_group(d + 'Volavg1_up')
        zavg_down = f.create_group(d + 'Zavg_down')
        zavg_up = f.create_group(d + 'Zavg_up')

        Zavg(directory=d,fname='zavg_down.dat').dump(zavg_down)
        Zavg(directory=d,fname='zavg_up.dat').dump(zavg_up)
        Tavg(directory=d,fname='volvg_down.dat').dump(vavg_down)
        Tavg(directory=d,fname='volvg_up.dat').dump(vavg_up)
        Tavg(directory=d,fname='volavgR_down.dat').dump(vavg1_down)
        Tavg(directory=d,fname='volavgR_up.dat').dump(vavg1_up)








def load_mesh(fname='RB.rea'):
    with open(fname,'r') as f:
        lines = f.readlines()

    for i,line in enumerate(lines):
        if 'MESH DATA' in line:
            start= i
        if 'CURVED SIDE' in line:
            end = i

    lines = lines[start+2:end]

    nelements = len(lines)//7



    mesh = []
    for i in range(nelements):
        dat =np.array( [x.strip().split() for x in lines[i*7:(i+1)*7][1:]]).astype(float)
        x = np.hstack((dat[0,:],dat[3,:]))
        y = np.hstack((dat[1,:],dat[4,:]))
        z = np.hstack((dat[2,:],dat[5,:]))
        pts = np.vstack((x,y,z)).T
        mesh.append(pts)
    return mesh

def find_point(point,mesh):
    for i,pts in enumerate(mesh):
        p1 = pts[0]
        p2 = pts[1]
        p3 = pts[3]
        p4 = pts[4]

        v1 = point - p1
        v2 = p2-p1
        v3 = p3-p1
        v4 = p4-p1

        res = (np.dot(v1,v2) <= np.dot(v2,v2))&(np.dot(v1,v2) >= 0)
        res &= (np.dot(v1,v3) <= np.dot(v3,v3))&(np.dot(v1,v3) >= 0)
        res &= (np.dot(v1,v4) <= np.dot(v4,v4))&(np.dot(v1,v4) >= 0)
        if res:
            print('Found ({:.3f},{:.3f},{:.3f}) '.format(*point) + 'in element {:d}'.format(i))
            print(pts)
            return i
    print('Point ({:.3f},{:.3f},{:.3f}) not in mesh!'.format(*point))
    return i












