import numpy as np
import matplotlib.pyplot as plt
import copy

class Disk():
    def __init__(self,fname='outputs/res.dat',p=0):
        fname = fname +'.{:d}'.format(p)
        dat = np.fromfile(fname)
        self.n,self.nm,self.mi,self.mf = dat[:4].astype(int)
        dat=dat[4:]

        self.TL = dat[:self.nm]
        dat=dat[self.nm:]
        self.TR = dat[:self.nm]
        dat=dat[self.nm:]
        self.r = dat[:self.n]
        dat=dat[self.n:]
        self.lamex = dat[:self.n*self.nm].reshape(self.nm,self.n).T
        dat=dat[self.n*self.nm:]
        self.lamdep = dat[:self.n*self.nm].reshape(self.nm,self.n).T
        dat=dat[self.n*self.nm:]
        self.drfw = dat[:self.n*self.nm].reshape(self.nm,self.n).T
        dat=dat[self.n*self.nm:]
        self.fw = dat[:self.n*self.nm].reshape(self.nm,self.n).T
        dat=dat[self.n*self.nm:]

        self.u = np.zeros((self.n,self.nm),dtype='complex')
        self.v = np.zeros((self.n,self.nm),dtype='complex')
        self.s = np.zeros((self.n,self.nm),dtype='complex')
        for indx,i in enumerate(np.arange(self.nm)+self.mi-1):
            dat = np.loadtxt('outputs/sol{:d}.dat.{:d}'.format(i+1,p))
            self.u[:,indx] = dat[:,0] + 1j*dat[:,1]
            self.v[:,indx] = dat[:,2] + 1j*dat[:,3]
            self.s[:,indx] = dat[:,4] + 1j*dat[:,5]


        self.dlr = np.diff(np.log(self.r))[0]
    def torque(self,m,logx=True,xlims=None,ax=None,integ=False):
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)

        i = m-1

        lamex = self.lamex[:,i].copy()
        lamdep = self.lamdep[:,i].copy()
        if integ:
            fw  = 2*np.pi*(self.fw[:,i])
            ind = self.r >= 1
            ilamex = np.zeros(lamex.shape)
            ilamex[ind] = (lamex*2*np.pi*self.r*self.r*self.dlr)[ind].cumsum()
            ilamdep= np.zeros(lamdep.shape)
            ilamdep[ind] = (lamdep*2*np.pi*self.r*self.r*self.dlr)[ind].cumsum()
            ind = self.r <= 1
            ilamex[ind] = -(lamex*2*np.pi*self.r*self.r*self.dlr)[ind][::-1].cumsum()[::-1]
            ilamdep[ind] = -(lamdep*2*np.pi*self.r*self.r*self.dlr)[ind][::-1].cumsum()[::-1]
            lamdep = ilamdep
            lamex = ilamex
        else:
            fw = self.drfw[:,i].copy()

        ax.plot(self.r,lamex,'-k')
        ax.plot(self.r,fw,'-r')
        ax.plot(self.r,lamdep,'-b')
        if logx:
            ax.set_xscale('log')
        if xlims is not None:
            ax.set_xlim(xlims)
    def plot(self,m,logx=True):
        fig,axes = plt.subplots(1,3,figsize=(15,5))
        i = m-1
        axes[0].plot(self.r,self.u.real[:,i],self.r,self.u.imag[:,i])
        axes[1].plot(self.r,self.v.real[:,i],self.r,self.v.imag[:,i])
        axes[2].plot(self.r,self.s.real[:,i],self.r,self.s.imag[:,i])

        if logx:
            for ax in axes:
                ax.set_xscale('log')
    def __add__(self,fld):
        out = copy.deepcopy(self)
        out.TL = np.hstack((self.TL,fld.TL))
        out.TR = np.hstack((self.TR,fld.TR))
        out.lamex = np.hstack((self.lamex,fld.lamex))
        out.lamdep = np.hstack((self.lamdep,fld.lamdep))
        out.drfw = np.hstack((self.drfw,fld.drfw))
        out.fw = np.hstack((self.fw,fld.fw))
        out.u = np.hstack((self.u,fld.u))
        out.v = np.hstack((self.v,fld.v))
        out.s = np.hstack((self.s,fld.s))

        out.mi = min(self.mi,fld.mi)
        out.mf = max(self.mf,fld.mf)
        out.nm = self.nm + fld.nm

        return out

class Mode():
    def __init__(self,fname='outputs/res.dat'):
        dat = np.fromfile(fname)
        n = int(dat[0])
        self.n = n
        dat = dat[1:]
        self.r = dat[:n]
        dat = dat[n:]
        self.u = dat[:n] + 1j*dat[n:2*n]
        dat = dat[2*n:]
        self.v = dat[:n] + 1j*dat[n:2*n]
        dat = dat[2*n:]
        self.s = dat[:n] + 1j*dat[n:2*n]
        dat = dat[2*n:]
        self.lamex = dat[:n]
        dat = dat[n:]
        self.lamdep = dat[:n]
        dat = dat[n:]
        self.drfw = dat[:n]
        dat = dat[n:]
        self.fw = dat[:n]
        self.phase = np.arctan2(self.s.imag,self.s.real)
        self.dlr = np.diff(np.log(self.r))[0]

        self.ilam =np.zeros(self.lamex.shape)
        self.ildep =np.zeros(self.lamdep.shape)
        ind = self.r>=1
        self.ilam[ind] = (self.lamex*self.r**2*self.dlr)[ind].cumsum()
        self.ilam[ind] -= self.ilam[ind][0]
        self.ildep[ind] = (self.lamdep*self.r**2*self.dlr)[ind].cumsum()
        self.ildep[ind] -= self.ildep[ind][0]
        ind = self.r<=1
        self.ilam[ind] = -(self.lamex*self.r**2*self.dlr)[ind][::-1].cumsum()[::-1]
        self.ilam[ind] -= self.ilam[ind][-1]
        self.ildep[ind] = -(self.lamdep*self.r**2*self.dlr)[ind][::-1].cumsum()[::-1]
        self.ildep[ind] -= self.ildep[ind][-1]

        #self.drfw = np.gradient(self.fw)/(self.r*self.dlr)

    def summary(self,logx=True):
        fig,axes = plt.subplots(1,3,figsize=(15,5))
        axes[0].plot(self.r,self.u.real,self.r,self.u.imag)
        axes[1].plot(self.r,self.v.real,self.r,self.v.imag)
        axes[2].plot(self.r,self.s.real,self.r,self.s.imag)

        if logx:
            for ax in axes:
                ax.set_xscale('log')

    def torque(self,logx=True):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(self.r,self.ilam,'-k')
        ax.plot(self.r,self.fw,'-r')
        ax.plot(self.r,self.ildep,'-b')
        if logx:
            ax.set_xscale('log')
    def dr_torque(self,logx=True):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(self.r,self.lamex,'-k')
        ax.plot(self.r,self.drfw,'-r')
        ax.plot(self.r,self.lamdep,'-b')
        if logx:
            ax.set_xscale('log')

    def phase_plot(self,logx=True):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        phase = np.arctan2(self.s.imag,self.s.real)
        phase = np.unwrap(phase,np.pi)
        dpdr = np.zeros(phase.shape)
        dlr = np.diff(np.log(self.r))[0]
        dpdr[1:-1] = (phase[2:]-phase[:-2])/(2*dlr)
        dpdr[0] = (phase[1]-phase[0])/dlr
        dpdr[-1] = (phase[-1] - phase[-2])/dlr
        dpdr /= self.r
        ax.plot(self.r,dpdr)

        if logx:
            ax.set_xscale('log')



