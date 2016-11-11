import numpy as np
import matplotlib.pyplot as plt



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
