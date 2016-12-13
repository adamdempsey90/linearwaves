import numpy as np
import matplotlib.pyplot as plt
import copy

class Disk():
    def __init__(self,fname='outputs/res.dat',p=0):
        dfname = fname + '.disk'
        fname = fname +'.{:d}'.format(p)
        dat = np.fromfile(fname)
        self.n,self.nm,self.mi,self.mf = dat[:4].astype(int)
        print(p,self.n,self.nm,self.mi,self.mf)
        dat=dat[4:]
        self.mvals = np.arange(self.mi,self.mf+1)

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
        self.pot = dat[:self.n*self.nm].reshape(self.nm,self.n).T
        dat=dat[self.n*self.nm:]

        self.u = np.zeros((self.n,self.nm),dtype='complex')
        self.v = np.zeros((self.n,self.nm),dtype='complex')
        self.s = np.zeros((self.n,self.nm),dtype='complex')
        for indx,i in enumerate(np.arange(self.nm)+self.mi):
            dat = np.loadtxt('outputs/sol{:d}.dat.{:d}'.format(i,p))
            self.u[:,indx] = dat[:,0] + 1j*dat[:,1]
            self.v[:,indx] = dat[:,2] + 1j*dat[:,3]
            self.s[:,indx] = dat[:,4] + 1j*dat[:,5]

        self.dlr = np.diff(np.log(self.r))[0]

        dat = np.fromfile(dfname)[1:]
        self.dbar = dat[:self.n]
        self.dlsdlr = dat[self.n:2*self.n]
        self.d2lsdlr = dat[self.n*2:3*self.n]

        self.lamex_tot = self.lamex[:,1:].sum(axis=1)
        self.lamdep_tot = self.lamdep[:,1:].sum(axis=1)
        self.drfw_tot = self.drfw[:,1:].sum(axis=1)
        self.fw_tot = self.fw[:,1:].sum(axis=1)
        self.pot_tot = self.pot[:,1:].sum(axis=1)
        self.drfd = self.drfw[:,0].copy()
        self.mdot_d = self.fw[:,0].copy()

        indR = self.r>=1
        indL = self.r<=1
        self.rR = self.r[indR]
        self.rL = self.r[indL]
        self.ilamexR = (self.lamex*2*np.pi*self.r[:,np.newaxis]**2 * self.dlr)[indR,:].cumsum(axis=0)
        self.ilamexL = -(self.lamex*2*np.pi*self.r[:,np.newaxis]**2 * self.dlr)[indL,:][::-1].cumsum(axis=0)[::-1]


    def torque_contours(self,ax=None,scaleH=True,norm=1,cmap='bwr',ylims=None,xlims=None):
        if ax is None:
            fig=plt.figure()
            ax=fig.add_subplot(111)

        rvals = np.hstack((self.rL,self.rR))
        if xlims is not None:
            indx = (rvals>=xlims[0])&(rvals<=xlims[1])
        else:
            indx = np.ones(rvals.shape).astype(bool)
        if ylims is not None:
            indy = (self.mvals>=ylims[0])&(self.mvals<=ylims[1])
        else:
            indy = np.ones(self.mvals.shape).astype(bool)
        yy,mm = np.meshgrid(rvals[indx],self.mvals[indy],indexing='ij')
        if scaleH:
            yy = (yy-1)/.05
        TT = np.vstack((self.ilamexL,self.ilamexR))
        TT = TT[indx,:][:,indy]
        img=ax.contour(yy,mm,TT/norm,100,cmap=cmap)
        return img

    def gt80_torque_density(self,r,mu=-.5,hp=.05,d=1.3):
        res = np.zeros(self.dbar.shape)
        ind = abs(r-1)>=d*hp
        res[ind] = np.sign(r[ind]-1)*2*np.pi*.4/abs(r[ind]-1)**4 * self.dbar[ind]
        return res
    def total_torque(self,logx=False,xlims=None,axes=None,cumulative=True,norm=(1./.05**2)):
        if axes is None:
            if cumulative:
                fig,axes = plt.subplots(1,2,figsize=(10,5))
            else:
                fig = plt.figure()
                axes = fig.add_subplot(111)

        if cumulative:
            axes[0].plot(self.mvals,self.TR/norm,'ok',label='Outer')
            axes[0].plot(self.mvals,self.TL/norm,'sr',label='Inner')
            axes[1].plot(self.mvals,self.TR.cumsum()/self.TR.sum(),'--ok')
            axes[1].plot(self.mvals,self.TL.cumsum()/self.TL.sum(),'-sr')
            for ax in axes:
                if logx:
                    ax.set_xscale('log')
                    ax.set_xlabel('m',fontsize=15)
            axes[0].legend(loc='upper right')
        else:
            axes.plot(self.mvals,self.TR/norm,'ok',label='Outer')
            axes.plot(self.mvals,self.TL/norm,'sr',label='Inner')
            if logx:
                axes.set_xscale('log')
                axes.set_xlabel('m',fontsize=15)
            axes.legend(loc='upper right')

    def torque(self,m,logx=True,xlims=None,ax=None,integ=False,from_inner=False,tot=False,**kargs):
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)

        i = m-1

        if tot:
            lamex = self.lamex_tot.copy()
            lamdep = self.lamdep_tot.copy()
        else:
            lamex = self.lamex[:,i].copy()
            lamdep = self.lamdep[:,i].copy()
        if integ:
            if tot:
                fw  = 2*np.pi*(self.fw_tot)
            else:
                fw  = 2*np.pi*(self.fw[:,i])

            if from_inner:
                lamex = 2*np.pi*(lamex*self.r**2*self.dlr).cumsum()
                lamdep = 2*np.pi*(lamdep*self.r**2*self.dlr).cumsum()
            else:
                #fr = self.fw[:,i][self.r>=1][0]
                #fl = self.fw[:,i][self.r<=1][-1]
                #rr = self.r[self.r>=1][0]
                #rl = self.r[self.r<=1][-1]
                #fp = fl + (fr-fl)/(np.log(rr/rl)) *(np.log(1./rl))
                #fw -= 2*np.pi*fp

                ind = self.r >= 1
                ilamex = np.zeros(lamex.shape)
                ilamex[ind] = (lamex*2*np.pi*self.r*self.r*self.dlr)[ind].cumsum()
                ilamex[ind] -= ilamex[ind][0]
                ilamdep= np.zeros(lamdep.shape)
                ilamdep[ind] = (lamdep*2*np.pi*self.r*self.r*self.dlr)[ind].cumsum()
                ilamdep[ind] -= ilamdep[ind][0]
                ind = self.r <= 1
                ilamex[ind] = -(lamex*2*np.pi*self.r*self.r*self.dlr)[ind][::-1].cumsum()[::-1]
                ilamdep[ind] = -(lamdep*2*np.pi*self.r*self.r*self.dlr)[ind][::-1].cumsum()[::-1]
                ilamex[ind] -= ilamex[ind][-1]
                ilamdep[ind] -= ilamdep[ind][-1]
                lamdep = ilamdep
                lamex = ilamex
        else:
            if tot:
                fw = self.drfw_tot.copy()
            else:
                fw = self.drfw[:,i].copy()

        ax.plot(self.r,lamex,c='k',**kargs)
        ax.plot(self.r,fw,c='r',**kargs)
        ax.plot(self.r,lamdep,c='b',**kargs)
        if logx:
            ax.set_xscale('log')
        if xlims is not None:
            ax.set_xlim(xlims)
            return lamex,lamdep,fw
    def plot(self,m,logx=True,axes=None,fig=None,xlims=None,ylims=None):
        if axes is None:
            fig,axes = plt.subplots(1,3,figsize=(15,5))
        if m == 0:
            axes[0].plot(self.r,self.dbar,'-k')
            axes[1].plot(self.r,self.dlsdlr,'-k')
            axes[2].plot(self.r,self.d2lsdlr,'-k')
        else:
            i = m-1
            axes[0].plot(self.r,self.u.real[:,i],self.r,self.u.imag[:,i])
            axes[1].plot(self.r,self.v.real[:,i],self.r,self.v.imag[:,i])
            axes[2].plot(self.r,self.s.real[:,i],self.r,self.s.imag[:,i])

        for ax in axes:
            if logx:
                ax.set_xscale('log')
            if xlims is not None:
                ax.set_xlim(xlims)
            if ylims is not None:
                ax.set_ylim(ylims)
    def fullplot(self,q='s',rlims=None,cartesian=True,ax=None,full=False):

        y = self.r.copy()
        phi = np.linspace(-np.pi,np.pi,2000)
        if rlims is not None:
            ind =(y>=rlims[0])&(y<=rlims[1])

        y = y[ind]
        yy,xx = np.meshgrid(y,phi,indexing='ij')

        if cartesian:
            x = np.cos(xx)*yy
            y = np.sin(xx)*yy
        else:
            x = xx
            y = yy



        if q[0] == 's':
            dat = self.s[ind,:]
            dat0 = self.dbar[ind]
        elif q[0] == 'u':
            dat = self.u[ind,:]
            dat0 = np.zeros(y.shape)
        elif q[0] == 'v':
            dat = self.v[ind,:]
            dat0 = y**(-.5)
        else:
            print('{} not a valid choice!'.format(q))
            return

        res = np.zeros(xx.shape)

        for i in range(len(y)):
            res[i,:] = sum([np.real( dat[i,j]*np.exp(1j*m*phi)) for j,m in enumerate(self.mvals)])

        if full:
            res += dat0[:,np.newaxis]


        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)

        ax.pcolormesh(x,y,dat,cmap='viridis')
        plt.colorbar()


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



