import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import copy

class Disk():
    def __init__(self,fname='outputs/res.dat',p=0,second=False):
        dfname = fname + '.disk'
        fname_d = fname +'.{:d}'.format(p)
        dat = np.fromfile(fname_d)
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
        print(fname)
        if second:
            dat = np.fromfile(fname+'.2',dtype='complex').reshape(self.nm,3,self.n).T
        dat = np.fromfile(fname,dtype='complex').reshape(self.nm,3,self.n).T

        self.u = dat[:,0,:]
        self.v = dat[:,1,:]
        self.s = dat[:,2,:]
        if second:
            dat = np.fromfile(fname+'.2',dtype='complex').reshape(self.nm,3,self.n).T
        self.u += dat[:,0,:]
        self.v += dat[:,1,:]
        self.s += dat[:,2,:]
        #for indx,i in enumerate(np.arange(self.nm)+self.mi):
        #    dat = np.loadtxt('outputs/sol{:d}.dat.{:d}'.format(i,p))
        #    self.u[:,indx] = dat[:,0] + 1j*dat[:,1]
        #    self.v[:,indx] = dat[:,2] + 1j*dat[:,3]
        #    self.s[:,indx] = dat[:,4] + 1j*dat[:,5]

        self.dlr = np.diff(np.log(self.r))[0]

        dat = np.fromfile(dfname)[1:]
        self.dbar = dat[:self.n]; dat=dat[self.n:]
        self.dlsdlr = dat[:self.n]; dat=dat[self.n:]
        self.d2lsdlr = dat[:self.n]; dat=dat[self.n:]
        self.omega = dat[:self.n]; dat=dat[self.n:]
        self.dlomdlr = dat[:self.n]; dat=dat[self.n:]
        self.kappa2 = dat[:self.n]; dat=dat[self.n:]

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
        self.ilam = np.vstack((self.ilamexL,self.ilamexR))


    def Dfunc(self,m):
        return self.kappa2 - m*(self.omega - 1)**2
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
    def total_torque(self,logx=False,xlims=None,axes=None,cumulative=True,fig=None,norm=(1./.05**2)):
        if axes is None:
            if cumulative:
                fig,axes = plt.subplots(1,2,figsize=(10,5))
            else:
                fig = plt.figure()
                axes = fig.add_subplot(111)

        if cumulative:
            axes[0].plot(self.mvals,self.TR/norm,'ok',label='Outer')
            axes[0].plot(self.mvals,self.TL/norm,'sr',label='Inner')
            axes[0].plot(self.mvals,self.TR/norm-self.TL/norm,'db',label='$\\Delta T$')
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
            axes.plot(self.mvals,self.TR/norm-self.TL/norm,'db',label='$\\Delta T$')
            if logx:
                axes.set_xscale('log')
                axes.set_xlabel('m',fontsize=15)
            axes.legend(loc='upper right')
        return fig,axes

    def itorque(self,m,logx=True,xlims=None,ax=None,plot_dep=False,tot=False,**kargs):
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

        ax.plot(self.r,ilamex,ls='-',**kargs)
        if plot_dep:
            ax.plot(self.r,ilamdep,ls='--',**kargs)
        ax.set_xlabel('$r/r_p$',fontsize=20)

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
                ilamdep[ind] -= ilamex[ind][0]
                ind = self.r <= 1
                ilamex[ind] = -(lamex*2*np.pi*self.r*self.r*self.dlr)[ind][::-1].cumsum()[::-1]
                ilamdep[ind] = -(lamdep*2*np.pi*self.r*self.r*self.dlr)[ind][::-1].cumsum()[::-1]
                ilamex[ind] -= ilamex[ind][-1]
                ilamdep[ind] -= ilamex[ind][-1]
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

        if m == 0:
            if axes is None:
                fig,axes = plt.subplots(2,3,figsize=(15,5))
            axes[0,0].plot(self.r,self.dbar,'-k')
            axes[0,1].plot(self.r,self.dlsdlr,'-k')
            axes[0,2].plot(self.r,self.d2lsdlr,'-k')
            axes[1,0].plot(self.r,self.omega,'-k')
            axes[1,1].plot(self.r,self.dlomdlr,'-k')
            axes[1,2].plot(self.r,self.kappa2,'-k')
            axes[0,0].text(.75,.8,'$\\Sigma$',fontsize=20,transform=axes[0,0].transAxes)
            axes[0,1].text(.75,.8,'$\\frac{d \\ln \\Sigma}{ d\\ln r}$',fontsize=20,transform=axes[0,1].transAxes)
            axes[0,2].text(.75,.8,'$\\frac{d^2 \\ln \\Sigma}{ d \\ln r^2}$',fontsize=20,transform=axes[0,2].transAxes)
            axes[1,0].text(.75,.8,'$\\Omega$',fontsize=20,transform=axes[1,0].transAxes)
            axes[1,1].text(.75,.8,'$\\frac{d \\ln \\Omega}{d \\ln r} $',fontsize=20,transform=axes[1,1].transAxes)
            axes[1,2].text(.75,.8,'$\\kappa^2$',fontsize=20,transform=axes[1,2].transAxes)
        else:
            if axes is None:
                fig,axes = plt.subplots(1,4,figsize=(20,5))
            i = m-1
            axes[0].plot(self.r,self.u.real[:,i],self.r,self.u.imag[:,i])
            axes[1].plot(self.r,self.v.real[:,i],self.r,self.v.imag[:,i])
            axes[2].plot(self.r,self.s.real[:,i],self.r,self.s.imag[:,i])
            self.torque(m,ax=axes[3],integ=True)

        for ax in axes.flatten():
            if logx:
                ax.set_xscale('log')
            if xlims is not None:
                ax.set_xlim(xlims)
            if ylims is not None:
                ax.set_ylim(ylims)
        return fig,axes
    def mod_pi(self,angle):
        from copy import copy
        p = copy(angle)
        try:
            while (p < -np.pi):
                p += 2*np.pi
            while (p > np.pi):
                p -= 2*np.pi
        except ValueError:
            for i in range(len(p)):
                while (p[i] < -np.pi):
                    p[i] += 2*np.pi
                while (p[i] > np.pi):
                    p[i] -= 2*np.pi

        return p

    def interp_2d(self,r_range,q='s',norm=1,h=.05,d=.5,mu=.5):
        from scipy.interpolate import RectBivariateSpline
        ind = (self.r>=r_range[0])&(self.r<=r_range[1])
        r = self.r[ind].copy()
        #dlr = np.diff(np.log(r))[0]

        if q[0] == 's':
            dat = self.s[ind,:]/norm
        elif q[0] == 'u':
            dat = self.u[ind,:]/norm
        elif q[0] == 'v':
            dat = self.v[ind,:]/norm
        else:
            print('{} not a valid choice!'.format(q))
            return
        #dr_dat = np.gradient(dat,axis=0)/(dlr*r[:,np.newaxis])

        nr = dat.shape[0]
        phi = np.linspace(-np.pi,np.pi,nr)
        nphi = len(phi)
        nm =len(self.mvals)

        dat_f = np.zeros((nr,nphi))
        #dr_dat_f = np.zeros((nr,nphi))
        dp_dat_f = np.zeros((nr,nphi))



        for j,m in enumerate(self.mvals):
            for i in range(nr):
                dat_f[i,:] += np.real(dat[i,j]*np.exp(1j*m*phi))/nm
                #dr_dat_f[i,:] += np.real(dr_dat[i,j]*np.exp(1j*m*phi))
                dp_dat_f[i,:] += np.real(1j*m*dat[i,j]*np.exp(1j*m*phi))/nm


        #S_dat =( 2**(5./4) * h**(3./2) / 3.) *(  (h*dr_dat_f + r[:,np.newaxis]**(.5) *(r[:,np.newaxis]**(-1.5)-1)*dp_dat_f)/(abs(r[:,np.newaxis]-1)**(1.5)*r[:,np.newaxis]**( .5*(5*d+mu)-11./4))) - (2.*h/3.)*dat_f*dp_dat_f
        g = 2**(.25) * h**(.5) * r**( .5*(1-mu-3*d))/np.sqrt(abs(r**(1.5)-1))

        func_chi = RectBivariateSpline(r,phi, dat_f*g[:,np.newaxis])
        func_s = RectBivariateSpline(r,phi, -(2*h/3)*dp_dat_f*dat_f*g[:,np.newaxis]**2)

        func = RectBivariateSpline(r,phi,dat_f)
        #func_r = RectBivariateSpline(r,phi,dr_dat_f)
        func_p = RectBivariateSpline(r,phi,dp_dat_f)
        return func,func_p,func_chi, func_s

    def get_ring(self,r,q='s',shift=0,norm=1,skiprange=None,skiplist=[],skip=False,skip_high=False):
        if skip:
            if skiprange is not None:
                skiplist = range(skiprange[0],skiprange[1]+1)
            if skip_high:
                mtop = 50
            else:
                mtop = 1e6
        else:
            mtop = 1e6
            skiplist=[]
        ind = np.argwhere(self.r>=r)[0][0]
        if q[0] == 's':
            dat = self.s[ind,:]
        elif q[0] == 'u':
            dat = self.u[ind,:]
        elif q[0] == 'v':
            dat = self.v[ind,:]
        else:
            print('{} not a valid choice!'.format(q))
            return

        ny = len(self.r)
        phi = np.linspace(-np.pi,np.pi,ny)
        nm = len(self.mvals)
        nm = max(1,nm)
        #res = sum([np.real( dat[j]*np.exp(1j*m*phi)) for j,m in enumerate(self.mvals) if m not in skiplist and m <= mtop])/nm

        res = np.zeros(phi.shape)
        for j,m in enumerate(self.mvals):
            if m not in skiplist and m <= mtop:
                res += np.real(dat[j]*np.exp(1j*m*phi))
        res /= nm

        power = np.real(dat*np.conj(dat))
        normp = sum(power)
        for m in skiplist:
            power[m] = 0
        if skip_high:
            power[50:] = 0

        if normp > 0:
            power /= normp

        if norm == 'max':
            res /= abs(res).max()
        else:
            res /= norm
        res -= shift


        return phi,res,power
    def plotslice_power(self,r,q='s',shift=0,norm=1,ax=None,fig=None,skiprange=None,skiplist=[],skip=False,skip_high=False,phicenter=False,**kargs):
        if ax is None:
            fig,ax = plt.subplots(figsize=(8,6))
        phi,res,power = self.get_ring(r,q=q,shift=shift,norm=norm,skiprange=skiprange,skiplist=skiplist,skip=skip,skip_high=skip_high)
        c = kargs.pop('c','k')
        marker = kargs.pop('marker','s')
        ax.plot(range(self.u.shape[1]),power,marker=marker,c=c,**kargs)
        ax.minorticks_on()
        ax.set_xlabel('$m$')
        ax.set_ylabel('$|A_m|^2$')
        ax.set_yscale('log')
        ax.set_xscale('log')
        ax.set_ylim((1e-4,ax.get_ylim()[1]))
        return fig,ax
    def plotslice(self,r,q='s',shift=0,norm=1,ax=None,fig=None,skiprange=None,skiplist=[],skip=False,skip_high=False,phicenter=False,**kargs):
        if ax is None:
            fig,ax = plt.subplots(figsize=(8,6))
        phi,res,_ = self.get_ring(r,q=q,shift=shift,norm=norm,skiprange=skiprange,skiplist=skiplist,skip=skip,skip_high=skip_high)
        c = kargs.pop('c','k')
        ls = kargs.pop('ls','-')
        if phicenter:
            phi -= 2./(3*.05) * (r**(1.5) - 1.5*np.log(r)-1)
        ax.plot(phi,res,ls=ls,c=c,**kargs)
        ax.minorticks_on()
        ax.set_xlabel('$\\phi$')
        ax.set_ylabel('$\\sigma(r={:.2f})$'.format(r))
        return fig,ax


    def fullplot(self,q='s',skip=1,rlims=None,cartesian=True,ax=None,fig=None,full=False,norm=1,shift=0,log=False,contour=False,skiplist = [],skiprange=None,skip_high=False,divbar=False,clrbar=True,**kargs):

        if skiprange is not None:
            skiplist = range(skiprange[0],skiprange[1]+1)
        if skip_high:
            mtop = 50
        else:
            mtop = 1e6

        r = self.r.copy()
        if rlims is not None:
            ind =(r>=rlims[0])&(r<=rlims[1])
        else:
            ind = np.ones(r.shape).astype(bool)

        r = r[ind][::skip]
        ny = len(r)
        phi = np.linspace(-np.pi,np.pi,ny)
        rr,pp = np.meshgrid(r,phi,indexing='ij')




        if q[0] == 's':
            dat = self.s[ind,:][::skip,:]
            dat0 = self.dbar[ind][::skip]
        elif q[0] == 'u':
            dat = self.u[ind,:][::skip,:]
        elif q[0] == 'v':
            dat = self.v[ind,:][::skip,:]
            dat0 = r**(-.5)
        else:
            print('{} not a valid choice!'.format(q))
            return
        res = np.zeros((ny,ny))
        nm = len(self.mvals)


        for i in range(len(r)):

            res[i,:] = sum([np.real( dat[i,j]*np.exp(1j*m*phi)) for j,m in enumerate(self.mvals) if m not in skiplist and m <= mtop])/nm

        if divbar:
            res /= dat0[:,np.newaxis]
        res = (res-shift)/norm


        if full:
            res += dat0[:,np.newaxis]

        vmin = kargs.pop('vmin',res.min())
        vmax = kargs.pop('vmax',res.max())
        cmap = kargs.pop('cmap','viridis')
        if log:
            norm = colors.LogNorm(vmin=vmin,vmax=vmax)
        else:
            norm = colors.Normalize(vmin=vmin,vmax=vmax)
        if ax is None:
            fig = plt.figure(figsize=(10,10))
            ax = fig.add_subplot(111,polar=cartesian)

        if cartesian:
            if contour:
                img = ax.contour(pp,rr,res,100,cmap=cmap)
            else:
                img = ax.pcolormesh(pp,rr,res,cmap=cmap,norm=norm)
        else:
            if contour:
                img = ax.contour(rr,pp,res,100,cmap=cmap)
            else:
                img = ax.pcolormesh(rr,pp,res,cmap=cmap,norm=norm)
        if clrbar:
            plt.colorbar(img,ax=ax)
        if cartesian:
            ax.set_yticklabels('')
            ax.set_xticklabels('')
            ax.set_xlabel('')
            ax.set_ylabel('')
        else:
            ax.minorticks_on()
            ax.set_xlabel('$r/a$')
            ax.set_ylabel('$\\phi-\\phi_p$')
        return fig,ax
    def burgers_t(self,r,mu=.5,d=.5,h=.05):
        from scipy.integrate import quad
        try :
            t = np.zeros(r.shape)
            for i,j in enumerate(r):
                t[i] =abs(quad(lambda s: abs(s**(1.5) - 1)**(1.5) * s**((5.*d+mu)/2 - 11./4), 1,j)[0])
        except AttributeError:
            t =abs(quad(lambda s: abs(s**(1.5) - 1)**(1.5) * s**((5.*d+mu)/2 - 11./4), 1,r)[0])
        return  (3./(2**(5./4) * h**(5./2))) * t

    def spiral_loc(self,r,d=.5,h=.05):
        if d == .5:
            return self.mod_pi( (np.sign(r-1)/h)*(2./3)*(1 - r**(1.5) + 1.5*np.log(r)) )
        elif d == -1:
            return self.mod_pi( (np.sign(r-1)/h)*(2./3)*(1 - r**(-1.5) - 1.5*np.log(r)) )
        else:
            return self.mod_pi( (np.sign(r-1)/h)*( r**(d-.5)/(d-.5) - r**(d+1)/(d+1) - 3./( (2*d-1)*(d+1))))

    def convert_coords(self,r,mu=.5,d=.5,h=.05,**kargs):
        from scipy.integrate import quad
        phi,sig,_ = self.get_ring(r,**kargs)

        chi = sig*( 2**(.25) * h**(.5) * r**((1. - mu -3*d)/2) / np.sqrt(abs(r**(1.5)-1)))

        if d == .5:
            phi0 = (np.sign(r-1)/h)*(2./3)*(1 - r**(1.5) + 1.5*np.log(r))
        elif d == -1:
            phi0 = (np.sign(r-1)/h)*(2./3)*(1 - r**(-1.5) - 1.5*np.log(r))
        else:
            phi0 = (np.sign(r-1)/h)*( r**(d-.5)/(d-.5) - r**(d+1)/(d+1) - 3./( (2*d-1)*(d+1)))

        #phi = np.array([self.mod_pi(p) for p in phi])

        eta = -1.5/h * (phi - self.mod_pi(phi0))
        inds = np.argsort(eta)


        t = (3./(2**(5./4) * h**(5./2))) * abs(quad(lambda s: abs(s**(1.5) - 1)**(1.5) * s**((5.*d+mu)/2 - 11./4), 1,r)[0])

        return t,eta[inds],chi[inds]
    def convert_coords_back(self,t,eta,chi,r=None,mu=.5,d=.5,h=.05,norm=1):
        from scipy.optimize import fsolve
        from scipy.integrate import quad
        if r is None:
            r = fsolve( lambda x: t - (3./(2**(5./4) * h**(5./2))) * abs(quad(lambda s: abs(s**(1.5) - 1)**(1.5) * s**((5.*d+mu)/2 - 11./4), 1,x)[0]),.8)[0]

        sig = chi/( 2**(.25) * h**(.5) * r**((1. - mu -3*d)/2) / np.sqrt(abs(r**(1.5)-1)))
        sig *= norm
        if d == .5:
            phi0 = (np.sign(r-1)/h)*(2./3)*(1 - r**(1.5) + 1.5*np.log(r))
        elif d == -1:
            phi0 = (np.sign(r-1)/h)*(2./3)*(1 - r**(-1.5) - 1.5*np.log(r))
        else:
            phi0 = (np.sign(r-1)/h)*( r**(d-.5)/(d-.5) - r**(d+1)/(d+1) - 3./( (2*d-1)*(d+1)))

        phi = eta/(-1.5/h) +  self.mod_pi(phi0)

        phi = np.array([self.mod_pi(p) for p in phi])

        inds = np.argsort(phi)
        #return phi,sig
        return phi[inds],sig[inds]
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




class Burgers():
    def __init__(self,x,u,cfl,out_times):
        self.x = x.copy()
        self.dx = np.diff(x)[0]
        self.t = out_times[0]
        self.output = np.zeros((len(u),len(out_times)))
        self.u = u.copy()
        self.output[:,0] = self.u
        self.cfl = cfl
        self.out_times = out_times
    def solve(self,tend):
        while self.t<tend:
            dt = self.update(self.u,self.cfl,self.dx,tend-self.t)
            self.t += dt
        return
    def solve_all(self):
        print("Solving from {:.2f} - {:.2f}".format(self.out_times[0],self.out_times[-1]))
        for i,tend in enumerate(self.out_times[1:]):
            while self.t<tend:
                dt = self.update(self.u,self.cfl,self.dx,tend-self.t)
                self.t += dt
            self.output[:,i+1] = self.u

        return self.output,self.out_times


    def riemann(self,u_left,u_right):
     # Exact Riemann solver for invscid Burgers equation
        vs = .5*(u_left + u_right)

        output = (u_left>u_right).astype(int)*((vs>=0).astype(int)*u_left \
                +(vs<0).astype(int)*u_right) \
                + (u_left<=u_right).astype(int)*((u_left>=0).astype(int)*u_left \
                + (u_right<=0).astype(int)*u_right)
        return output


    def get_timestep(self,u,cfl,dx,max_dt):

        smax = abs(u).max()
        if smax==0:
            smax = 1
        out_dt = cfl*dx/smax
        if out_dt > max_dt:
            out_dt = max_dt
        return out_dt

    def update(self,u,cfl,dx,max_dt):


        dt=self.get_timestep(u,cfl,dx,max_dt)
        dtdx = dt/dx
        flux = self.riemann(np.roll(u,1),u)
        #flux = self.riemann(u[:-1],u[1:])
        flux = .5*flux**2

        #u[1:-1] -= dtdx*(flux[1:] - flux[:-1])
        u -= dtdx*(np.roll(flux,-1) -flux)

        return dt



class BurgersL():
    def __init__(self,x,u,a,s,cfl,out_times,dt_min=1e-8):
        """
            u_t + ( u+ a(x,t) ) u_x + s(x,t) = 0
            and u is periodic.
        """

        self.x = x.copy()
        self.dx = np.diff(x)[0]
        self.t = out_times[0]
        self.output = np.zeros((len(u),len(out_times)))
        self.u = u.copy()
        self.output[:,0] = self.u
        self.cfl = cfl
        self.out_times = out_times
        self.afunc = a
        self.sfunc = s
        self.dt_min = dt_min

    def solve(self,tend):
        while self.t < tend:
            self.t = self.update(self.u,self.cfl,self.dx,self.t,tend-self.t)
        return
    def solve_all(self):
        for i,tend in enumerate(self.out_times[1:]):
            while self.t<tend:
                self.t = self.update(self.u,self.cfl,self.dx,self.t,tend-self.t)
            self.output[:,i+1] = self.u

        return self.output,self.out_times


    def riemann(self,u_left,u_right,a_left,a_right):

        ind =  ~np.isclose(u_left,u_right)

        if np.any( np.isclose(u_right-u_left,a_left-a_right)):
            print('Contact disconinuity!')

        output = u_left.copy()

        u_L= u_left[ind].copy()
        u_R = u_right[ind].copy()
        a_R = a_right[ind].copy()
        a_L = a_left[ind].copy()



        vs = .5*(u_L + u_R) +(a_R * u_R - a_L * u_L)/(u_R -u_L)

        output[ind] = (u_L>u_R).astype(int)*((vs>=0).astype(int)*u_L \
                +(vs<0).astype(int)*u_R) \
                + (u_L<u_R).astype(int)*(((u_L+a_L)>=0).astype(int)*u_L \
                + ((u_R+a_R)<=0).astype(int)*u_R \
                + ( ((u_R+a_R)>0)&((u_L+a_L)<0)).astype(int)*-.5*(a_L+a_R))
        return output
    def source_step(self,u,t,dt):
        du = np.array( [self.afunc(i,t+dt) - self.afunc(i,t) for i in self.x ])
        du += (dt/6.)*np.array([self.sfunc(i,t)+self.sfunc(i,t+.5*dt)+ self.sfunc(i,t+dt) for i in self.x])
        return -du

    def get_timestep(self,u,a,cfl,dx,max_dt):

        smax = abs(u+a).max()
        if smax==0:
            smax = 1
        return min(cfl*dx/smax,max_dt)

    def update(self,u,cfl,dx,t,max_dt):
        a = np.array([self.afunc(i,t) for i in self.x])
        dt=self.get_timestep(u,a,cfl,dx,max_dt)
        if dt < self.dt_min:
            print('Integrator not converging timestep of {:.2e} below {:.2e}'.format(dt, self.dt_min))
            return None
        u += self.source_step(u,t,.5*dt)
        a = np.array([self.afunc(i,t+.5*dt) for i in self.x])
        #a = np.zeros(u.shape)

        flux = self.riemann(np.roll(u,1),u,np.roll(a,1),a)
        flux = .5*flux**2 + .5*(a + np.roll(a,1))*flux

        u -= (dt/dx)*(np.roll(flux,-1)- flux)
        u += self.source_step(u,t+.5*dt,t + dt)
        t += dt
        return t



