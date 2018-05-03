import numpy as np

class loc_emiss:
    def __init__(self):

        self.energies=np.array([ 56.56, 72.39, 92.65, 118.57, 151.74, 194.2 , 248.54, 318.08, 407.08, 520.99, 666.76, 853.32, 1092.08, 1397.65, 1788.72, 2289.2 , 2929.73, 3749.47, 4798.58, 6141.23, 7859.57, 10058.69, 12873.13, 16475.06, 23852.91, 39068.54])
        self.ebounds=np.logspace(np.log10(50),np.log10(50000.),29)
        self.ebounds=np.delete(self.ebounds,-4)
        self.ebounds=np.delete(self.ebounds,-2)
        self.energies=np.append(self.energies,120000.)
        self.ebounds=np.append(self.ebounds,np.power(10,2*np.log10(120000)-np.log10(self.ebounds[-1])))
        
        self.values=np.array([ 4.07e-25, 6.28e-25, 8.19e-25, 1.09e-24, 1.38e-24, 1.73e-24, 2.02e-24, 2.25e-24, 2.40e-24, 2.49e-24, 2.49e-24, 2.43e-24, 2.32e-24, 2.15e-24, 1.99e-24, 1.72e-24, 1.44e-24, 1.21e-24, 1.04e-24, 8.39e-25, 6.87e-25, 5.57e-25, 4.60e-25, 3.65e-25, 3.16e-25, 2.11e-25,8.2418e-26])

        stat_err=np.array([ 1.73e-26, 1.73e-26, 1.55e-26, 1.48e-26, 1.32e-26, 1.29e-26, 1.21e-26, 1.30e-26, 1.29e-26, 1.31e-26, 1.40e-26, 1.44e-26, 1.50e-26, 1.57e-26, 1.60e-26, 1.69e-26, 1.58e-26, 1.78e-26, 1.86e-26, 1.87e-26, 2.01e-26, 2.13e-26, 2.32e-26, 2.28e-26, 2.03e-26, 1.83e-26,0])
        sys=np.array([ 7.82e-26, 1.02e-25, 1.09e-25, 1.26e-25, 1.41e-25, 1.55e-25, 1.60e-25, 1.58e-25, 1.46e-25, 1.33e-25, 1.34e-25, 1.40e-25, 1.43e-25, 1.42e-25, 1.42e-25, 1.35e-25, 1.28e-25, 1.23e-25, 1.16e-25, 1.03e-25, 8.91e-26, 7.45e-26, 5.84e-26, 4.52e-26, 3.26e-26, 2.25e-26,0])

        self.values=self.values/(self.energies*self.energies)
        self.errs=(np.sqrt(np.power(stat_err,2)+np.power(sys,2)))/(self.energies*self.energies)


    def PLInt(self,prefactor,s_prefactor,gamma,s_gamma,emin,emax):

        integral=prefactor/(-gamma+1)*(np.power(emax,-gamma+1)-np.power(emin,-gamma+1))

        dIdgamma=-integral*(-gamma+1)/np.power(-gamma+1,2)-prefactor/(-gamma+1)*(np.power(emax,-gamma+1)*np.log(emax)-np.power(emin,-gamma+1)*np.log(emin))
        s_I=np.sqrt(np.power(integral/prefactor*s_prefactor,2)+np.power(dIdgamma*s_gamma,2))

        return integral, s_I

    def PLCoeff(self,e1,f1,s_f1,e2,f2,s_f2):

        x1=np.log(e1)
        y1=np.log(f1)
        x2=np.log(e2)
        y2=np.log(f2)

        gamma=-(y2-y1)/(x2-x1)
        prefactor=y1+gamma*x1
        prefactor=np.exp(prefactor)

        s_gamma=1/np.log(e2/e1)*np.sqrt(np.power(s_f1/f1,2)+np.power(s_f2/f2,2))
        s_prefactor=prefactor*np.sqrt(np.power(s_f1/f1,2)+np.power(np.log(e1)*s_gamma/gamma,2))

        return prefactor, s_prefactor, gamma, s_gamma

    def int_emiss(self,Emin,Emax):
        tot=0.
        en_emiss=0
        errs=np.array([])
        bins=0

        for s, value in enumerate(self.values):
            if self.ebounds[s]>Emin and self.ebounds[s+1]<Emax:
                tot+=value*(self.ebounds[s+1]-self.ebounds[s])
                en_emiss+=self.energies[s]*value*(self.ebounds[s+1]-self.ebounds[s])
                errs=np.append(errs,self.errs[s]*(self.ebounds[s+1]-self.ebounds[s]))
                bins+=1
            else:
                if self.ebounds[s+1]>Emin and self.ebounds[s]<Emin:
                    N0, sN0, g, sg = self.PLCoeff(self.energies[s],self.values[s],self.errs[s],self.energies[s+1],self.values[s+1],self.errs[s+1])
                    a=value*(self.ebounds[s+1]-self.ebounds[s])#integral over the whole energy bin
                    sa=self.errs[s]*(self.ebounds[s+1]-self.ebounds[s])
                    b, sb=self.PLInt(N0,sN0,g,sg,self.ebounds[s],Emin)
                    tot+=a-b
                    en_emiss+=(a-b)*self.energies[s]
                    errs=np.append(errs,np.sqrt(sa**2+sb**2))
                    bins+=1
                if self.ebounds[s]<Emax and self.ebounds[s+1]>Emax:
                    N0, sN0, g, sg = self.PLCoeff(self.energies[s-1],self.values[s-1],self.errs[s-1],self.energies[s],self.values[s],self.errs[s])
                    a=value*(self.ebounds[s+1]-self.ebounds[s])#integral over the whole energy bin
                    sa=self.errs[s]*(self.ebounds[s+1]-self.ebounds[s])
                    b, sb=self.PLInt(N0,sN0,g,sg,Emax,self.ebounds[s+1])
                    tot+=a-b
                    en_emiss+=(a-b)*self.energies[s]
                    errs=np.append(errs,np.sqrt(sa**2+sb**2))
                    bins+=1

        errtot=np.sqrt(np.sum(np.power(errs,2)))
        return tot, errtot








