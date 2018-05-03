import matplotlib.pyplot as plt
import os
import numpy as np
import P7REP_emiss as loc_emiss
import pdb
import scipy.stats as sp

datadir="/users-data/tsalgues/Fermi/Gasmaps/10GeV/"
nreg=4

regions=['Outer region','Loc-Inter','Sgr-Perseus','Scutum Tan']
cols=['b','g','r','m']

#retrieve emissivity normalization coefficients from output files
nebins=len([dirname for dirname, dirnames, filenames in os.walk(datadir) if 'bin_' in dirname])
emiss=np.zeros([4,nreg,nebins])#four types of emissivities are HI, CO, DNM, HII
emiss_err=np.zeros([4,nreg,nebins])
emin=np.zeros(nebins)
emax=np.zeros(nebins)

s=0
for dirname, dirnames, filenames in os.walk(datadir):
    if "bin_" in dirname:
        try:
            results=open(dirname+"/output.dat").readlines()
#            print s, "reading file", dirname+"/output.dat"
            emin[s]=float(results[0].split(" ")[0])
            emax[s]=float(results[0].split(" ")[1])
#            print emin[s], emax[s]
            for entry in results:
                for reg in range(nreg):
                    if "A-"+str(reg) in entry:
                        emiss[0,reg,s]=float(entry.split(" ")[1])
                        emiss_err[0,reg,s]=float(entry.split(" ")[2])
#                        print emiss[0,reg,s]
                    if "B-"+str(reg) in entry:
                        emiss[1,reg,s]=float(entry.split(" ")[1])
                        emiss_err[1,reg,s]=float(entry.split(" ")[2])
                    if "C-"+str(reg) in entry:
                        emiss[2,reg,s]=float(entry.split(" ")[1])
                        emiss_err[2,reg,s]=float(entry.split(" ")[2])
##                  if "D-"+str(reg) in entry:
##                      emiss[3,reg,s]=float(entry.split(" ")[1])
##                      emiss_err[3,reg,s]=float(entry.split(" ")[2])
            s+=1
        except:
            pass

emean=np.sqrt(emin*emax)
emin=emin[emean>0]  
emax=emax[emean>0]    
emiss=emiss[:,:,emean>0]
emiss_err=emiss_err[:,:,emean>0]
emean=emean[emean>0]

emin=emin[emean<10000]  
emax=emax[emean<10000]    
emiss=emiss[:,:,emean<10000]
emiss_err=emiss_err[:,:,emean<10000]
emean=emean[emean<10000]   

#convert normalization into absolute emissivities
locem=loc_emiss.loc_emiss()
locem_val=np.zeros(len(emean))
for s in range(len(emean)):
    locem_val[s] = locem.int_emiss(emin[s],emax[s])[0]

emiss*=locem_val[np.newaxis,np.newaxis,:]
emiss_err*=locem_val[np.newaxis,np.newaxis,:]

emiss[1:3,:,:]*=1.e20
emiss_err[1:3,:,:]*=1.e20

#print 'emiss[0,:,:]\n'
#print emiss[0,1:,:]
#print '\n\n'
#print 'emiss[1,:,:]\n'
#print emiss[1,1:,:]
#print '\n\n'
#print 'emiss[2,:,:]\n'
#print emiss[2,1:,:]
#print '\n\n'


#print 'Weighted Average Emissivity HI'
#for reg in range(nreg):
#    print ': ' + str(reg)
#    print np.sqrt((np.sum(emiss_err[0,reg,:]**(-2)))**(-1))
#    print np.average(emiss[0,reg,:], weights=(emiss_err[0,reg,:])**(-2))
#
#print '\n \n'
#
#print 'Weighted Average Emissivity CO'
#for reg in range(1,nreg):
#    print 'region: ' + str(reg)
#    print np.sqrt((np.sum(emiss_err[1,reg,:]**(-2)))**(-1))
#    print np.average(emiss[1,reg,:], weights=(emiss_err[1,reg,:])**(-2))
#
#print '\n \n'
#
#print 'Weighted Average Emissivity Dust'
#for reg in range(1,nreg):
#    print 'region: ' + str(reg)
#    print np.sqrt((np.sum(emiss_err[2,reg,:]**(-2)))**(-1))
#    print np.average(emiss[2,reg,:], weights=(emiss_err[2,reg,:])**(-2))
#    
#print '\n \n'

xco=emiss[1,1:,:]/(2*emiss[0,1:,:]) # Plot X_co
xco*=1.e-20
xco_err=xco*np.sqrt(np.power(emiss_err[0,1:,:]/emiss[0,1:,:],2)+np.power(emiss_err[1,1:,:]/emiss[1,1:,:],2)) # erreur sans covariance
#print '\n\n Xco:'
#print xco
#print '\n'
#print xco_err
#print '\n'

xdust=emiss[0,1:,:]/emiss[2,1:,:] # Plot Opacite specifique
xdust*=1.e20
xdust_err=xdust*np.sqrt(np.power(emiss_err[0,1:,:]/emiss[0,1:,:],2)+np.power(emiss_err[2,1:,:]/emiss[2,1:,:],2)) # erreur sans covariance
#print xdust
#print '\n'
#print xdust_err
#print '\n'

energy_bins=[str(emin[0]) + ' to ' + str(emax[1]),str(emin[1]) + ' to ' + str(emax[2]),str(emin[2]) + ' to ' + str(emax[3]),str(emin[3]) + ' to ' + str(emax[0]),str(emin[4]) + ' to ' + str(emax[4])]
#print energy_bins
#print '\n'

#print xco[0][1]
#print '\n'
#print xco[0,1]
#print '\n'
#for i in range(3):
#    print 'ligne no.', i+1
#    for j in range(5):
#        print xco[i][j]
#    print '\n'


#plot Xco
#print 'Conversion Factor'
ax0=plt.figure('Xco')
for reg in range(nreg-1):
    xco_mean=[np.average(xco[reg,:],weights=(xco_err[reg,:])**(-2))]*len(emean)
    print 'Region: ' + str(reg+1)
    print xco_mean
    print np.sqrt((np.sum(xco_err[reg,:]**(-2)))**(-1))
    print '\n'
    
    ax0=plt.errorbar(emean,xco[reg,:],yerr=xco_err[reg,:],xerr=(emax[:]-emin[:])/2.,fmt='o',markersize=0,linewidth=2,capsize=0,color=cols[reg+1],label="{}".format(regions[reg+1]))
    ax0=plt.plot(emean,xco_mean, '--', color=cols[reg+1])

plt.xscale('log')
plt.xlim((0,20000))
plt.xlabel("$E$ (MeV)",fontsize=30)
plt.tick_params(axis='both', labelsize=30)
plt.ylabel("$X_{CO}$ ($10^{20} H_2.cm^{-2}.(K.km.s^{-1})^{-1}$)",fontsize=30)
plt.legend()

print '\n \n'

#plot Dust opacity
#print 'Dust Opacity'
ax0=plt.figure('Dust opacity')
for reg in range(nreg-1):
    xdust_mean=[np.average(xdust[reg,:],weights=(xdust_err[reg,:])**(-2))]*len(emean)
    print 'Region: ' + str(reg+1)
    print xdust_mean
    print np.sqrt((np.sum(xdust_err[reg,:]**(-2)))**(-1))
    print '\n'
    
    ax0=plt.errorbar(emean,xdust[reg,:],yerr=xdust_err[reg,:],xerr=(emax[:]-emin[:])/2.,fmt='o',markersize=0,linewidth=2,capsize=0,color=cols[reg+1],label="{}".format(regions[reg+1]))
    ax0=plt.plot(emean,xdust_mean, '--', color=cols[reg+1])

plt.xscale('log')
plt.xlim(xmin=-100,xmax=20000)
plt.xlabel("$E$ (MeV)",fontsize=30)
plt.tick_params(axis='both', labelsize=30)
plt.ylabel("$\sigma$ ($10^{-26} cm^2.H^{-1}$)",fontsize=30)
plt.legend()

#print 'Weighted Average Xco'
#for reg in range(nreg-1):
#    print 'region: ' + str(reg+1)
#    print np.average(xco[reg,:],weights=(xco_err[reg,:])**(-2))
#
#print '\n \n'
#
#print 'Weighted Average Xdust'
#for reg in range(nreg-1):
#    print 'region: ' + str(reg+1)
#    print np.average(xdust[reg,:],weights=(xdust_err[reg,:])**(-2))

		
#plot HI emissivities
fig0=plt.figure('HI emissivities')
fig0.subplots_adjust(left=0.12,bottom=0.12,right=0.96,top=0.96)
ax0=fig0.add_subplot(111)
ax0.errorbar(locem.energies,locem.energies*locem.energies*locem.values,yerr=locem.energies*locem.energies*locem.errs,xerr=[locem.energies-locem.ebounds[:-1],locem.ebounds[1:]-locem.energies],fmt='o',markersize=0,linewidth=3,capsize=0,color='lightgrey',label='loc emissivity',zorder=0)
for reg in range(nreg):
    ax0.errorbar(emean[(emiss_err[0,reg,:]/emiss[0,reg,:])<1.],(emean*emean*emiss[0,reg,:]/(emax-emin))[(emiss_err[0,reg,:]/emiss[0,reg,:])<1.],yerr=(emean*emean*emiss_err[0,reg,:]/(emax-emin))[(emiss_err[0,reg,:]/emiss[0,reg,:])<1.],xerr=[(emean-emin)[(emiss_err[0,reg,:]/emiss[0,reg,:])<1.],(emax-emean)[(emiss_err[0,reg,:]/emiss[0,reg,:])<1.]],fmt='o',markersize=0,linewidth=2,capsize=0,color=cols[reg],label="{}".format(regions[reg]),zorder=1)
ax0.set_xscale('log')
ax0.set_yscale('log')
ax0.set_xlabel("$E$ (MeV)",fontsize=40)
ax0.xaxis.set_tick_params(labelsize=40)
ax0.set_ylabel("$E^2$ * q$_{HI}$ (MeV s$^{-1}$ sr$^{-1}$ HI$^{-1}$)",fontsize=40)
ax0.yaxis.set_tick_params(labelsize=40)


leg0=ax0.legend(loc=1, prop={'size': 30})
plt.show()

#plot scatter of other emissivities (CO,DNM,HII) vs HI
#fig1, ax1 = plt.subplots(2, sharex=True,figsize=(6,11))
#fig1.subplots_adjust(left=0.18,bottom=0.07,right=0.95,top=0.97,hspace=0.1,wspace=0.2)
#for s in range(1,3):
#    for reg in range(nreg):
#        if np.any(emiss[s,reg,:]>0):
#            ax1[s-1].errorbar(emiss[0,reg,:][(emiss_err[s,reg,:]/emiss[s,reg,:])<1.],emiss[s,reg,:][(emiss_err[s,reg,:]/emiss[s,reg,:])<1.],yerr=emiss_err[s,reg,:][(emiss_err[s,reg,:]/emiss[s,reg,:])<1.],xerr=emiss_err[0,reg,:][(emiss_err[s,reg,:]/emiss[s,reg,:])<1.],fmt='o',markersize=0,linewidth=2,capsize=0,color=cols[reg],label="region {}".format(reg))
#    if s==1:
#        ax1[s-1].set_ylabel(r"per $W_\mathrm{CO}$ (cm$^{-2}$ sr$^{-1}$ K$^{-1}$ km$^{-1}$)",fontsize=14)
#    elif s==2:
#        ax1[s-1].set_ylabel(r"per DNM A$_\mathrm{Ks}$ (cm$^{-2}$ s$^{-1}$ sr$^{-1}$ mag$^{-1}$)",fontsize=14)
#    elif s==3:
#        ax1[s-1].set_ylabel(r"per HII (s$^{-1}$ sr$^{-1}$ HII$^{-1}$)",fontsize=14)
#        ax1[s-1].set_xlabel(r"emissivity per HI (s$^{-1}$ sr$^{-1}$ HI$^{-1}$)",fontsize=14)
#    ax1[s-1].yaxis.set_tick_params(labelsize=14)
#    ax1[s-1].xaxis.set_tick_params(labelsize=14)
#    ax1[s-1].set_xscale('log')
#    ax1[s-1].set_yscale('log')
#    leg=ax1[s-1].legend(loc=4)

#fig1=plt.figure('CO vs HI')
#fig1.subplots_adjust(left=0.12,bottom=0.12,right=0.96,top=0.96)
#ax1=fig1.add_subplot(111)
#ax1.set_ylabel(r"per $W_\mathrm{CO}$ (cm$^{-2}$ sr$^{-1}$ K$^{-1}$ km$^{-1}$)",fontsize=14)
#ax1.set_xlabel(r"emissivity per HI (s$^{-1}$ sr$^{-1}$ HI$^{-1}$)",fontsize=14)
#ax1.yaxis.set_tick_params(labelsize=14)
#ax1.xaxis.set_tick_params(labelsize=14)
#ax1.set_xscale('log')
#ax1.set_yscale('log')
#s=1
#for reg in range(nreg):
#    if np.any(emiss[s,reg,:]>0):
#        ax1.errorbar(emiss[0,reg,:][(emiss_err[s,reg,:]/emiss[s,reg,:])<1.],emiss[s,reg,:][(emiss_err[s,reg,:]/emiss[s,reg,:])<1.],yerr=emiss_err[s,reg,:][(emiss_err[s,reg,:]/emiss[s,reg,:])<1.],xerr=emiss_err[0,reg,:][(emiss_err[s,reg,:]/emiss[s,reg,:])<1.],fmt='o',markersize=0,linewidth=2,capsize=0,color=cols[reg],label="region {}".format(reg))
#    
#leg=ax1.legend(loc=4)

#plt.show()

