import os
import analysis
import powerlaw
import numpy as np
import scipy as sp
from scipy import constants
from scipy.stats import maxwell
from astropy.table import Table
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib as mpl

# Create directories for the images if they don't already exist

if not os.path.exists('../Graphs'):
     os.mkdir('../Graphs')
if not os.path.exists('../Maps'):
     os.mkdir('../Maps')
if not os.path.exists('../Parameters'):
     os.mkdir('../Parameters')
if not os.path.exists('../Data'):
     os.mkdir('../Data')
if not os.path.exists('./Disk'):
     os.mkdir('./Disk')

########################################################################################################################################
#                                                           Galaxy Analysis                                                            #
########################################################################################################################################

# Load the data file and display it in a few graphs.

#NGC0628
analysis.graphs ('NGC0628','NGC0628_co21_cube_pbcor_props_clfind.fits')
analysis.mapgmc ('NGC0628','NGC0628_co21_cube_pbcor_props_clfind.fits','NGC0628_co21_cube_tpeak.fits')
NGC0628_mass_disk,NGC0628_rad_disk,NGC0628_alpha_disk,NGC0628_sigma_disk,NGC0628_lwo_disk,NGC0628_lw_disk = analysis.param ('NGC0628','NGC0628_co21_cube_pbcor_props_clfind.fits')
NGC0628_distance,NGC0628_inclination = analysis.data('NGC0628','NGC0628_co21_cube_pbcor_props_clfind.fits',n_bins=5)

#NGC1672
analysis.graphs ('NGC1672','NGC1672_co21_cube_pbcor_props_clfind.fits')
analysis.mapgmc ('NGC1672','NGC1672_co21_cube_pbcor_props_clfind.fits','NGC1672_co21_cube_tpeak.fits')
NGC1672_mass_disk,NGC1672_rad_disk,NGC1672_alpha_disk,NGC1672_sigma_disk,NGC1672_lwo_disk,NGC1672_lw_disk = analysis.param ('NGC1672','NGC1672_co21_cube_pbcor_props_clfind.fits',r_nuc=2)
NGC1672_distance,NGC1672_inclination = analysis.data('NGC1672','NGC1672_co21_cube_pbcor_props_clfind.fits',n_bins=3,r_nuc=2)

#NGC3351
analysis.graphs ('NGC3351','NGC3351_co21_cube_pbcor_props_clfind.fits')
analysis.mapgmc ('NGC3351','NGC3351_co21_cube_pbcor_props_clfind.fits','NGC3351_co21_cube_tpeak.fits')
NGC3351_mass_disk,NGC3351_rad_disk,NGC3351_alpha_disk,NGC3351_sigma_disk,NGC3351_lwo_disk,NGC3351_lw_disk = analysis.param ('NGC3351','NGC3351_co21_cube_pbcor_props_clfind.fits',r_nuc=1)
NGC3351_distance,NGC3351_inclination = analysis.data('NGC3351','NGC3351_co21_cube_pbcor_props_clfind.fits',n_bins=3,r_nuc=1)

#NGC3627
analysis.graphs ('NGC3627','NGC3627_co21_cube_pbcor_props_clfind.fits')
analysis.mapgmc ('NGC3627','NGC3627_co21_cube_pbcor_props_clfind.fits','NGC3627_co21_cube_tpeak.fits')
NGC3627_mass_disk,NGC3627_rad_disk,NGC3627_alpha_disk,NGC3627_sigma_disk,NGC3627_lwo_disk,NGC3627_lw_disk = analysis.param ('NGC3627','NGC3627_co21_cube_pbcor_props_clfind.fits',r_nuc=1)
NGC3627_distance,NGC3627_inclination = analysis.data('NGC3627','NGC3627_co21_cube_pbcor_props_clfind.fits',n_bins=5,r_nuc=1)

#NGC4254
analysis.graphs ('NGC4254','NGC4254_co21_cube_pbcor_props_clfind.fits')
analysis.mapgmc ('NGC4254','NGC4254_co21_cube_pbcor_props_clfind.fits','NGC4254_co21_cube_tpeak.fits')
NGC4254_mass_disk,NGC4254_rad_disk,NGC4254_alpha_disk,NGC4254_sigma_disk,NGC4254_lwo_disk,NGC4254_lw_disk = analysis.param ('NGC4254','NGC4254_co21_cube_pbcor_props_clfind.fits',r_nuc=2)
NGC4254_distance,NGC4254_inclination = analysis.data('NGC4254','NGC4254_co21_cube_pbcor_props_clfind.fits',n_bins=4,r_nuc=2)

#NGC4303
analysis.graphs ('NGC4303','NGC4303_co21_cube_pbcor_props_clfind.fits')
analysis.mapgmc ('NGC4303','NGC4303_co21_cube_pbcor_props_clfind.fits','NGC4303_co21_cube_tpeak.fits')
NGC4303_mass_disk,NGC4303_rad_disk,NGC4303_alpha_disk,NGC4303_sigma_disk,NGC4303_lwo_disk,NGC4303_lw_disk = analysis.param ('NGC4303','NGC4303_co21_cube_pbcor_props_clfind.fits',r_nuc=1)
NGC4303_distance,NGC4303_inclination = analysis.data('NGC4303','NGC4303_co21_cube_pbcor_props_clfind.fits',n_bins=4,r_nuc=1)

#NGC4321
analysis.graphs ('NGC4321','NGC4321_co21_cube_pbcor_props_clfind.fits')
analysis.mapgmc ('NGC4321','NGC4321_co21_cube_pbcor_props_clfind.fits','NGC4321_co21_cube_tpeak.fits')
NGC4321_mass_disk,NGC4321_rad_disk,NGC4321_alpha_disk,NGC4321_sigma_disk,NGC4321_lwo_disk,NGC4321_lw_disk = analysis.param ('NGC4321','NGC4321_co21_cube_pbcor_props_clfind.fits',r_nuc=1)
NGC4321_distance,NGC4321_inclination = analysis.data('NGC4321','NGC4321_co21_cube_pbcor_props_clfind.fits',n_bins=4,r_nuc=1)

#NGC4535
analysis.graphs ('NGC4535','NGC4535_co21_cube_pbcor_props_clfind.fits')
analysis.mapgmc ('NGC4535','NGC4535_co21_cube_pbcor_props_clfind.fits','NGC4535_co21_cube_tpeak.fits')
NGC4535_mass_disk,NGC4535_rad_disk,NGC4535_alpha_disk,NGC4535_sigma_disk,NGC4535_lwo_disk,NGC4535_lw_disk = analysis.param ('NGC4535','NGC4535_co21_cube_pbcor_props_clfind.fits',r_nuc=1)
NGC4535_distance,NGC4535_inclination = analysis.data('NGC4535','NGC4535_co21_cube_pbcor_props_clfind.fits',n_bins=3,r_nuc=1)

#NGC5068
analysis.graphs ('NGC5068','NGC5068_co21_cube_pbcor_props_clfind.fits')
analysis.mapgmc ('NGC5068','NGC5068_co21_cube_pbcor_props_clfind.fits','NGC5068_co21_cube_tpeak.fits')
NGC5068_mass_disk,NGC5068_rad_disk,NGC5068_alpha_disk,NGC5068_sigma_disk,NGC5068_lwo_disk,NGC5068_lw_disk = analysis.param ('NGC5068','NGC5068_co21_cube_pbcor_props_clfind.fits')
NGC5068_distance,NGC5068_inclination = analysis.data('NGC5068','NGC5068_co21_cube_pbcor_props_clfind.fits',n_bins=4)

#NGC6744
analysis.graphs ('NGC6744','NGC6744_co21_cube_pbcor_props_clfind.fits')
analysis.mapgmc ('NGC6744','NGC6744_co21_cube_pbcor_props_clfind.fits','NGC6744_co21_cube_tpeak.fits')
NGC6744_mass_disk,NGC6744_rad_disk,NGC6744_alpha_disk,NGC6744_sigma_disk,NGC6744_lwo_disk,NGC6744_lw_disk = analysis.param ('NGC6744','NGC6744_co21_cube_pbcor_props_clfind.fits')
NGC6744_distance,NGC6744_inclination = analysis.data('NGC6744','NGC6744_co21_cube_pbcor_props_clfind.fits',n_bins=3)

########################################################################################################################################
#                                                      Anderson-Darling Analysis                                                       #
########################################################################################################################################

# Create tuples of the galactic disk GMC and galaxy properties.

distance = (NGC0628_distance,NGC1672_distance,NGC3351_distance,NGC3627_distance,NGC4254_distance,
            NGC4303_distance,NGC4321_distance,NGC4535_distance,NGC5068_distance,NGC6744_distance)
inclination = (NGC0628_inclination,NGC1672_inclination,NGC3351_inclination,NGC3627_inclination,NGC4254_inclination,
               NGC4303_inclination,NGC4321_inclination,NGC4535_inclination,NGC5068_inclination,NGC6744_inclination)
mass_dict = (NGC0628_mass_disk,NGC1672_mass_disk,NGC3351_mass_disk,NGC3627_mass_disk,NGC4254_mass_disk,
             NGC4303_mass_disk,NGC4321_mass_disk,NGC4535_mass_disk,NGC5068_mass_disk,NGC6744_mass_disk)
alpha_dict = (NGC0628_alpha_disk,NGC1672_alpha_disk,NGC3351_alpha_disk,NGC3627_alpha_disk,NGC4254_alpha_disk,
              NGC4303_alpha_disk,NGC4321_alpha_disk,NGC4535_alpha_disk,NGC5068_alpha_disk,NGC6744_alpha_disk)
sigma_dict = (NGC0628_sigma_disk,NGC1672_sigma_disk,NGC3351_sigma_disk,NGC3627_sigma_disk,NGC4254_sigma_disk,
              NGC4303_sigma_disk,NGC4321_sigma_disk,NGC4535_sigma_disk,NGC5068_sigma_disk,NGC6744_sigma_disk)
lwo_dict = (NGC0628_lwo_disk,NGC1672_lwo_disk,NGC3351_lwo_disk,NGC3627_lwo_disk,NGC4254_lwo_disk,
            NGC4303_lwo_disk,NGC4321_lwo_disk,NGC4535_lwo_disk,NGC5068_lwo_disk,NGC6744_lwo_disk)

# Create a template for the comparison tables.

column_names = ['Galaxies','NGC0628','NGC1672','NGC3351','NGC3627','NGC4254','NGC4303','NGC4321','NGC4535','NGC5068','NGC6744']
galaxy_names = ['NGC0628','NGC1672','NGC3351','NGC3627','NGC4254','NGC4303','NGC4321','NGC4535','NGC5068','NGC6744']
column_types = ['S7','f4','f4','f4','f4','f4','f4','f4','f4','f4','f4']
alpha_table = Table(names=column_names,dtype=column_types)
sigma_table = Table(names=column_names,dtype=column_types)
lwo_table = Table(names=column_names,dtype=column_types)

# Fill the tables with the results of the Anderson-Darling k-sample tests.

for i in range(len(galaxy_names)):
    alpha_table.add_row()
    sigma_table.add_row()
    lwo_table.add_row()
    alpha_table[-1]['Galaxies'] = galaxy_names[i]
    sigma_table[-1]['Galaxies'] = galaxy_names[i]
    lwo_table[-1]['Galaxies'] = galaxy_names[i]
    a = 0         #dummy variables for the unnecessary outputs from the tests
    b = [0]
    ca = 0
    cs = 0
    cl = 0
    for j in range(len(galaxy_names)):
         ai = alpha_dict[i]
         aj = alpha_dict[j]
         si = sigma_dict[i]
         sj = sigma_dict[j]
         li = lwo_dict[i]
         lj = lwo_dict[j]

         try:
              stat_a,b,c = sp.stats.anderson_ksamp([ai,aj])
         except OverflowError:
              stat_a,b,c = 50,None,None
         try:
              stat_s,b,c = sp.stats.anderson_ksamp([si,sj])
         except OverflowError:
              stat_s,b,c = 50,None,None
         try:
              stat_l,b,c = sp.stats.anderson_ksamp([li,lj])
         except OverflowError:
              stat_l,b,c = 50,None,None
         alpha_table[-1][galaxy_names[j]] = stat_a
         sigma_table[-1][galaxy_names[j]] = stat_s
         lwo_table[-1][galaxy_names[j]] = stat_l

alpha_table.write('./Disk/alpha comparison.FITS',overwrite=True)
sigma_table.write('./Disk/sigma comparison.FITS',overwrite=True)
lwo_table.write('./Disk/lwo comparison.FITS',overwrite=True)

########################################################################################################################################
#                                                      Galactic Disk GMC Analysis                                                      #
########################################################################################################################################

# Concatenate all the galactic disk data.

mass_disk = np.concatenate((NGC0628_mass_disk,NGC1672_mass_disk,NGC3351_mass_disk,NGC3627_mass_disk,NGC4254_mass_disk,
                            NGC4303_mass_disk,NGC4321_mass_disk,NGC4535_mass_disk,NGC5068_mass_disk,NGC6744_mass_disk))
rad_disk = np.concatenate((NGC0628_rad_disk,NGC1672_rad_disk,NGC3351_rad_disk,NGC3627_rad_disk,NGC4254_rad_disk,
                           NGC4303_rad_disk,NGC4321_rad_disk,NGC4535_rad_disk,NGC5068_rad_disk,NGC6744_rad_disk))
alpha_disk = np.concatenate((NGC0628_alpha_disk,NGC1672_alpha_disk,NGC3351_alpha_disk,NGC3627_alpha_disk,NGC4254_alpha_disk,
                             NGC4303_alpha_disk,NGC4321_alpha_disk,NGC4535_alpha_disk,NGC5068_alpha_disk,NGC6744_alpha_disk))
sigma_disk = np.concatenate((NGC0628_sigma_disk,NGC1672_sigma_disk,NGC3351_sigma_disk,NGC3627_sigma_disk,NGC4254_sigma_disk,
                             NGC4303_sigma_disk,NGC4321_sigma_disk,NGC4535_sigma_disk,NGC5068_sigma_disk,NGC6744_sigma_disk))
lwo_disk = np.concatenate((NGC0628_lwo_disk,NGC1672_lwo_disk,NGC3351_lwo_disk,NGC3627_lwo_disk,NGC4254_lwo_disk,
                           NGC4303_lwo_disk,NGC4321_lwo_disk,NGC4535_lwo_disk,NGC5068_lwo_disk,NGC6744_lwo_disk))
lw_disk = np.concatenate((NGC0628_lw_disk,NGC1672_lw_disk,NGC3351_lw_disk,NGC3627_lw_disk,NGC4254_lw_disk,
                           NGC4303_lw_disk,NGC4321_lw_disk,NGC4535_lw_disk,NGC5068_lw_disk,NGC6744_lw_disk))
alpha_disk_alike = np.concatenate((NGC4303_alpha_disk,NGC4321_alpha_disk,NGC4535_alpha_disk,NGC5068_alpha_disk,NGC6744_alpha_disk))
sigma_disk_alike = np.concatenate((NGC1672_sigma_disk,NGC4303_sigma_disk,NGC4321_sigma_disk,NGC5068_sigma_disk,NGC6744_sigma_disk))
lwo_disk_alike = np.concatenate((NGC3351_lwo_disk,NGC4254_lwo_disk,NGC4321_lwo_disk,NGC6744_lwo_disk))

# Histogram of the parameters.

parameter_hist = plt.figure(figsize=(18,6))

#alpha
ax = plt.subplot(1,3,1)
n,bins,patches = plt.hist(alpha_disk,bins=50,normed=True,alpha=0.6,facecolor='green')
alpha_mu = alpha_disk.mean()
alpha_std = alpha_disk.std()
y = mlab.normpdf(bins,alpha_mu,alpha_std)
plt.plot(bins,y,'k--')
mpl.rc('xtick',labelsize=16)
mpl.rc('ytick',labelsize=16)
plt.text(0.01,0.99,'(a)',ha='left',va='top',transform=ax.transAxes,fontsize=16)
plt.text(0.99,0.99,'$\mu$ = %5.4f\n$\sigma$ = %5.4f'%(alpha_mu,alpha_std),ha='right',va='top',transform=ax.transAxes,fontsize=16)
#plt.title('Galactic Disk Virial Parameter',fontsize=20)
plt.xlabel(r'$log_\mathrm{10}(\alpha)$',fontsize=20)

#sigma
ax = plt.subplot(1,3,2)
n,bins,patches = plt.hist(sigma_disk,bins=50,normed=True,alpha=0.6,facecolor='green')
sigma_mu = sigma_disk.mean()
sigma_std = sigma_disk.std()
y = mlab.normpdf(bins,sigma_mu,sigma_std)
plt.plot(bins,y,'k--')
mpl.rc('xtick',labelsize=16)
mpl.rc('ytick',labelsize=16)
plt.text(0.01,0.99,'(b)',ha='left',va='top',transform=ax.transAxes,fontsize=16)
plt.text(0.99,0.99,'$\mu$ = %5.4f\n$\sigma$ = %5.4f'%(sigma_mu,sigma_std),ha='right',va='top',transform=ax.transAxes,fontsize=16)
#plt.title('Galactic Disk Surface density',fontsize=20)
plt.xlabel(r'$log_\mathrm{10}(\Sigma)$',fontsize=20)

#lwo
ax = plt.subplot(1,3,3)
n,bins,patches = plt.hist(lwo_disk,bins=50,normed=True,alpha=0.6,facecolor='green')
mw_param = maxwell.fit(lwo_disk)
lwo_mu = maxwell.mean(*mw_param)
lwo_std = maxwell.std(*mw_param)
y = np.linspace(0,np.max(lwo_disk),1000)
z = maxwell.pdf(y,*mw_param)
plt.plot(y,z,'k--')
mpl.rc('xtick',labelsize=16)
mpl.rc('ytick',labelsize=16)
plt.text(0.01,0.99,'(c)',ha='left',va='top',transform=ax.transAxes,fontsize=16)
plt.text(0.99,0.99,'$\mu$ = %5.4f\n$\sigma$ = %5.4f'%(lwo_mu,lwo_std),ha='right',va='top',transform=ax.transAxes,fontsize=16)
#plt.title('Galactic Disk Normalized Linewidth',fontsize=20)
plt.xlabel(r'$\sigma_{o}$',fontsize=20)

plt.tight_layout()
plt.savefig('./Disk/parameter_hist.png')
plt.close()

# Star formation rate and depletion time histograms.

SFR = []
depletion = []
T = 20
mu = 0.61
alpha_disk = 10**alpha_disk

mach = 2.0727*(lw_disk)*((mu/0.61*20/T)**(0.5))
t_ff = (16.572*10**6)*(rad_disk**(1.5))*((mass_disk)**(-0.5))
SFR = 0.014*mass_disk/t_ff*((alpha_disk/1.3)**(-0.68))*((mach/100)**(-0.32))
depletion = mass_disk/SFR
SFR = np.delete(SFR,np.isnan(SFR))
SFR = np.log10(SFR)
depletion = np.delete(depletion,np.isnan(depletion))
depletion = np.log10(depletion)

formation_hist = plt.figure(figsize=(12,6))

ax = plt.subplot(1,2,1)
n,bins,patches = plt.hist(SFR,bins=50,normed=True,log=False,alpha=0.6,facecolor='blue')
SFR_mu = SFR.mean()
SFR_std = SFR.std()
y = mlab.normpdf(bins,SFR_mu,SFR_std)
plt.plot(bins,y,'r--')
mpl.rc('xtick',labelsize=16)
mpl.rc('ytick',labelsize=16)
plt.text(0.01,0.99,'(a)',ha='left',va='top',transform=ax.transAxes,fontsize=16)
plt.text(0.99,0.99,'$\mu$ = %5.4f\n$\sigma$ = %5.4f'%(SFR_mu,SFR_std),ha='right',va='top',transform=ax.transAxes,fontsize=16)
#plt.title('Galactic Disk Star Formation Rate',fontsize=20)
plt.xlabel(r'$log_\mathrm{10}(SFR) (\frac{M_\mathrm{\odot}}{yr})$',fontsize=20)

ax = plt.subplot(1,2,2)
n,bins,patches = plt.hist(depletion,bins=50,normed=True,log=False,alpha=0.6,facecolor='blue')
depletion_mu = depletion.mean()
depletion_std = depletion.std()
y = mlab.normpdf(bins,depletion_mu,depletion_std)
plt.plot(bins,y,'r--')
mpl.rc('xtick',labelsize=16)
mpl.rc('ytick',labelsize=16)
plt.text(0.01,0.99,'(b)',ha='left',va='top',transform=ax.transAxes,fontsize=16)
plt.text(0.99,0.99,'$\mu$ = %5.4f\n$\sigma$ = %5.4f'%(depletion_mu,depletion_std),ha='right',va='top',transform=ax.transAxes,fontsize=16)
#plt.title('Galactic Disk Depletion Rate',fontsize=20)
plt.xlabel(r'$log_\mathrm{10}(t_\mathrm{d}) (yr)$',fontsize=20)

plt.tight_layout()
plt.savefig('./Disk/formation_hist.png')
plt.close()

# Create a table of some galaxy properties.

r_nuc = [0,1,1,1,2,1,1,1,0,0]
properties = ['Galaxy','Distance (Mpc)','Inclination (deg)','Number of GMCs','Number of disk GMCs','R_nuc (kpc)']
column_types = ['S7','f4',int,int,int,int]
galaxy_table = Table(names=properties,dtype=column_types)

for i in range(len(galaxy_names)):
     datatable = Table.read(galaxy_names[i]+'_co21_cube_pbcor_props_clfind.fits')
     mass = datatable['MASS_EXTRAP']
     galaxy_table.add_row()
     galaxy_table[-1]['Galaxy'] = galaxy_names[i]
     galaxy_table[-1]['Distance (Mpc)'] = distance[i]
     galaxy_table[-1]['Inclination (deg)'] = inclination[i]
     galaxy_table[-1]['Number of GMCs'] = len(mass)
     galaxy_table[-1]['Number of disk GMCs'] = len(mass_dict[i])
     galaxy_table[-1]['R_nuc (kpc)'] = r_nuc[i]

galaxy_table.write('./Disk/galaxy_table.FITS',overwrite=True)

# Fit all of the galactic disk GMC data in the power law and truncated-power-law distributions.

myfit = powerlaw.Fit(mass_disk)
R, p = myfit.distribution_compare('power_law','truncated_power_law')
fig = myfit.truncated_power_law.plot_ccdf(label='Truncated\nPower Law')
myfit.power_law.plot_ccdf(label='Power Law',ax=fig)
myfit.plot_ccdf(drawstyle='steps',label='Data',ax=fig)

plt.legend(loc=0)
#plt.title('Galactic Disk GMC Mass Distribution')
mpl.rc('xtick',labelsize=16)
mpl.rc('ytick',labelsize=16)
plt.xlabel(r'$M_\mathrm{\odot}$',fontsize=20)
plt.ylabel('CCDF',fontsize=20)
plt.text(0.01,0.01,'$R$ = %5.4f\n$p$ = %5.4f'%(R,p),ha='left',va='bottom',transform=fig.transAxes,fontsize=16)
plt.savefig('./Disk/disk_power_law_mass.png')
plt.close()
