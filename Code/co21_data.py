import os
import mkgraph
import powerlaw
import numpy as np
import scipy as sp
from scipy import constants
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

r_nuc = [0,2,1,1,2,1,1,1,0,0]
n_bins = [5,3,3,5,4,4,4,3,4,3]

# Load the data file and display it in a few graphs.
#NGC0628
#mkgraph.graphs ('NGC0628','NGC0628_co21_cube_pbcor_props_clfind.fits')
#mkgraph.mapgmc ('NGC0628','NGC0628_co21_cube_pbcor_props_clfind.fits','NGC0628_co21_cube_tpeak.fits')
NGC0628_distance,NGC0628_mass_disk,NGC0628_rad_disk,NGC0628_alpha_disk,NGC0628_sigma_disk,NGC0628_lwo_disk,NGC0628_lw_disk = mkgraph.param ('NGC0628','NGC0628_co21_cube_pbcor_props_clfind.fits')
#mkgraph.data('NGC0628','NGC0628_co21_cube_pbcor_props_clfind.fits',n_bins=5)

#NGC1672
#mkgraph.graphs ('NGC1672','NGC1672_co21_cube_pbcor_props_clfind.fits')
#mkgraph.mapgmc ('NGC1672','NGC1672_co21_cube_pbcor_props_clfind.fits','NGC1672_co21_cube_tpeak.fits')
NGC1672_distance,NGC1672_mass_disk,NGC1672_rad_disk,NGC1672_alpha_disk,NGC1672_sigma_disk,NGC1672_lwo_disk,NGC1672_lw_disk = mkgraph.param ('NGC1672','NGC1672_co21_cube_pbcor_props_clfind.fits',r_nuc=2)
#mkgraph.data('NGC1672','NGC1672_co21_cube_pbcor_props_clfind.fits',n_bins=3,r_nuc=2)

#NGC3351
#mkgraph.graphs ('NGC3351','NGC3351_co21_cube_pbcor_props_clfind.fits')
#mkgraph.mapgmc ('NGC3351','NGC3351_co21_cube_pbcor_props_clfind.fits','NGC3351_co21_cube_tpeak.fits')
NGC3351_distance,NGC3351_mass_disk,NGC3351_rad_disk,NGC3351_alpha_disk,NGC3351_sigma_disk,NGC3351_lwo_disk,NGC3351_lw_disk = mkgraph.param ('NGC3351','NGC3351_co21_cube_pbcor_props_clfind.fits',r_nuc=1)
#mkgraph.data('NGC3351','NGC3351_co21_cube_pbcor_props_clfind.fits',n_bins=3,r_nuc=1)

#NGC3627
#mkgraph.graphs ('NGC3627','NGC3627_co21_cube_pbcor_props_clfind.fits')
#mkgraph.mapgmc ('NGC3627','NGC3627_co21_cube_pbcor_props_clfind.fits','NGC3627_co21_cube_tpeak.fits')
NGC3627_distance,NGC3627_mass_disk,NGC3627_rad_disk,NGC3627_alpha_disk,NGC3627_sigma_disk,NGC3627_lwo_disk,NGC3627_lw_disk = mkgraph.param ('NGC3627','NGC3627_co21_cube_pbcor_props_clfind.fits',r_nuc=1)
#mkgraph.data('NGC3627','NGC3627_co21_cube_pbcor_props_clfind.fits',n_bins=5,r_nuc=1)

#NGC4254
#mkgraph.graphs ('NGC4254','NGC4254_co21_cube_pbcor_props_clfind.fits')
#mkgraph.mapgmc ('NGC4254','NGC4254_co21_cube_pbcor_props_clfind.fits','NGC4254_co21_cube_tpeak.fits')
NGC4254_distance,NGC4254_mass_disk,NGC4254_rad_disk,NGC4254_alpha_disk,NGC4254_sigma_disk,NGC4254_lwo_disk,NGC4254_lw_disk = mkgraph.param ('NGC4254','NGC4254_co21_cube_pbcor_props_clfind.fits',r_nuc=2)
#mkgraph.data('NGC4254','NGC4254_co21_cube_pbcor_props_clfind.fits',n_bins=4,r_nuc=2)

#NGC4303
#mkgraph.graphs ('NGC4303','NGC4303_co21_cube_pbcor_props_clfind.fits')
#mkgraph.mapgmc ('NGC4303','NGC4303_co21_cube_pbcor_props_clfind.fits','NGC4303_co21_cube_tpeak.fits')
NGC4303_distance,NGC4303_mass_disk,NGC4303_rad_disk,NGC4303_alpha_disk,NGC4303_sigma_disk,NGC4303_lwo_disk,NGC4303_lw_disk = mkgraph.param ('NGC4303','NGC4303_co21_cube_pbcor_props_clfind.fits',r_nuc=1)
#mkgraph.data('NGC4303','NGC4303_co21_cube_pbcor_props_clfind.fits',n_bins=4,r_nuc=1)

#NGC4321
#mkgraph.graphs ('NGC4321','NGC4321_co21_cube_pbcor_props_clfind.fits')
#mkgraph.mapgmc ('NGC4321','NGC4321_co21_cube_pbcor_props_clfind.fits','NGC4321_co21_cube_tpeak.fits')
NGC4321_distance,NGC4321_mass_disk,NGC4321_rad_disk,NGC4321_alpha_disk,NGC4321_sigma_disk,NGC4321_lwo_disk,NGC4321_lw_disk = mkgraph.param ('NGC4321','NGC4321_co21_cube_pbcor_props_clfind.fits',r_nuc=1)
#mkgraph.data('NGC4321','NGC4321_co21_cube_pbcor_props_clfind.fits',n_bins=4,r_nuc=1)

#NGC4535
#mkgraph.graphs ('NGC4535','NGC4535_co21_cube_pbcor_props_clfind.fits')
#mkgraph.mapgmc ('NGC4535','NGC4535_co21_cube_pbcor_props_clfind.fits','NGC4535_co21_cube_tpeak.fits')
NGC4535_distance,NGC4535_mass_disk,NGC4535_rad_disk,NGC4535_alpha_disk,NGC4535_sigma_disk,NGC4535_lwo_disk,NGC4535_lw_disk = mkgraph.param ('NGC4535','NGC4535_co21_cube_pbcor_props_clfind.fits',r_nuc=1)
#mkgraph.data('NGC4535','NGC4535_co21_cube_pbcor_props_clfind.fits',n_bins=3,r_nuc=1)

#NGC5068
#mkgraph.graphs ('NGC5068','NGC5068_co21_cube_pbcor_props_clfind.fits')
#mkgraph.mapgmc ('NGC5068','NGC5068_co21_cube_pbcor_props_clfind.fits','NGC5068_co21_cube_tpeak.fits')
NGC5068_distance,NGC5068_mass_disk,NGC5068_rad_disk,NGC5068_alpha_disk,NGC5068_sigma_disk,NGC5068_lwo_disk,NGC5068_lw_disk = mkgraph.param ('NGC5068','NGC5068_co21_cube_pbcor_props_clfind.fits')
#mkgraph.data('NGC5068','NGC5068_co21_cube_pbcor_props_clfind.fits',n_bins=4)

#NGC6744
#mkgraph.graphs ('NGC6744','NGC6744_co21_cube_pbcor_props_clfind.fits')
#mkgraph.mapgmc ('NGC6744','NGC6744_co21_cube_pbcor_props_clfind.fits','NGC6744_co21_cube_tpeak.fits')
NGC6744_distance,NGC6744_mass_disk,NGC6744_rad_disk,NGC6744_alpha_disk,NGC6744_sigma_disk,NGC6744_lwo_disk,NGC6744_lw_disk = mkgraph.param ('NGC6744','NGC6744_co21_cube_pbcor_props_clfind.fits')
#mkgraph.data('NGC6744','NGC6744_co21_cube_pbcor_props_clfind.fits',n_bins=3)

########################################################################################################################################
#                                                      Anderson-Darling Analysis                                                       #
#                                                           <<NOT WORKING>>                                                            #
########################################################################################################################################

# Create tuples of the galactic disk GMC and galaxy properties.
distance = (NGC0628_distance,NGC1672_distance,NGC3351_distance,NGC3627_distance,NGC4254_distance,
            NGC4303_distance,NGC4321_distance,NGC4535_distance,NGC5068_distance,NGC6744_distance)
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
         #print i,j
         #a,b,ca = sp.stats.anderson_ksamp([ai,aj])
         a,b,cs = sp.stats.anderson_ksamp([si,sj])
         #a,b,cl = sp.stats.anderson_ksamp([li,lj])
         alpha_table[-1][galaxy_names[j]] = ca
         sigma_table[-1][galaxy_names[j]] = cs
         lwo_table[-1][galaxy_names[j]] = cl

alpha_table.write('./Disk/alpha comparison.FITS',overwrite=True)
sigma_table.write('./Disk/sigma comparison.FITS',overwrite=True)
lwo_table.write('./Disk/lwo comparison.FITS',overwrite=True)

##print '\n'
##print alpha_table
##print '\n'
##print sigma_table
##print '\n'
##print lwo_table
##print '\n'

########################################################################################################################################
#                                                      Galactic Disk GMC Analysis                                                      #
#                                                           <<NOT WORKING>>                                                            #
########################################################################################################################################

# Concatenate all the galactic disk data.
inclination = (0,0,0,0,0,0,0,0,0,0)
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

# Histogram of the parameters.

parameter_hist = plt.figure(figsize=(18,6))

#alpha
ax = plt.subplot(1,3,1)
n,bins,patches = plt.hist(alpha_disk,bins=50,normed=True,alpha=0.6,facecolor='green')
alpha_mu = alpha_disk.mean()
alpha_std = alpha_disk.std()
y = mlab.normpdf(bins,alpha_mu,alpha_std)
plt.plot(bins,y,'k--')
mpl.rc('xtick',labelsize=18)
mpl.rc('ytick',labelsize=18)
plt.text(0.7,0.95,'$\mu$ = %5.4f\n$\sigma$ = %5.4f'%(alpha_mu,alpha_std),ha='left',va='top',transform=ax.transAxes)
plt.title('Galactic Disk Virial Parameter')
plt.xlabel(r'$log_\mathrm{10}(\alpha)$',fontsize=20)

#sigma
ax = plt.subplot(1,3,2)
n,bins,patches = plt.hist(sigma_disk,bins=50,normed=True,alpha=0.6,facecolor='green')
sigma_mu = sigma_disk.mean()
sigma_std = sigma_disk.std()
y = mlab.normpdf(bins,sigma_mu,sigma_std)
plt.plot(bins,y,'k--')
mpl.rc('xtick',labelsize=18)
mpl.rc('ytick',labelsize=18)
plt.text(0.7,0.95,'$\mu$ = %5.4f\n$\sigma$ = %5.4f'%(sigma_mu,sigma_std),ha='left',va='top',transform=ax.transAxes)
plt.title('Galactic Disk Surface density')
plt.xlabel(r'$log_\mathrm{10}(\Sigma)$',fontsize=20)

#lwo
ax = plt.subplot(1,3,3)
n,bins,patches = plt.hist(lwo_disk,bins=50,normed=True,alpha=0.6,facecolor='green')
lwo_mu = lwo_disk.mean()
lwo_std = lwo_disk.std()
y = mlab.normpdf(bins,lwo_mu,lwo_std)
plt.plot(bins,y,'k--')
mpl.rc('xtick',labelsize=18)
mpl.rc('ytick',labelsize=18)
plt.text(0.7,0.95,'$\mu$ = %5.4f\n$\sigma$ = %5.4f'%(lwo_mu,lwo_std),ha='left',va='top',transform=ax.transAxes)
plt.title('Galactic Disk Normalized Linewidth')
plt.xlabel(r'$\sigma_{o}$',fontsize=20)

plt.tight_layout()
plt.savefig('./Disk/parameter_hist.png')
plt.close()

# Star formation rate and depletion time histograms.

SFR = [0]
depletion = [0]
T = 20
mu = 0.71
##m_H = 1.672621898   #*10**-27
##k = 1.38064852      #*10**-23
##G = 6.67408         #*10**-11
alpha_disk = 10**alpha_disk

for i in range(len(mass_disk)):
     mach = 207.27*(lw_disk[i]/100)*((mu/0.71*20/T)**(0.5))
     t_ff = (5.02406*10**11)*(rad_disk[i]**(1.5))*((mass_disk[i])**(-0.5))
     rate_form = 0.014*mass_disk[i]/t_ff*((alpha_disk[i]/1.3)**(-0.68))*((mach/100)**(-0.32))
     SFR = np.append(SFR,rate_form)
     rate_depl = mass_disk[i]/rate_form
     depletion = np.append(SFR,rate_depl)
     print i,rate_form,rate_depl
SFR = np.delete(SFR,0)
SFR = np.delete(SFR,np.isnan(SFR))
depletion = np.delete(depletion,0)
depletion = np.delete(depletion,np.isnan(depletion))
print np.nanmax(SFR),np.nanmin(SFR),np.nanmax(depletion),np.nanmin(depletion)

formation_hist = plt.figure(figsize=(12,6))

ax = plt.subplot(1,2,1)
n,bins,patches = plt.hist(SFR,bins=50,normed=True,log=True,alpha=0.6,facecolor='blue')
SFR_mu = SFR.mean()
SFR_std = SFR.std()
y = mlab.normpdf(bins,SFR_mu,SFR_std)
plt.plot(bins,y,'r--')
mpl.rc('xtick',labelsize=18)
mpl.rc('ytick',labelsize=18)
plt.text(0.5,0.95,'$\mu$ = %5.4f\n$\sigma$ = %5.4f'%(SFR_mu,SFR_std),ha='left',va='top',transform=ax.transAxes)
plt.title('Galactic Disk Star Formation Rate')
plt.xlabel(r'$SFR (s)$',fontsize=20)

ax = plt.subplot(1,2,2)
n,bins,patches = plt.hist(depletion,bins=50,normed=True,log=True,alpha=0.6,facecolor='blue')
depletion_mu = depletion.mean()
depletion_std = depletion.std()
y = mlab.normpdf(bins,depletion_mu,depletion_std)
plt.plot(bins,y,'r--')
mpl.rc('xtick',labelsize=18)
mpl.rc('ytick',labelsize=18)
plt.text(0.5,0.95,'$\mu$ = %5.4f\n$\sigma$ = %5.4f'%(depletion_mu,depletion_std),ha='left',va='top',transform=ax.transAxes)
plt.title('Galactic Disk Depletion Rate')
plt.xlabel(r'$t_d (s)$',fontsize=20)

plt.tight_layout()
plt.savefig('./Disk/formation_hist.png')
plt.close()

# Create a table of some galaxy properties.

properties = ['Galaxy','Distance (Mpc)','Inclination (deg)','Number of disk GMCs','R_nuc (kpc)']
column_types = ['S7','f4','f4','f4','f4']
galaxy_table = Table(names=properties,dtype=column_types)

for i in range(len(galaxy_names)):
    galaxy_table.add_row()
    galaxy_table[-1]['Galaxy'] = galaxy_names[i]
    galaxy_table[-1]['Distance (Mpc)'] = distance[i]
    galaxy_table[-1]['Inclination (deg)'] = inclination[i]
    galaxy_table[-1]['Number of disk GMCs'] = len(mass_dict[i])
    galaxy_table[-1]['R_nuc (kpc)'] = r_nuc[i]

galaxy_table.write('./Disk/galaxy_table.FITS',overwrite=True)
##print galaxy_table

# Fit the data in the power law and truncated-power-law distributions.
myfit = powerlaw.Fit(mass_disk)
R, p = myfit.distribution_compare('power_law','truncated_power_law')
fig = myfit.truncated_power_law.plot_ccdf(label='Truncated\nPower Law')
myfit.power_law.plot_ccdf(label='Power Law',ax=fig)
myfit.plot_ccdf(drawstyle='steps',label='Data',ax=fig)

# Format the plot.
plt.legend(loc=0)
plt.title('Galactic Disk GMC Mass Distribution')
plt.xlabel(r'$M_\mathrm{\odot}$')
plt.ylabel('CCDF')
plt.text(0.45,0.05,'$R$ = %5.4f\n$p$ = %5.4f'%(R,p),ha='left',va='bottom',transform=fig.transAxes)
plt.savefig('./Disk/disk_power_law_mass.png')
plt.close()
