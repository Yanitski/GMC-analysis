def graphs (galaxyname,data):
# This is a module to save data in a series of graphs.

     # import the libraries that are needed.
     from astropy.table import Table
     import matplotlib.pyplot as plt

     # Load all of the data.
     datatable = Table.read(data)

     # Extract different data sets.
     lw = datatable['VRMS_EXTRAP']
     vmass = datatable['VIRMASS_EXTRAP_DECONV']
     lmass = datatable['MASS_EXTRAP']
     radius = datatable['RADRMS_EXTRAP_DECONV']

     # Calculate the canonical values according to the Milky Way.
     mcan = 170*radius**2
     lwcan = 0.7*radius**0.5

     # Initialize a figure that can easily show three plots
     #arranged horizontally.
     figure = plt.figure(figsize=(20,6))

     # Create a subplot for the GMC luminous and virial masses.
     plt.subplot(1,3,1)
     lines = plt.loglog (lmass,vmass,lmass,lmass)
     l1,l2 = lines
     plt.setp(l1,marker='s',markersize=5,linestyle='None')
     plt.setp(l2,color='k',linestyle='solid',linewidth=3,alpha=0.5)
     plt.legend([l1,l2],[galaxyname,r'$M_{\mathrm{lum}}$'],loc=2)
     plt.xlabel(r'$M_{\mathrm{lum}}\ (M_{\odot})$',fontsize=20)
     plt.ylabel(r'$M_{\mathrm{vir}}\ (M_{\odot})$',fontsize=20)
     plt.title('(a)',fontsize=16)

     # Create a subplot for the GMC luminous radius and mass.
     plt.subplot(1,3,2)
     lines = plt.loglog(radius,lmass,radius,mcan)
     l1,l2 = lines
     plt.setp(l1,marker='s',markersize=5,linestyle='None')
     plt.setp(l2,color='k',linestyle='solid',linewidth=3,alpha=0.5)
     plt.legend([l1,l2],[galaxyname,'Canonical Mass'],loc=2)
     plt.xlabel(r'$R\ (\mathrm{pc})$',fontsize=20)
     plt.ylabel(r'$M_{\mathrm{lum}}\ (M_{\odot})$',fontsize=20)
     plt.title('(b)',fontsize=16)

     # Create a subplot for the GMC radius and linewidth.
     plt.subplot(1,3,3)
     lines = plt.loglog(radius,lw,radius,lwcan)
     l1,l2 = lines
     plt.setp(l1,marker='s',markersize=5,linestyle='None')
     plt.setp(l2,color='k',linestyle='solid',linewidth=3,alpha=0.5)
     plt.legend([l1,l2],[galaxyname,'Canonical Linewidth'],loc=2)
     plt.xlabel(r'$R\ (\mathrm{pc})$',fontsize=20)
     plt.ylabel(r'$\sigma_{v}\ (\frac{\mathrm{km}}{\mathrm{s}})$',fontsize=20)
     plt.title('(c)',fontsize=16)

     # Reduce the size of the labels and save the figure
     plt.tight_layout()
     plt.savefig(galaxyname+'_trend.png')
     plt.close()

def mapgmc (galaxyname,data,imgfile):
# This is a module to plot the location of the GMCs in a given galaxy on a map of the galaxy.

     # import the libraries that are needed.
     import aplpy
     from astropy.table import Table
     from astropy.io import fits
     import matplotlib.pyplot as plt
     import numpy as np

     # Modify the FITS image file so it is compatible with APLpy
     imgdata = fits.open(imgfile)
     imgheader = imgdata[0].header
     imgheader['NAXIS'] = 3
     imgheader['NAXIS3'] = 1

     # Load all of the data.
     datatable = Table.read(data)
     img = aplpy.FITSFigure(imgfile,slices=[4,5])

     # Extract data for the GMC locations in the galaxy.
     xval = datatable['XPOS']
     yval = datatable['YPOS']

     # Display the image of the galaxy with markers showing the locations of the GMCs.
     img.show_colorscale(cmap='gist_heat')
     img.show_markers(xval,yval,edgecolor='white',facecolor='none',marker='D',s=10,alpha=0.3)

     # Save the figure.
     img.save(galaxyname+'_map.png')

def param (galaxyname,data):
# This is a module that will plot graphs of the clouds' virial parameter, surface density,
# and linewidth on a one-parsec scale to its galactocentric radius.

     # import the libraries that are needed.
     from scipy import stats as sps
     from galaxies import Galaxy
     from astropy.table import Table
     import numpy as np
     import astropy.units as u
     import matplotlib.pyplot as plt

     # Load all of the data.
     gxy = Galaxy(galaxyname)
     cpropstable = Table.read(data)

     # Extract different data set
     lw = cpropstable['VRMS_EXTRAP']
     xval = cpropstable['XPOS']
     yval = cpropstable['YPOS']
     vmass = cpropstable['VIRMASS_EXTRAP_DECONV']
     lmass = cpropstable['MASS_EXTRAP']
     radrms = cpropstable['RADRMS_EXTRAP_DECONV']
     rad = gxy.radius(ra=xval,dec=yval)
     rad = rad.to(u.kpc)

     # Calculate the virial parameter, surface density, and linewidth on a one-parsec scale.
     lwo = lw/((radrms)**0.5)
     alpha = vmass/lmass
     sigma = lmass/(radrms**2)

     # Separate the data into bins, then calculate the mean values.
     binaxis = [0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5]
     bin_means,bin_edges,binnumber = sps.binned_statistic(rad,lmass,statistic='mean',bins=[0,1,2,3,4,5,6,7,8,9,10,11])
     bin_medians,bin_edges,binnumber = sps.binned_statistic(rad,lmass,statistic='median',bins=[0,1,2,3,4,5,6,7,8,9,10,11])
     lwo_medians,bin_edges,binnumber = sps.binned_statistic(rad,lwo,statistic=np.nanmedian,bins=[0,1,2,3,4,5,6,7,8,9,10,11])
     alpha_medians,bin_edges,binnumber = sps.binned_statistic(rad,alpha,statistic=np.nanmedian,bins=[0,1,2,3,4,5,6,7,8,9,10,11])
     sigma_medians,bin_edges,binnumber = sps.binned_statistic(rad,sigma,statistic=np.nanmedian,bins=[0,1,2,3,4,5,6,7,8,9,10,11])

     # Plot the parameters with respect to the galactocentric radius.
     figure1 = plt.figure(figsize=(24,8))
     
     plt.subplot(1,3,1)
     l1 = plt.plot (rad,alpha)
     plt.yscale('log')
     plt.setp (l1,marker='D',markersize=5,linestyle='None',color='b',label='GMC')
     plt.errorbar(binaxis,alpha_medians,xerr=0.5,linestyle='None',marker='o',color='k',lw=5,alpha=0.7,label='Median')
     plt.legend(loc=0)
     plt.xlabel(r'$R\ (\mathrm{kpc})$',fontsize=20)
     plt.ylabel(r'$\alpha$',fontsize=20)
     plt.title('(a)',fontsize=20)
     
     plt.subplot(1,3,2)
     l1 = plt.plot(rad,sigma)
     plt.yscale('log')
     plt.setp (l1,marker='D',markersize=5,linestyle='None',color='b',label='GMC')
     plt.errorbar(binaxis,sigma_medians,xerr=0.5,linestyle='None',marker='o',color='k',lw=5,alpha=0.7,label='Median')
     plt.legend(loc=0)
     plt.xlabel(r'$R\ (\mathrm{kpc})$',fontsize=20)
     plt.ylabel(r'$\Sigma\ (\frac{M_\mathrm{\odot}}{\mathrm{pc}^{2}})$',fontsize=20)
     plt.title('(b)',fontsize=20)
     
     plt.subplot(1,3,3)
     l1 = plt.plot(rad,lwo)
     plt.setp (l1,marker='D',markersize=5,linestyle='None',color='b',label='GMC')
     plt.errorbar(binaxis,lwo_medians,xerr=0.5,linestyle='None',marker='o',color='k',lw=5,alpha=0.7,label='Median')
     plt.legend(loc=0)
     plt.xlabel(r'$R\ (\mathrm{kpc})$',fontsize=20)
     plt.ylabel(r'$\sigma_{o}\ (\frac{\mathrm{km}}{\mathrm{s}})$',fontsize=20)
     plt.title('(c)',fontsize=20)
     
     plt.tight_layout()
     plt.savefig(galaxyname+'_param.png')
     plt.close()

     # Make a figure of the median mass for a given radius.
     figure2 = plt.figure(figsize=(8,8))
     
     plt.errorbar(binaxis,bin_means,xerr=0.5,linestyle='None',marker='D',color='b',label='Mean')
     plt.errorbar(binaxis,bin_medians,xerr=0.5,linestyle='None',marker='D',color='g',label='Median')
     plt.legend(loc=0)
     plt.xlabel(r'$R\ (\mathrm{kpc})$',fontsize=20)
     plt.ylabel(r'$M_{\mathrm{lum}}\ (M_{\odot})$',fontsize=20)
     plt.title('GMC Mass Distribution in '+galaxyname,fontsize=20)
     
     plt.tight_layout()
     plt.savefig(galaxyname+'_dist.png')
     plt.close()

# This is some code from the powerlaw template.
def data(galaxyname,data,n_bins):
     
     #import libraries
     from galaxies import Galaxy
     from astropy.table import Table
     from astropy.table import Column
     import astropy
     import astropy.table
     import powerlaw
     import numpy as np
     import astropy.units as u
     import matplotlib.pyplot as plt

     #get info about m83
     mygalaxy = Galaxy(galaxyname)
     print mygalaxy

     #load fits file
     t = Table.read(data)

     #find cloud's galactocentric distance
     rgal = mygalaxy.radius(ra=(t['XPOS']), dec=(t['YPOS']))
     rgal = rgal.to(u.pc)
     rpgal = np.asarray(rgal)

     #add those distances to the fits table
     col_rgal = Column(name='RADIUS_PC',data=(rgal))
     t.add_column(col_rgal)

     # Sort the masses according to galactocentric radius.
     mass = t['MASS_EXTRAP'].data
     i_sorted = np.argsort(rgal)
     rgal_sorted = np.asarray(rgal[i_sorted])
     mass_sorted = mass[i_sorted]

     # Create a loop to calculate the bin boundaries.
     totmass = np.sum(mass)/n_bins
     edge_f = 1.1*np.max(rgal_sorted)
     edges = np.zeros(n_bins-1)
     start = 0
     e = 0
     for i in range(len(mass_sorted)):
          if np.sum(mass_sorted[start:i])>totmass:
               edges[e] = 0.5*(rgal_sorted[i]+rgal_sorted[i-1])
               start = i
               e = e+1
     inneredge = np.array([0, edges[:]])
     outeredge = np.array([edges[:], edge_f])
     print outeredge, inneredge

     # Create a template for a new FITS table.
     column_names = ['Inner edge (pc)', 'Outer edge (pc)', 'GMC index', 'R', 'p', 'Truncation mass ($M_\mathrm{\odot}$)', 'Largest cloud ($M_\mathrm{\odot}$)', '5th largest cloud ($M_\mathrm{\odot}$)']
     column_types = ['f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4']
     table = Table(names=column_names, dtype=column_types)

     # Fill the FITS table.
     for inneredge, outeredge in zip(inneredge, outeredge):
          
         idx = np.where((t['RADIUS_PC']>=inneredge)&(t['RADIUS_PC']<outeredge))
         mass = t['MASS_EXTRAP'][idx].data
         print len(mass)
         fit = powerlaw.Fit(mass)
         fit_subset = powerlaw.Fit(mass, xmin=3e5)
         R,p = fit_subset.distribution_compare('power_law', 'truncated_power_law')
         table.add_row()
         table[-1]['R'] = R
         table[-1]['p'] = p
         table[-1]['GMC index'] = -fit_subset.alpha
         table[-1]['Inner edge (pc)'] = inneredge
         table[-1]['Outer edge (pc)'] = outeredge
         table[-1]['Largest cloud ($M_\mathrm{\odot}$)'] = np.nanmax(mass)
         table[-1]['Truncation mass ($M_\mathrm{\odot}$)'] = 1/fit_subset.truncated_power_law.parameter2
         table[-1]['5th largest cloud ($M_\mathrm{\odot}$)'] = np.sort(t['MASS_EXTRAP'][idx])[-5]

         #print(table)
         #print(-fit.alpha, -fit_subset.alpha, R, p, 1/fit_subset.truncated_power_law.parameter2)
         
     table.write(galaxyname+'_data.fits', overwrite=True)
     print table

# This is a module to analyze whether a galaxy's GMC distribution is truncated.
def pwrlaw(galaxyname, data):

     from astropy.table import Table
     import powerlaw
     import numpy as np
     import matplotlib.pyplot as plt

     t = Table.read(data)
     n_bins = len(t['Inner Edge (pc)'])
     mass = t['MASS_EXTRAP'].data
     myfit = powerlaw.Fit(mass)
     R, p = myfit.distribution_compare('power_law','truncated_power_law')
     fig = myfit.truncated_power_law.plot_ccdf(label='Truncated\nPower Law')
     myfit.power_law.plot_ccdf(label='Power Law',ax=fig)
     myfit.plot_ccdf(drawstyle='steps',label='Data',ax=fig)
     plt.legend(loc=0)
     plt.title('GMC Mass Distribution')
     plt.xlabel(r'$M_\mathrm{\odot}$')
     plt.ylabel('CCDF')
     plt.savefig(galaxyname+'_power_law.png')
     plt.close()
