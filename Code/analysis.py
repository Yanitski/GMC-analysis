########################################################################################################################################
#                                                           GMC Property Analysis                                                      #
########################################################################################################################################

def graphs (galaxyname,data):
# This is a module to save data in a series of graphs.

     # import the libraries that are needed.
     from astropy.table import Table
     import matplotlib.pyplot as plt
     import matplotlib as mpl

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
     figure = plt.figure(figsize=(15,5))

     # Create a subplot for the GMC luminous and virial masses.
     ax = plt.subplot(1,3,1)
     lines = plt.loglog (lmass,vmass,lmass,lmass)
     l1,l2 = lines
     plt.setp(l1,marker='s',markersize=5,linestyle='None')
     plt.setp(l2,color='k',linestyle='solid',linewidth=3,alpha=0.5)
     plt.legend([l1,l2],[galaxyname,r'$M_{\mathrm{lum}}$'],loc=2)
     mpl.rc('xtick',labelsize=16)
     mpl.rc('ytick',labelsize=16)
     plt.text(0.01,0.5,'(a)',ha='left',va='center',transform=ax.transAxes,fontsize=16)
     plt.xlabel(r'$M_{\mathrm{lum}}\ (M_{\odot})$',fontsize=20)
     plt.ylabel(r'$M_{\mathrm{vir}}\ (M_{\odot})$',fontsize=20)
     #plt.title('(a)',fontsize=20)

     # Create a subplot for the GMC luminous radius and mass.
     ax = plt.subplot(1,3,2)
     lines = plt.loglog(radius,lmass,radius,mcan)
     l1,l2 = lines
     plt.setp(l1,marker='s',markersize=5,linestyle='None')
     plt.setp(l2,color='k',linestyle='solid',linewidth=3,alpha=0.5)
     plt.legend([l1,l2],[galaxyname,'Canonical Mass'],loc=2)
     mpl.rc('xtick',labelsize=16)
     mpl.rc('ytick',labelsize=16)
     plt.text(0.01,0.5,'(a)',ha='left',va='center',transform=ax.transAxes,fontsize=16)
     plt.xlabel(r'$R\ (\mathrm{pc})$',fontsize=20)
     plt.ylabel(r'$M_{\mathrm{lum}}\ (M_{\odot})$',fontsize=20)
     #plt.title('(b)',fontsize=20)

     # Create a subplot for the GMC radius and linewidth.
     ax = plt.subplot(1,3,3)
     lines = plt.loglog(radius,lw,radius,lwcan)
     l1,l2 = lines
     plt.setp(l1,marker='s',markersize=5,linestyle='None')
     plt.setp(l2,color='k',linestyle='solid',linewidth=3,alpha=0.5)
     plt.legend([l1,l2],[galaxyname,'Canonical Linewidth'],loc=2)
     mpl.rc('xtick',labelsize=16)
     mpl.rc('ytick',labelsize=16)
     plt.text(0.01,0.5,'(a)',ha='left',va='center',transform=ax.transAxes,fontsize=16)
     plt.xlabel(r'$R\ (\mathrm{pc})$',fontsize=20)
     plt.ylabel(r'$\sigma_{v}\ (\frac{\mathrm{km}}{\mathrm{s}})$',fontsize=20)
     #plt.title('(c)',fontsize=20)

     # Reduce the size of the labels and save the figure
     plt.tight_layout()
     plt.savefig('../Graphs/'+galaxyname+'_trend.png')
     plt.close()

########################################################################################################################################
#                                                           Galaxy GMC Mapping                                                         #
########################################################################################################################################

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
     img.show_markers(xval,yval,edgecolor='cyan',marker='D',s=10,alpha=0.7)

     # Save the figure.
     img.save('../Maps/'+galaxyname+'_map.png')

########################################################################################################################################
#                                                           GMC Parameter Analysis                                                     #
########################################################################################################################################)

def param (galaxyname,data,r_nuc=0):
# This is a module that will plot graphs of the clouds' virial parameter, surface density,
# and linewidth on a one-parsec scale to its galactocentric radius.

     # import the libraries that are needed.
     from scipy import stats as sps
     from galaxies import Galaxy
     from astropy.table import Table
     import numpy as np
     import astropy.units as u
     import matplotlib.pyplot as plt
     import matplotlib as mpl

     # Load all of the data.
     gxy = Galaxy(galaxyname)
     cpropstable = Table.read(data)

     # Extract different data set
     lw = cpropstable['VRMS_EXTRAP']
     xval = cpropstable['XPOS']
     yval = cpropstable['YPOS']
     vmass = cpropstable['VIRMASS_EXTRAP_DECONV']
     mass = cpropstable['MASS_EXTRAP']
     radrms = cpropstable['RADRMS_EXTRAP_DECONV']
     rgal = gxy.radius(ra=xval,dec=yval)
     rgal = rgal.to(u.kpc)

     # Sort the masses according to galactocentric radius.
     i_sorted = np.argsort(rgal)
     rgal_sorted = np.asarray(rgal[i_sorted])
     mass_sorted = np.asarray(mass[i_sorted])

     # Calculate the virial parameter, surface density, and linewidth on a one-parsec scale.
     alpha = vmass/mass
     sigma = mass/(radrms**2)
     lwo = lw/((radrms)**0.5)

     # Initialize the arrays to extract the glactic disk cloud properties.
     rad_sorted = radrms[i_sorted]
     alpha_sorted = np.log10(alpha[i_sorted])
     sigma_sorted = np.log10(sigma[i_sorted])
     lwo_sorted = np.log10(lwo[i_sorted])
     #lwo_sorted = lwo[i_sorted]
     lw_sorted = lw[i_sorted]
     mass_disk = []      #mass of each disk GMC
     rad_disk = []       #radius of each disk GMC
     alpha_disk = []     #virial parameter of each disk GMC
     sigma_disk = []     #suface density of each disk GMC
     lwo_disk = []       #normalized linewidth of each disk GMC
     lw_disk = []        #linewidth of each disk GMC

     # Create a loop to extract the galactic disk GMC properties.
     for i in range(len(mass_sorted)):
          if rgal_sorted[i]>r_nuc:
               if not np.isnan(alpha_sorted[i]) and not np.isnan(sigma_sorted[i]) and not np.isnan(lwo_sorted[i]) and not np.isnan(lw_sorted[i]):
                    mass_disk = np.append(mass_disk,mass_sorted[i])
                    rad_disk = np.append(rad_disk,rad_sorted[i])
                    alpha_disk = np.append(alpha_disk,alpha_sorted[i])
                    sigma_disk = np.append(sigma_disk,sigma_sorted[i])
                    lwo_disk = np.append(lwo_disk,lwo_sorted[i])
                    lw_disk = np.append(lw_disk,lw_sorted[i])

     # Separate the data into bins, then calculate the mean values.
     binaxis = [0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5]
     bin_means,bin_edges,binnumber = sps.binned_statistic(rgal,mass,statistic='mean',bins=[0,1,2,3,4,5,6,7,8,9,10,11])
     bin_medians,bin_edges,binnumber = sps.binned_statistic(rgal,mass,statistic='median',bins=[0,1,2,3,4,5,6,7,8,9,10,11])
     lwo_medians,bin_edges,binnumber = sps.binned_statistic(rgal,lwo,statistic=np.nanmedian,bins=[0,1,2,3,4,5,6,7,8,9,10,11])
     alpha_medians,bin_edges,binnumber = sps.binned_statistic(rgal,alpha,statistic=np.nanmedian,bins=[0,1,2,3,4,5,6,7,8,9,10,11])
     sigma_medians,bin_edges,binnumber = sps.binned_statistic(rgal,sigma,statistic=np.nanmedian,bins=[0,1,2,3,4,5,6,7,8,9,10,11])

     # Plot the parameters with respect to the galactocentric radius.
     figure1 = plt.figure(figsize=(15,5))
     
     ax = plt.subplot(1,3,1)
     l1 = plt.plot (rgal,alpha)
     plt.yscale('log')
     plt.setp (l1,marker='D',markersize=5,linestyle='None',color='b',label='GMC')
     plt.errorbar(binaxis,alpha_medians,xerr=0.5,linestyle='None',marker='o',color='k',lw=5,alpha=0.7,label='Median')
     plt.legend(loc=0)
     mpl.rc('xtick',labelsize=16)
     mpl.rc('ytick',labelsize=16)
     plt.text(0.5,0.99,'(a)',ha='center',va='top',transform=ax.transAxes,fontsize=16)
     plt.xlabel(r'$R_g\ (\mathrm{kpc})$',fontsize=20)
     plt.ylabel(r'$\alpha$',fontsize=20)
     #plt.title('(a)',fontsize=20)
     
     ax = plt.subplot(1,3,2)
     l1 = plt.plot(rgal,sigma)
     plt.yscale('log')
     plt.setp (l1,marker='D',markersize=5,linestyle='None',color='b',label='GMC')
     plt.errorbar(binaxis,sigma_medians,xerr=0.5,linestyle='None',marker='o',color='k',lw=5,alpha=0.7,label='Median')
     plt.legend(loc=0)
     mpl.rc('xtick',labelsize=16)
     mpl.rc('ytick',labelsize=16)
     plt.text(0.5,0.99,'(b)',ha='center',va='top',transform=ax.transAxes,fontsize=16)
     plt.xlabel(r'$R_g\ (\mathrm{kpc})$',fontsize=20)
     plt.ylabel(r'$\Sigma\ (\frac{M_\mathrm{\odot}}{\mathrm{pc}^{2}})$',fontsize=20)
     #plt.title('(b)',fontsize=20)
     
     ax = plt.subplot(1,3,3)
     l1 = plt.plot(rgal,lwo)
     plt.setp (l1,marker='D',markersize=5,linestyle='None',color='b',label='GMC')
     plt.errorbar(binaxis,lwo_medians,xerr=0.5,linestyle='None',marker='o',color='k',lw=5,alpha=0.7,label='Median')
     plt.legend(loc=0)
     mpl.rc('xtick',labelsize=16)
     mpl.rc('ytick',labelsize=16)
     plt.text(0.5,0.99,'(c)',ha='center',va='top',transform=ax.transAxes,fontsize=16)
     plt.xlabel(r'$R_g\ (\mathrm{kpc})$',fontsize=20)
     plt.ylabel(r'$\sigma_{o}\ (\frac{\mathrm{km}}{\mathrm{s}})$',fontsize=20)
     #plt.title('(c)',fontsize=20)

     plt.tight_layout()
     plt.savefig('../Parameters/'+galaxyname+'_param.png')
     plt.close()

     # Make a figure of the median mass for a given radius.
     figure2 = plt.figure(figsize=(8,8))
     
     plt.errorbar(binaxis,bin_means,xerr=0.5,linestyle='None',marker='D',color='b',label='Mean')
     plt.errorbar(binaxis,bin_medians,xerr=0.5,linestyle='None',marker='D',color='g',label='Median')
     plt.legend(loc=0)
     mpl.rc('xtick',labelsize=16)
     mpl.rc('ytick',labelsize=16)
     plt.xlabel(r'$R\ (\mathrm{kpc})$',fontsize=20)
     plt.ylabel(r'$M_{\mathrm{lum}}\ (M_{\odot})$',fontsize=20)
     #plt.title('GMC Mass Distribution in '+galaxyname,fontsize=20)
     
     plt.tight_layout()
     plt.savefig('../Parameters/'+galaxyname+'_dist.png')
     plt.close()

     return [mass_disk,rad_disk,alpha_disk,sigma_disk,lwo_disk,lw_disk]

########################################################################################################################################
#                                                           GMC Power-law Fit                                                          #
########################################################################################################################################

# This is some code from the powerlaw template.
def data (galaxyname,data,n_bins=1,r_nuc=0):

     # Import the libraries.
     from galaxies import Galaxy
     from astropy.table import Table
     from astropy.table import Column
     import astropy
     import powerlaw
     import numpy as np
     import astropy.table
     import astropy.units as u
     import matplotlib.pyplot as plt
     import matplotlib as mpl

     # Load its FITS file.
     t = Table.read(data)

     # Load the information about the galaxy.
     gxy = Galaxy(galaxyname)

     # Calculate the galaxy's properties.
     distance = np.asarray(gxy.distance)
     inclination = np.asarray(gxy.inclination)
     rgal = gxy.radius(ra=(t['XPOS']), dec=(t['YPOS']))
     rgal = rgal.to(u.kpc)
     rpgal = np.asarray(rgal)

     # Append these to the FITS table.
     col_rgal = Column(name='RADIUS_KPC',data=(rgal))
     t.add_column(col_rgal)

     # Sort the masses according to galactocentric radius.
     mass = t['MASS_EXTRAP'].data
     i_sorted = np.argsort(rgal)
     rgal_sorted = np.asarray(rgal[i_sorted])
     mass_sorted = np.asarray(mass[i_sorted])

     # Initiate a loop to calculate the bin boundaries and the indeces of these boundaries in the sorted list.
     totmass = np.sum(mass)/n_bins
     edge_f = 1.1*np.max(rgal_sorted)
     edges = np.zeros(n_bins-1)
     start = 0
     mass_equiv = [0]    #indeces for the sorted mass bins of equal mass
     mass_area = [0]     #indeces for the sorted mass bins of equal area
     rgal_equiv = [0,2,8**0.5,12**0.5,4,20**0.5,24**0.5,28**0.5,32**0.5,6]
     r = 1               #equal-area radial index
     e = 0               #edge index
     f = 0               #loop flag to skip the mass_area loop
     c = 0               #loop counter
     for i in range(len(mass_sorted)):
          #Find the indeces for bins of equal mass (equivalent to totmass)
          if np.sum(mass_sorted[start:i])>totmass:
               edges[e] = 0.5*(rgal_sorted[i]+rgal_sorted[i-1])
               start = i
               mass_equiv = np.append(mass_equiv,i)
               e = e+1
          #Find the indeces for bins of equal area (4pi kpc^2)
          if rgal_sorted[i]>rgal_equiv[r] and not f:
               if rgal_sorted[i]<rgal_equiv[3] and not f:     
                    mass_area = np.append(mass_area,i)
                    r = r+1
                    c = 0
               if rgal_sorted[i]>rgal_equiv[3]:
                    f=1
                    mass_area = np.append(mass_area,i)
     mass_equiv = np.append(mass_equiv,i)
     inneredge = np.concatenate(([0.000000],edges))
     outeredge = np.concatenate((edges,[edge_f]))

     # Create a template for a new table.
     column_names = ['Inner edge (kpc)', 'Outer edge (kpc)', 'GMC index', 'R', 'p', 'Truncation mass ($M_\mathrm{\odot}$)', 'Largest cloud ($M_\mathrm{\odot}$)', '5th largest cloud ($M_\mathrm{\odot}$)']
     column_types = ['f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4']
     table = Table(names=column_names, dtype=column_types)

     # Fill the table.
     for inneredge, outeredge in zip(inneredge, outeredge):
          
         idx = np.where((t['RADIUS_KPC']>=inneredge)&(t['RADIUS_KPC']<outeredge))
         mass = t['MASS_EXTRAP'][idx].data
         fit = powerlaw.Fit(mass)
         fit_subset = powerlaw.Fit(mass, xmin=3e5)
         R,p = fit_subset.distribution_compare('power_law', 'truncated_power_law')
         table.add_row()
         table[-1]['R'] = R
         table[-1]['p'] = p
         table[-1]['GMC index'] = -fit_subset.alpha
         table[-1]['Inner edge (kpc)'] = inneredge
         table[-1]['Outer edge (kpc)'] = outeredge
         table[-1]['Largest cloud ($M_\mathrm{\odot}$)'] = np.nanmax(mass)
         table[-1]['Truncation mass ($M_\mathrm{\odot}$)'] = 1/fit_subset.truncated_power_law.parameter2
         table[-1]['5th largest cloud ($M_\mathrm{\odot}$)'] = np.sort(t['MASS_EXTRAP'][idx])[-5]

     # Write the data to a FITS file.
     table.write('../Data/'+galaxyname+'_data.fits', overwrite=True)

     # Plot the mass distribution trends for equal-mass bins.
     t = Table.read('../Data/'+galaxyname+'_data.fits')
     inneredge = t['Inner edge (kpc)'].data
     outeredge = t['Outer edge (kpc)'].data
     subplot_label = ('(a)','(b)','(c)','(d)','(e)','(f)')
     for i in range(len(mass_equiv)-1):

          binmass = mass_sorted[mass_equiv[i]:mass_equiv[i+1]]
          myfit = powerlaw.Fit(binmass)
          R, p = myfit.distribution_compare('power_law','truncated_power_law')
          fig = myfit.truncated_power_law.plot_ccdf(label='Truncated\nPower Law')
          myfit.power_law.plot_ccdf(label='Power Law',ax=fig)
          myfit.plot_ccdf(drawstyle='steps',label='Data',ax=fig)

          # Format the plot.
          plt.legend(loc=0)
          #plt.title(galaxyname+'Equal-mass Mass Distribution, bin '+repr(i+1))
          plt.ylim(ymin=10**-3)
          plt.xlabel(r'$M_\mathrm{\odot}$',fontsize=20)
          plt.ylabel('CCDF',fontsize=20)
          mpl.rc('xtick',labelsize=16)
          mpl.rc('ytick',labelsize=16)
          plt.text(0.01,0.5,subplot_label[i],ha='left',va='center',transform=fig.transAxes,fontsize=16)
          plt.text(0.35,0.01,r'$M_{bin} = %e M_\mathrm{\odot}$'%(totmass)+'\n'+r'$R_{gal} = %5.4f kpc to %5.4f kpc$'%(inneredge[i],outeredge[i])+
                   '\n'+r'$R = %5.4f, p = %5.4f$'%(R,p)+'\n'+r'$\alpha = %5.4f, M_\mathrm{o} = %5.4eM_\mathrm{\odot}$'%(myfit.alpha,myfit.xmin),
                   ha='left',va='bottom',transform=fig.transAxes,fontsize=16)
          plt.savefig('../Data/'+galaxyname+'_power_law_equal_mass_'+repr(i+1)+'.png')
          plt.close()

     # Plot the mass distribution trend for equal-area bins.
     for i in range(len(mass_area)-1):

          f = 0     #loop flag
          if mass_area[i+1]-mass_area[i]<3:
               f = 1
          if not f:
               binmass = mass_sorted[mass_area[i]:mass_area[i+1]]
               myfit = powerlaw.Fit(binmass)
               R, p = myfit.distribution_compare('power_law','truncated_power_law')
               fig = myfit.truncated_power_law.plot_ccdf(label='Truncated\nPower Law')
               myfit.power_law.plot_ccdf(label='Power Law',ax=fig)
               myfit.plot_ccdf(drawstyle='steps',label='Data',ax=fig)

               # Format the plot.
               plt.legend(loc=0)
               #plt.title(galaxyname+'Equal-area Mass Distribution, bin '+repr(i+1))
               plt.ylim(ymin=10**-3)
               plt.xlabel(r'$M_\mathrm{\odot}$',fontsize=20)
               plt.ylabel('CCDF',fontsize=20)
               mpl.rc('xtick',labelsize=16)
               mpl.rc('ytick',labelsize=16)
               plt.text(0.01,0.5,subplot_label[i],ha='left',va='center',transform=fig.transAxes,fontsize=16)
               plt.text(0.35,0.01,r'$R_{gal} = %5.4f kpc to %5.4f kpc$'%(rgal_equiv[i],rgal_equiv[i+1])+'\n'+r'$R = %5.4f, p = %5.4f$'%(R,p)+'\n'+
                        r'$\alpha = %5.4f, M_\mathrm{o} = %5.4eM_\mathrm{\odot}$'%(myfit.alpha,myfit.xmin),ha='left',va='bottom',transform=fig.transAxes,fontsize=16)
               plt.savefig('../Data/'+galaxyname+'_power_law_equal_area_'+repr(i+1)+'.png')
               plt.close()

     return [distance,inclination]
