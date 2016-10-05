def graphs (galaxyname,data,ofile):
# This is a module to save data in a series of graphs.

     # import the libraries that are needed.
     from astropy.table import Table
     import matplotlib.pyplot as plt

     # Load all of the data.
     datatable = Table.read(data)

     # Extract different data sets.
     vmass = datatable['VIRMASS_EXTRAP_DECONV']
     lmass = datatable['MASS_EXTRAP']
     radius = datatable['RADRMS_EXTRAP_DECONV']
     lw = datatable['VRMS_EXTRAP']

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
     plt.setp (l1,marker='s',markersize=5,linestyle='None')
     plt.setp (l2,color='k',linestyle='solid',linewidth=3,alpha=0.5)
     plt.legend ([l1,l2],[galaxyname,r'$M_{\mathrm{lum}}$'],loc=2)
     plt.xlabel(r'$M_{\mathrm{lum}}\ (M_{\odot})$',fontsize=20)
     plt.ylabel(r'$M_{\mathrm{vir}}\ (M_{\odot})$',fontsize=20)
     plt.title('Luminous Mass vs. Virial Mass for GMCs in '+galaxyname,fontsize=16)

     # Create a subplot for the GMC luminous radius and mass.
     plt.subplot(1,3,2)
     lines = plt.loglog(radius,lmass,radius,mcan)
     l1,l2 = lines
     plt.setp(l1,marker='s',markersize=5,linestyle='None')
     plt.setp(l2,color='k',linestyle='solid',linewidth=3,alpha=0.5)
     plt.legend ([l1,l2],[galaxyname,'Canonical Mass'],loc=2)
     plt.xlabel(r'$R\ (\mathrm{pc})$',fontsize=20)
     plt.ylabel(r'$M_{\mathrm{lum}}\ (M_{\odot})$',fontsize=20)
     plt.title('Radius vs. Luminous Mass for GMCs in '+galaxyname,fontsize=16)

     # Create a subplot for the GMC radius and linewidth.
     plt.subplot(1,3,3)
     lines = plt.loglog(radius,lw,radius,lwcan)
     l1,l2 = lines
     plt.setp(l1,marker='s',markersize=5,linestyle='None')
     plt.setp(l2,color='k',linestyle='solid',linewidth=3,alpha=0.5)
     plt.legend ([l1,l2],[galaxyname,'Canonical Linewidth'],loc=2)
     plt.xlabel(r'$R\ (\mathrm{pc})$',fontsize=20)
     plt.ylabel(r'$\sigma_{v}\ (\frac{\mathrm{km}}{\mathrm{s}})$',fontsize=20)
     plt.title('Radius vs. Linewidth for GMCs in '+galaxyname,fontsize=16)

     # Reduce the size of the labels and save the figure
     plt.tight_layout()
     plt.savefig(ofile)

def mapgmc (galaxyname,data,imgfile,ofile):
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
     imgheader['NAXIS']=3
     imgheader['NAXIS3']=1

     # Load all of the data.
     datatable = Table.read(data)
     img = aplpy.FITSFigure(imgfile,slices=[4,5])

     # Extract data for the GMC locations in the galaxy.
     xval = datatable['XPOS']
     yval = datatable['YPOS']

     # Display the image of the galaxy with markers showing the locations of the GMCs.
     img.show_colorscale(cmap='gist_heat')
     img.show_markers(xval,yval,edgecolor='white',facecolor='none',marker='D',s=10,alpha=0.7)

     # Save the figure.
     img.save(ofile)

def param (galaxyname,data,ofile1,ofile2):
# This is a module that will plot graphs of the clouds' virial parameter, surface density,
# and linewidth on a one-parsec scale to its galactocentric radius.

     # import the libraries that are needed.
     from astropy.table import Table
     import astropy.units as u
     from galaxies import Galaxy
     import matplotlib.pyplot as plt
     from scipy import stats as sps
     import numpy as np

     # Load all of the data.
     cpropstable = Table.read(data)
     gxy = Galaxy(galaxyname)

     # Extract different data set
     vmass = cpropstable['VIRMASS_EXTRAP_DECONV']
     lmass = cpropstable['MASS_EXTRAP']
     radrms = cpropstable['RADRMS_EXTRAP_DECONV']
     lw = cpropstable['VRMS_EXTRAP']
     xval = cpropstable['XPOS']
     yval = cpropstable['YPOS']
     rad = gxy.radius(ra=xval,dec=yval)
     rad = rad.to(u.kpc)
     print np.amax(rad)

     # Calculate the virial parameter, surface density, and linewidth on a one-parsec scale.
     alpha = vmass/lmass
     sigma = lmass/(radrms**2)
     lwo = lw/((radrms*1000)**0.5)

     # Separate the data into bins, then calculate the mean values.
     binaxis = [0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5]
     bin_means,bin_edges,binnumber = sps.binned_statistic(rad,lmass,statistic='mean',bins=[0,1,2,3,4,5,6,7,8,9,10,11])
     bin_medians,bin_edges,binnumber = sps.binned_statistic(rad,lmass,statistic='median',bins=[0,1,2,3,4,5,6,7,8,9,10,11])
     alpha_medians,bin_edges,binnumber = sps.binned_statistic(rad,alpha,statistic=np.nanmedian,bins=[0,1,2,3,4,5,6,7,8,9,10,11])
     sigma_medians,bin_edges,binnumber = sps.binned_statistic(rad,sigma,statistic=np.nanmedian,bins=[0,1,2,3,4,5,6,7,8,9,10,11])
     lwo_medians,bin_edges,binnumber = sps.binned_statistic(rad,lwo,statistic=np.nanmedian,bins=[0,1,2,3,4,5,6,7,8,9,10,11])

     # Plot the parameters with respect to the galactocentric radius.
     figure1 = plt.figure(figsize=(24,8))
     
     plt.subplot(1,3,1)
     l1 = plt.plot (rad,alpha)
     plt.setp (l1,marker='D',markersize=5,linestyle='None',color='b',label='GMC')
     plt.errorbar(binaxis,alpha_medians,xerr=0.5,linestyle='None',marker='o',color='k',lw=5,alpha=0.7,label='Median')
#     plt.legend (loc=2)
     plt.xlabel(r'$R\ (\mathrm{kpc})$',fontsize=20)
     plt.ylabel(r'$\alpha$',fontsize=20)
     plt.title('Virial Parameter for GMCs in '+galaxyname,fontsize=20)
     
     plt.subplot(1,3,2)
     l1 = plt.plot(rad,sigma)
     plt.setp (l1,marker='D',markersize=5,linestyle='None',color='b',label='GMC')
     plt.errorbar(binaxis,sigma_medians,xerr=0.5,linestyle='None',marker='o',color='k',lw=5,alpha=0.7,label='Median')
#     plt.legend (loc=2)
     plt.xlabel(r'$R\ (\mathrm{kpc})$',fontsize=20)
     plt.ylabel(r'$\Sigma\ (\frac{M_\mathrm{\odot}}{\mathrm{pc}^{2}})$',fontsize=20)
     plt.title('Luminous Mass for GMCs in '+galaxyname,fontsize=20)
     
     plt.subplot(1,3,3)
     l1 = plt.plot(rad,lwo)
     plt.setp (l1,marker='D',markersize=5,linestyle='None',color='b',label='GMC')
     plt.errorbar(binaxis,lwo_medians,xerr=0.5,linestyle='None',marker='o',color='k',lw=5,alpha=0.7,label='Median')
#     plt.legend (loc=2)
     plt.xlabel(r'$R\ (\mathrm{kpc})$',fontsize=20)
     plt.ylabel(r'$\sigma_{o}\ (\frac{\mathrm{km}}{\mathrm{s}})$',fontsize=20)
     plt.title('Linewidth for GMCs in '+galaxyname,fontsize=20)
     
     plt.tight_layout()
     plt.savefig(ofile1)

     # Make a figure of the median mass for a given radius.
     figure2 = plt.figure(figsize=(8,8))
     
     plt.errorbar(binaxis,bin_means,xerr=0.5,linestyle='None',marker='D',color='b',label='Mean')
     plt.errorbar(binaxis,bin_medians,xerr=0.5,linestyle='None',marker='D',color='g',label='Median')
     plt.legend (loc=9)
     plt.xlabel(r'$R\ (\mathrm{kpc})$',fontsize=20)
     plt.ylabel(r'$M_{\mathrm{lum}}\ (M_{\odot})$',fontsize=20)
     plt.title('GMC Mass Distribution in '+galaxyname,fontsize=20)
     
     plt.tight_layout()
     plt.savefig(ofile2)
