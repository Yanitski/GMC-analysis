# This is a module to display data in a series of graphs
def graphs(data,ofile):
	from astropy.table import Table
	import matplotlib.pyplot as plt
	mytable = Table.read(data)
	figure = plt.figure(figsize=(16,6))
	# Create a subplot for the GMC virial and luminous masses
	plt.subplot(1,3,1)
	plt.loglog(mytable['VIRMASS_EXTRAP_DECONV'],mytable['MASS_EXTRAP'],marker='s',linestyle='None')
	plt.xlabel(r'$M_{\mathrm{vir}}\ (M_{\odot})$')
	plt.ylabel(r'$M_{\mathrm{lum}}\ (M_{\odot})$') 
	plt.title('Virial Mass v Luminous Mass')
	# Create a subplot for the GMC luminous mass and radius
	plt.subplot(1,3,2)
	plt.loglog(mytable['MASS_EXTRAP'],mytable['RADRMS_EXTRAP_DECONV'],marker='s',linestyle='None')
	plt.xlabel(r'$M_{\mathrm{lum}}\ (M_{\odot})$') 
	plt.ylabel(r'$\mathrm{Radius}\ (\mathrm{pc})$')
	plt.title('Luminous Mass v Radius')
	# Create a subplot for the GMC linewidth and radius
	plt.subplot(1,3,3)
	plt.loglog(mytable['VRMS_EXTRAP'],mytable['RADRMS_EXTRAP_DECONV'],marker='s',linestyle='None')
	plt.xlabel(r'$\mathrm{Linewidth}\ (\frac{km}{s})$') 
	plt.ylabel(r'$\mathrm{Radius}\ (\mathrm{pc})$')
	plt.title('Linewidth v Radius')
	plt.tight_layout() 	
	plt.savefig(ofile)
