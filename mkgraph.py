# This is a module to save data in a series of graphs.
def graphs (galaxy,data,ofile):

        # import the libraries that are needed.
	from astropy.table import Table
	import matplotlib.pyplot as plt

	# Load all of the data.
	mytable = Table.read(data)

	# Extract different data sets.
	vmass = mytable['VIRMASS_EXTRAP_DECONV']
	lmass = mytable['MASS_EXTRAP']
	radius = mytable['RADRMS_EXTRAP_DECONV']
	lw = mytable['VRMS_EXTRAP']

	# Calculate the canonical values according to the Milky Way.
	mcan = 170*radius**2
	lwcan = 0.7*radius**0.5

	# Initialize a figure that can easily show three plots
        #arranged horizontally.
	figure = plt.figure(figsize=(16,6))

	# Create a subplot for the GMC luminous and virial masses.
	plt.subplot(1,3,1)
	lines = plt.loglog (lmass,vmass,lmass,lmass)
	l1,l2 = lines
	plt.setp (l1,marker='s',markersize=5,linestyle='None')
	plt.setp (l2,color='k',linestyle='solid',linewidth=3,alpha=0.5)
	plt.legend ([l1,l2],[galaxy,r'$M_{\mathrm{lum}}$'],loc=2)
	plt.xlabel(r'$M_{\mathrm{lum}}\ (M_{\odot})$')
	plt.ylabel(r'$M_{\mathrm{vir}}\ (M_{\odot})$')
	plt.title('Luminous Mass vs. Virial Mass for GMCs in '+galaxy)

	# Create a subplot for the GMC luminous radius and mass.
	plt.subplot(1,3,2)
	lines = plt.loglog(radius,lmass,radius,mcan)
	l1,l2 = lines
	plt.setp(l1,marker='s',markersize=5,linestyle='None')
	plt.setp(l2,color='k',linestyle='solid',linewidth=3,alpha=0.5)
	plt.legend ([l1,l2],[galaxy,'Canonical Mass'],loc=2)
	plt.xlabel(r'$\mathrm{Radius}\ (\mathrm{pc})$')
	plt.ylabel(r'$M_{\mathrm{lum}}\ (M_{\odot})$')
	plt.title('Radius vs. Luminous Mass for GMCs in '+galaxy)

	# Create a subplot for the GMC radius and linewidth.
	plt.subplot(1,3,3)
	lines = plt.loglog(radius,lw,radius,lwcan)
        l1,l2 = lines
	plt.setp(l1,marker='s',markersize=5,linestyle='None')
	plt.setp(l2,color='k',linestyle='solid',linewidth=3,alpha=0.5)
	plt.legend ([l1,l2],[galaxy,'Canonical Linewidth'],loc=2)
	plt.xlabel(r'$\mathrm{Radius}\ (\mathrm{pc})$')
	plt.ylabel(r'$\mathrm{Linewidth}\ (\mathrm{km})$')
	plt.title('Radius vs. Linewidth for GMCs in '+galaxy)


	plt.tight_layout()
	plt.savefig(ofile)

