import mkgraph
import os

# Create directories for the images if they don't already exist
if not os.path.exists('../Graphs'):
     os.mkdir('../Graphs')
if not os.path.exists('../Maps'):
     os.mkdir('../Maps')
if not os.path.exists('../Parameters'):
     os.mkdir('../Parameters')
if not os.path.exists('../Data'):
     os.mkdir('../Data')

# Load the data file and display it in a few graphs.
#NGC0628
mkgraph.graphs ('NGC0628','NGC0628_co21_cube_pbcor_props_clfind.fits')
mkgraph.mapgmc ('NGC0628','NGC0628_co21_cube_pbcor_props_clfind.fits','NGC0628_co21_cube_tpeak.fits')
mkgraph.param ('NGC0628','NGC0628_co21_cube_pbcor_props_clfind.fits')
mkgraph.data('NGC0628','NGC0628_co21_cube_pbcor_props_clfind.fits',5)

#NGC1672
mkgraph.graphs ('NGC1672','NGC1672_co21_cube_pbcor_props_clfind.fits')
mkgraph.mapgmc ('NGC1672','NGC1672_co21_cube_pbcor_props_clfind.fits','NGC1672_co21_cube_tpeak.fits')
mkgraph.param ('NGC1672','NGC1672_co21_cube_pbcor_props_clfind.fits')
mkgraph.data('NGC1672','NGC1672_co21_cube_pbcor_props_clfind.fits',3)

#NGC3351
mkgraph.graphs ('NGC3351','NGC3351_co21_cube_pbcor_props_clfind.fits')
mkgraph.mapgmc ('NGC3351','NGC3351_co21_cube_pbcor_props_clfind.fits','NGC3351_co21_cube_tpeak.fits')
mkgraph.param ('NGC3351','NGC3351_co21_cube_pbcor_props_clfind.fits')
mkgraph.data('NGC3351','NGC3351_co21_cube_pbcor_props_clfind.fits',3)

#NGC3627
mkgraph.graphs ('NGC3627','NGC3627_co21_cube_pbcor_props_clfind.fits')
mkgraph.mapgmc ('NGC3627','NGC3627_co21_cube_pbcor_props_clfind.fits','NGC3627_co21_cube_tpeak.fits')
mkgraph.param ('NGC3627','NGC3627_co21_cube_pbcor_props_clfind.fits')
mkgraph.data('NGC3627','NGC3627_co21_cube_pbcor_props_clfind.fits',5)

#NGC4254
mkgraph.graphs ('NGC4254','NGC4254_co21_cube_pbcor_props_clfind.fits')
mkgraph.mapgmc ('NGC4254','NGC4254_co21_cube_pbcor_props_clfind.fits','NGC4254_co21_cube_tpeak.fits')
mkgraph.param ('NGC4254','NGC4254_co21_cube_pbcor_props_clfind.fits')
mkgraph.data('NGC4254','NGC4254_co21_cube_pbcor_props_clfind.fits',4)

#NGC4303
mkgraph.graphs ('NGC4303','NGC4303_co21_cube_pbcor_props_clfind.fits')
mkgraph.mapgmc ('NGC4303','NGC4303_co21_cube_pbcor_props_clfind.fits','NGC4303_co21_cube_tpeak.fits')
mkgraph.param ('NGC4303','NGC4303_co21_cube_pbcor_props_clfind.fits')
mkgraph.data('NGC4303','NGC4303_co21_cube_pbcor_props_clfind.fits',4)

#NGC4321
mkgraph.graphs ('NGC4321','NGC4321_co21_cube_pbcor_props_clfind.fits')
mkgraph.mapgmc ('NGC4321','NGC4321_co21_cube_pbcor_props_clfind.fits','NGC4321_co21_cube_tpeak.fits')
mkgraph.param ('NGC4321','NGC4321_co21_cube_pbcor_props_clfind.fits')
mkgraph.data('NGC4321','NGC4321_co21_cube_pbcor_props_clfind.fits',4)

#NGC4535
mkgraph.graphs ('NGC4535','NGC4535_co21_cube_pbcor_props_clfind.fits')
mkgraph.mapgmc ('NGC4535','NGC4535_co21_cube_pbcor_props_clfind.fits','NGC4535_co21_cube_tpeak.fits')
mkgraph.param ('NGC4535','NGC4535_co21_cube_pbcor_props_clfind.fits')
mkgraph.data('NGC4535','NGC4535_co21_cube_pbcor_props_clfind.fits',3)

#NGC5068
mkgraph.graphs ('NGC5068','NGC5068_co21_cube_pbcor_props_clfind.fits')
mkgraph.mapgmc ('NGC5068','NGC5068_co21_cube_pbcor_props_clfind.fits','NGC5068_co21_cube_tpeak.fits')
mkgraph.param ('NGC5068','NGC5068_co21_cube_pbcor_props_clfind.fits')
mkgraph.data('NGC5068','NGC5068_co21_cube_pbcor_props_clfind.fits',4)

#NGC6744
mkgraph.graphs ('NGC6744','NGC6744_co21_cube_pbcor_props_clfind.fits')
mkgraph.mapgmc ('NGC6744','NGC6744_co21_cube_pbcor_props_clfind.fits','NGC6744_co21_cube_tpeak.fits')
mkgraph.param ('NGC6744','NGC6744_co21_cube_pbcor_props_clfind.fits')
mkgraph.data('NGC6744','NGC6744_co21_cube_pbcor_props_clfind.fits',3)