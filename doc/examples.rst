Examples
========


Readin and plot a MODIS granule
-------------------------------

.. code-block:: python

   import csat2.MODIS
   import matplotlib.pyplot as plt
   import cartopy.crs as ccrs
   import numpy as np

   # Defaults to collection 6.1
   gran = csat2.MODIS.Granule.fromtext('2015011.2110A')

   # Downloads the file into the location specified by the machine file
   # This requires you to have setup an Earthdata account and put your
   # key in your csat2 config folder
   gran.download('06_L2')

   cer = gran.get_variable('06_L2', ['Cloud_Effective_Radius'])

   # Get the lat and lon for the data array
   lon, lat = gran.get_lonlat()

   # Mask missing data
   toplot = cer['Cloud_Effective_Radius']
   toplot = np.ma.array(toplot, mask=np.isnan(toplot))

   plt.subplot(111, projection=ccrs.PlateCarree())
   plt.pcolormesh(lon, lat, toplot)
   plt.gca().coastlines()
   plt.colorbar(shrink=0.7, label=r'Effective Radius ($\mu$m)')
   plt.title(str(gran))
   
   fig = plt.gcf()
   fig.savefig('cer_{}.png'.format(gran.astext()), bbox_inches='tight')
   plt.show()


Readin and plot and FCI image
-----------------------------

.. code-block:: python

   from satpy.writers import to_image
   from csat2.FCI import Granule

   gran = Granule.fromtext("MTG1.2025024.1050.FD.FDHSI")
   gran.download()
   scn = gran.get_satpy_scene()
   scn.load(['ash'], upper_right_corner='NE')
   new_scn = scn.resample("eurol1")
   img = to_image(new_scn['ash'])
   img.stretch("crude", min_stretch=[-4,-4,243], max_stretch=[2, 5, 303])

   # %%
   import cartopy.crs as ccrs
   import matplotlib.pyplot as plt
   fig, ax = plt.subplots(subplot_kw={'projection': new_scn['ash'].attrs['area'].to_cartopy_crs()})
   img.data.plot.imshow(ax=ax)
   ax.coastlines()
   ax.gridlines(draw_labels=True)
   ax.set_extent([-3, 0, 50, 53], crs=ccrs.PlateCarree())
   ax.set_title("2025-01-24 10:50 UTC")
   # img.save("ash_contrail_radar.png")
