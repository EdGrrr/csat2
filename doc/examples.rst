Examples
========


Readin and plot a MODIS granule
-------------------------------

.. code-block:: python

   import csat2.MODIS
   import matplotlib.pyplot as plt
   import cartopy.crs as ccrs
   import numpy as np
   
   gran = csat2.MODIS.Granule.fromtext('2015011.2110A')

   cer = gran.get_variable('06_L2', ['Cloud_Effective_Radius'])

   lon, lat = gran.get_lonlat()

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
