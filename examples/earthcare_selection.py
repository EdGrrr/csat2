from csat2 import EarthCARE

gran = EarthCARE.Granule('06063', 'E', stream=True)
gran.datetime()
files = EarthCARE.download.granule_files(gran)

top_level = None
top_level = 'A'

proc_files = [[a['properties']['product:type'],
               a['properties']['version'],
               a['assets']['enclosure_h5']['href'],
               a['assets']['product']['href']] for a in files]

if top_level:
    proc_files = [a for a in proc_files if a[1].startswith(top_level)]

recent = {}
for pf in proc_files:
    if pf[1] > recent.get(pf[0], ''):
        recent[pf[0]] = pf[1]

pd.DataFrame(proc_files).groupby(0).max()
