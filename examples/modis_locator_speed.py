from csat2 import MODIS
import matplotlib.pyplot as plt
import copy
import numpy as np
import time

gran = MODIS.Granule.fromtext("2015001.1220A")
gran.download('021KM')

pattern_locs = np.array(
    np.meshgrid(np.arange(5, 2100, 100),
                np.arange(5, 1350, 65))).T
geo_locs = gran.geolocate(pattern_locs).reshape((-1, 2))

numbers = [1, 3, 10, 30, 100]
locators = ['BallTree',
            'BallTree1km',
            'FullSearch',
            'SphereRemap',
            'SphereRemap2km'
            ]

timing = {}
for locator_type in locators:
    print(locator_type)
    timing[locator_type] = []
    for number in numbers:
        gran_test = copy.deepcopy(gran)

        stime = time.time()
        calc_plocs = gran_test.locate(
            geo_locs[:number],
            locator_type=locator_type
        )
        etime = time.time() - stime
        timing[locator_type].append([number, etime])
        print(number, etime)

print('{name: <15} '.format(name='NumPoints')+' '.join([f'{n: >6}' for n in numbers]))

print('')
for name in locators:
    print(f'{name: <15} '+' '.join(['{n:6.2f}'.format(n=timing[name][n][1]) for n in range(len(numbers))]))
