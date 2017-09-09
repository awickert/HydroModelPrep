r.cell.area output=cellArea_km2 units=km2

grass.run_command('r.watershed', elev='srtm', flow='cellArea_km2', accumulation='flowAccum', drainage='drainageDirection', flags='s', overwrite=True)

r.streams.topological

grass.run_command('r.stream.basins', direction='drainageDirection', stream_rast='streams', basins='basins', overwrite=True)
grass.run_command('r.to.vect', input='basins', output='basins', type='area', flags='v', overwrite=True)

v.no.offmap.flow # define as needing the above inputs?
                 # create my own stream module set?
