#! /usr/bin/env python

#!/usr/bin/env python
############################################################################
#
# MODULE:       r.streams.topological
#
# AUTHOR(S):    Andrew Wickert
#
# PURPOSE:      Compute a topologically consistent vector stream network
#               from a flow accumulation raster
#
# COPYRIGHT:    (c) 2017 Andrew Wickert
#
#               This program is free software under the GNU General Public
#               License (>=v2). Read the file COPYING that comes with GRASS
#               for details.
#
#############################################################################
#
# REQUIREMENTS:
#      -  r.stream.basins
# 
# NOTES:
# Based off of r.prms.py
# (A. Wickert, last update = 10 December, 2015)
# But updated to include functions for MODFLOW
# And hopefully something to solve this drainage basin delineation problem!
 
#%module
#% description: Input topology and attributes for GSFLOW
#% keyword: raster
#% keyword: hydrology
#%end

#%flag
#%  key: l
#%  description: Allows running in lat/lon: dx is f(lat) at grid N-S midpoint
#%end

#%option G_OPT_R_INPUT
#%  key: input
#%  type: string
#%  description: Flow accumulation
#%  required : yes
#%end

#%option
#%  key: threshold
#%  type: double
#%  description: Threshold drainage area to define streams
#%  required : yes
#%end

##################
# IMPORT MODULES #
##################

# PYTHON
import os
import glob
import numpy as np
# GRASS
import grass.script as grass
from grass.script import array as garray
from grass.pygrass.vector import VectorTopo

def main():
    """
    Input for GSFLOW
    """

    reg = grass.region()

    options, flags = grass.parser()

    basin_mouth_E = options['E']
    basin_mouth_N = options['N']

    accum_thresh = options['threshold']

    # Create drainage direction, flow accumulation, and rivers

    grass.mapcalc('streams_unthinned = flowAccum > '+str(accum_thresh), overwrite=True)
    grass.run_command('r.null', map='streams_unthinned', setnull=0)
    grass.run_command('r.thin', input='streams_unthinned', output='streams', overwrite=True)
    grass.run_command('r.to.vect', input='streams', output='streams_raw', type='line', overwrite=True)
    # Clean with a 1-cell threshold to remove loops created when diagonal
    # streams intersect one another
    grass.run_command('v.clean', input='streams_raw', output='streams', tool='snap', threshold=1.42*(grass.region()['nsres'] + grass.region()['ewres'])/2., flags='c', overwrite=True)
    grass.run_command('v.to.rast', input='streams', output='streams_unthinned', use='val', val=1, overwrite=True)
    grass.run_command('r.thin', input='streams_unthinned', output='streams', overwrite=True)
    grass.run_command('r.to.vect', input='streams', output='streams', type='line', overwrite=True)
    grass.run_command('v.to.rast', input='streams', output='streams', use='cat', overwrite=True)
    

    # Include this?

    ###############
    # PLACEHOLDER #
    ###################################################################
    # To do in near future: limit to this basin
    ###################################################################

    # Next, get the order of basins the old-fashioned way: coordinates of endpoints of lines
    # Because I can't use GRASS to query multiple points
    #grass.run_command('v.extract', input='streams', output='streamSegments', type='line', overwrite=True)
    # Maybe I don't even need nodes! 9/4/16 -- nope, doesn't seem so.
    grass.run_command('g.copy', rast='streams,streamSegments')
    grass.run_command('v.db.addcolumn', map='streamSegments', columns='z double precision, flow_accum double precision, x1 double precision, y1 double precision, x2 double precision, y2 double precision')
    grass.run_command('v.to.db', map='streamSegments', option='start', columns='x1, y1')
    grass.run_command('v.to.db', map='streamSegments', option='end', columns='x2, y2')

    colNames = np.array(grass.vector_db_select('streamSegments')['columns'])
    colValues = np.array(grass.vector_db_select('streamSegments')['values'].values())
    cats = colValues[:,colNames == 'cat'].astype(int).squeeze()
    xy1 = colValues[:,(colNames == 'x1') + (colNames == 'y1')].astype(float)
    xy2 = colValues[:,(colNames == 'x2') + (colNames == 'y2')].astype(float)
    xy  = np.vstack((xy1, xy2))

    # xy1: UPSTREAM
    # xy2: DOWNSTREAM
    # (I checked.)
    # So now can use this information to find headwaters and mouths

    # Not sure that thsi is necessary
    nsegs_at_point_1 = []
    nsegs_at_point_2 = []
    for row in xy1:
      nsegs_at_point_1.append(np.sum( np.prod(xy == row, axis=1)))
    for row in xy2:
      nsegs_at_point_2.append(np.sum( np.prod(xy == row, axis=1)))
    nsegs_at_point_1 = np.array(nsegs_at_point_1)
    nsegs_at_point_2 = np.array(nsegs_at_point_2)


    # HRU's have same numbers as their enclosed segments
    # NOT TRUE IN GENERAL -- JUST FOR THIS CASE WITH SUB-BASINS -- WILL NEED TO FIX IN FUTURE



