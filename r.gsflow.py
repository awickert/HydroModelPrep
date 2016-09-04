#! /usr/bin/env python

# Based off of r.prms.py
# (A. Wickert, last update = 10 December, 2015)
# But updated to include functions for MODFLOW
# And hopefully something to solve this drainage basin delineation problem!

import os
import glob
import numpy as np
#from Scientific.IO.NetCDF import NetCDFFile
from scipy.io import loadmat
import time
import grass.script as grass
from grass.script import array as garray
from scipy.interpolate import interp1d
import multiprocessing
import re
from grass.pygrass.vector import VectorTopo

#basin_mouth_E = 209254.357775
#basin_mouth_N = 7284428.06681
basin_mouth_E = 226803.784295
basin_mouth_N = 7247299.0249

# Get r.stream.extract, r.stream.basins

# Create srtm from r.in.srtm.region (get this too)

# Then define the area for which drainage calculations will be made.
grass.run_command('g.region', rast='srtm')
#grass.run_command('g.region', n=7440000, s=7230000, w=120000, e=250000)
grass.run_command('g.region', n=7350000, s=7200000, w=170000, e=260000)

grass.run_command('r.relief', input='srtm', output='srtm.shaded', overwrite=True)

# Made to work on a projected coordinate system
reg = grass.region()

grass.mapcalc('cellArea_km2 = '+str(reg['nsres'] * reg['ewres'])+' / 10.^6', overwrite=True)

# Create drainage direction, flow accumulation, and rivers
# Creates sub-basins, river headwaters, and river segments, all with the same ID's
grass.run_command('r.watershed', elev='srtm', flow='cellArea_km2', accumulation='flowAccum', drainage='drainageDirection', stream='streams', threshold=100, flags='s', overwrite=True)
grass.run_command('r.stream.basins', direction='drainageDirection', stream_rast='streams', basins='basins', overwrite=True)
# Now have to convert streams to vector, or so I think ... maybe it doesn't matter.
# If so, omitting r.stream.extract solves basins problem!
# But now have rivers with different cat flowing next to each other, and r.to.vect combines them (unfortunately.
# This change of course breaks the "streams"-related functionalities below; will have to be fixed, including the combination of points and links in the streams network.
# I should check the code I wrote for Kelly Monteleone's paper -- this has river identification and extraction, including intersection points.

# Vectorize drainage basins
grass.run_command('r.to.vect', input='basins', output='basins', type='area', flags='sv', overwrite=True)

# Then remove all sub-basins and segments that have negative flow accumulation
# (i.e. have contributions from outside the map)

###################################################################
# Intermediate step: Remove all basins that have offmap flow
# i.e., those containing cells with negative flow accumulation
###################################################################

# Method 3 -- even easier
grass.mapcalc("has_offmap_flow = (flowAccum < 0)", overwrite=True)
grass.run_command('r.null', map='has_offmap_flow', setnull=0)
grass.run_command('r.to.vect', input='has_offmap_flow', output='has_offmap_flow', type='point', overwrite=True)
grass.run_command('r.to.vect', input='has_offmap_flow', output='has_offmap_flow', type='point', overwrite=True)
grass.run_command('v.db.addcolumn', map='has_offmap_flow', columns='badbasin_cats integer')
grass.run_command('v.what.vect', map='has_offmap_flow', column='badbasin_cats', query_map='basins', query_column='cat', dmax=60)
colNames = np.array(grass.vector_db_select('has_offmap_flow', layer=1)['columns'])
# offmap incoming flow points
colValues = np.array(grass.vector_db_select('has_offmap_flow', layer=1)['values'].values())
badcats = colValues[:,colNames == 'badbasin_cats'].squeeze()
badcats = badcats[badcats != '']
badcats = badcats.astype(int)
badcats = list(set(list(badcats)))
# basins for full cat list
colNames = np.array(grass.vector_db_select('basins', layer=1)['columns'])
colValues = np.array(grass.vector_db_select('basins', layer=1)['values'].values())
allcats = colValues[:,colNames == 'cat'].astype(int).squeeze()
allcats = list(set(list(allcats)))
# xor to goodcats
#goodcats = set(badcats).symmetric_difference(allcats)
# but better in case somehow there are badcats that are not allcats to do NOT
goodcats = list(set(allcats) - set(badcats))
goodcats_str = ''
for cat in goodcats:
  goodcats_str += str(cat) + ','
goodcats_str = goodcats_str[:-1] # super inefficient but quick
grass.run_command('g.rename', vect='basins,tmp', overwrite=True)
grass.run_command('v.extract', input='tmp', output='basins', cats=goodcats_str)
grass.run_command('g.rename', vect='streams,tmp', overwrite=True)
grass.run_command('v.extract', input='tmp', output='streams', cats=goodcats_str)


# Clean small areas -- pixellation
# 20 might be a bit much!
# ~~~~~~~~~~~~ !!!!!!!!! ~~~~~~~~~~~~
# NEED A BETTER PERMANENT SOLUTION HERE -- ONE THAT VECTORIZES THROUGH
# ANGLED 1-PIEXEL-WIDE PIECES
# ~~~~~~~~~~~~ !!!!!!!!! ~~~~~~~~~~~~
reg = grass.region()
grass.run_command('g.rename', vect='basins,basins_messy', overwrite=True)
grass.run_command('v.clean', input='basins_messy', output='basins', tool='rmarea', threshold=reg['nsres']*reg['ewres'], overwrite=True)


# Optional -- choose a subset of the region in which to do the PRMS calculation
grass.run_command( 'r.water.outlet', input='drainageDirection', output='studyBasin', coordinates=str(basin_mouth_E)+','+str(basin_mouth_N) , overwrite=True)
# Vectorize
grass.run_command( 'r.to.vect', input='studyBasin', output='studyBasin', type='area', overwrite=True)
# If there are dangling areas (single-pixel?), just drop them. Not sure if this is the best way to do it
# No check for two equal areas -- if we have this, there are more fundamental problems in defining 
# a watershed in contiguous units

#"""
# ONLY IF MORE THAN ONE STUDY BASIN -- REMOVE DANLGE!
# BUT MAYBE TAKEN CARE OF BY A SIMPLER V.CLEAN?
grass.run_command( 'v.db.addcolumn', map='studyBasin', columns='area_m2 double precision' )
grass.run_command( 'v.db.dropcolumn', map='studyBasin', columns='label' )
grass.run_command( 'v.to.db', map='studyBasin', columns='area_m2', option='area', units='meters')
drainageAreasRaw = sorted( grass.parse_command( 'v.db.select', map='studyBasin', flags='c').keys() ) # could update to grass.vector_db_select
drainageAreasList = []
for row in drainageAreasRaw:
  # cat, area
  drainageAreasList.append(row.split('|'))
drainageAreasOnly = np.array(drainageAreasList).astype(float)
catsOnly = drainageAreasOnly[:,0].astype(int)
drainageAreasOnly = drainageAreasOnly[:,1]
row_with_max_drainage_area = (drainageAreasOnly == np.max(drainageAreasOnly)).nonzero()[0][0]
cat_with_max_drainage_area = catsOnly[row_with_max_drainage_area]
grass.run_command('g.rename', vect='studyBasin,tmp')
grass.run_command('v.extract', input='tmp', output='studyBasin', cats=cat_with_max_drainage_area, overwrite=True)
grass.run_command('g.remove', type='vector', name='tmp', flags='f')
grass.run_command('v.to.rast', input='studyBasin', output='studyBasin', use='val', value=1, overwrite=True)
#"""


###############
# PLACEHOLDER #
###################################################################
# To do in near future: limit to this basin
###################################################################

# Next, get the order of basins the old-fashioned way: coordinates of endpoints of lines
# Because I can't use GRASS to query multiple points
grass.run_command('v.extract', input='streams', output='streamSegments', type='line', overwrite=True)
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



#############
# Now, let's copy/rename the sub-basins to HRU and the streamSegments to segment and give them attributes
###########################################################################################################

# Attributes (in order given in manual)

# HRU
hru_columns = []
# Self ID
hru_columns.append('id integer') # nhru
# Basic Physical Attributes (Geometry)
hru_columns.append('hru_area double precision') # acres (!!!!)
hru_columns.append('hru_aspect double precision') # Mean aspect [degrees]
hru_columns.append('hru_elev double precision') # Mean elevation
hru_columns.append('hru_lat double precision') # Latitude of centroid
hru_columns.append('hru_slope double precision') # Mean slope [percent]
# Basic Physical Attributes (Other)
#hru_columns.append('hru_type integer') # 0=inactive; 1=land; 2=lake; 3=swale; almost all will be 1
#hru_columns.append('elev_units integer') # 0=feet; 1=meters. 0=default. I think I will set this to 1 by default.
# Measured input
hru_columns.append('outlet_sta integer') # Index of streamflow station at basin outlet:
                                     #   station number if it has one, 0 if not
#    Note that the below specify projections and note lat/lon; they really seem
#    to work for any projected coordinates, with _x, _y, in meters, and _xlong, 
#    _ylat, in feet (i.e. they are just northing and easting). The meters and feet
#    are not just simple conversions, but actually are required for different
#    modules in the code, and are hence redundant but intentional.
hru_columns.append('hru_x double precision') # Easting [m]
hru_columns.append('hru_xlong double precision') # Easting [feet]
hru_columns.append('hru_y double precision') # Northing [m]
hru_columns.append('hru_ylat double precision') # Northing [feet]
# Streamflow and lake routing
hru_columns.append('K_coef double precision') # Travel time of flood wave to next downstream segment;
                                              #   this is the Muskingum storage coefficient
                                              #   1.0 for reservoirs, diversions, and segments flowing
                                              #   out of the basin
hru_columns.append('x_coef double precision') # Amount of attenuation of flow wave;
                                              #   this is the Muskingum routing weighting factor
                                              #   range: 0.0--0.5; default 0.2
                                              #   0 for all segments flowing out of the basin
hru_columns.append('hru_segment integer') # ID of stream segment to which flow will be routed
                                          #   this is for non-cascade routing (flow goes directly
                                          #   from HRU to stream segment)
hru_columns.append('obsin_segment integer') # Index of measured streamflow station that replaces
                                            #   inflow to a segment

# Segments
segment_columns = []
# Self ID
segment_columns.append('id integer') # nsegment
# Streamflow and lake routing
segment_columns.append('tosegment integer') # Index of downstream segment to which a segment
                                            #   flows (thus differentiating it from hru_segment,
                                            #   which is for HRU's, though segment and HRU ID's
                                            #   are the same when HRU's are sub-basins

# PRODUCE THE DATA TABLES
##########################

# Create strings
hru_columns = ",".join(hru_columns)
segment_columns = ",".join(segment_columns)

#"""
# Copy
grass.run_command('g.copy', vect='basins,HRU', overwrite=True)
grass.run_command('g.copy', vect='streamSegments,segment', overwrite=True)
#"""

# Rename / subset
"""
# OR GO BACK TO HRU_messy
grass.run_command('v.overlay', ainput='basins', binput='studyBasin', operator='and', output='HRU_messy', overwrite=True)
grass.run_command('v.overlay', ainput='streamSegments', binput='studyBasin', operator='and', output='segment_messy', overwrite=True)
# And clean as well
grass.run_command('v.clean', input='HRU_messy', output='HRU', tool='rmarea', threshold=reg['nsres']*reg['ewres']*40, overwrite=True)
grass.run_command('v.clean', input='segment_messy', output='segment', tool='rmdangle', threshold=reg['nsres']*2, overwrite=True)
# And now that the streams and HRU's no longer have the same cat values, fix 
# this.
grass.run_command('v.db.droptable', map='HRU', flags='f')
grass.run_command('v.db.droptable', map='segment', flags='f')
#grass.run_command('v.category', input='HRU', option='del', cat='-1', out='tmp', overwrite=True)
#grass.run_command('v.category', input='tmp', option='add', out='HRU' overwrite=True)
grass.run_command('v.db.addtable', map='HRU')
grass.run_command('v.db.addtable', map='segment')

grass.run_comm


v.clean HRU
v.clean
v
v.what.vect 
"""

#grass.run_command('v.clean', input='segment_messy', output='HRU', tool='rmarea', threshold=reg['nsres']*reg['ewres']*20, overwrite=True)


# Add columns to tables
grass.run_command('v.db.addcolumn', map='HRU', columns=hru_columns)
grass.run_command('v.db.addcolumn', map='segment', columns=segment_columns)


# Produce the data table entries
##################################

"""
# ID numbers
# There should be a way to do this all at once, but...
for i in range(len(cats)):
  grass.run_command('v.db.update', map='HRU', column='id', value=nhru[i], where='cat='+str(cats[i]))
nsegment = nhru.copy() # ONLY FOR THIS SPECIAL CASE -- will be different in general
for i in range(len(cats)):
  grass.run_command('v.db.update', map='segment', column='id', value=nsegment[i], where='cat='+str(cats[i]))
"""

nhru = np.arange(1, xy1.shape[0]+1)
nhrut = []
for i in range(len(nhru)):
  nhrut.append( (nhru[i], cats[i]) )
# Access the HRU's 
hru = VectorTopo('HRU')
# Open the map with topology:
hru.open('rw')
# Create a cursor
cur = hru.table.conn.cursor()
# Use it to loop across the table
cur.executemany("update HRU set id=? where cat=?", nhrut)
# Commit changes to the table
hru.table.conn.commit()
# Close the table
hru.close()

# if you want to append to table
# cur.executemany("update HRU(id) values(?)", nhrut) # "insert into" will add rows

# Same for segments
nsegment = nhru.copy() # ONLY FOR THIS SPECIAL CASE -- will be different in general
nsegmentt = nhrut # ONLY FOR THIS SPECIAL CASE -- will be different in general

# Somehow only works after I v.clean, not right after v.overlay
segment = VectorTopo('segment')
segment.open('rw')
cur = segment.table.conn.cursor()
cur.executemany("update segment set id=? where cat=?", nsegmentt)
segment.table.conn.commit()
segment.close()

#hru_columns.append('hru_area double precision')
grass.run_command('v.to.db', map='HRU', option='area', columns='hru_area', units='acres')

# GET MEAN VALUES FOR THESE NEXT ONES, ACROSS THE BASIN

# hru_columns.append('hru_aspect double precision') # Mean aspect [degrees]
# hru_columns.append('hru_slope double precision') # Mean slope [percent]
# Slope
grass.run_command('r.slope.aspect', elevation='srtm', slope='tmp', aspect='aspect', format='percent', overwrite=True) # zscale=0.01 also works to make percent be decimal 0-1
grass.mapcalc('slope = tmp / 100.', overwrite=True)
grass.run_command('v.rast.stats', map='HRU', raster='slope', method='average', column_prefix='tmp', flags='c')
grass.run_command('v.db.update', map='HRU', column='hru_slope', query_column='tmp_average')
grass.run_command('v.db.dropcolumn', map='HRU', column='tmp_average')
# Dealing with conversion from degrees (no good average) to something I can
# average -- x- and y-vectors
# Geographic coordinates, so sin=x, cos=y.... not that it matters so long 
# as I am consistent in how I return to degrees
grass.mapcalc('aspect_x = sin(aspect)', overwrite=True)
grass.mapcalc('aspect_y = cos(aspect)', overwrite=True)
#grass.run_command('v.db.addcolumn', map='HRU', columns='aspect_x_sum double precision, aspect_y_sum double precision, ncells_in_hru integer')
grass.run_command('v.rast.stats', map='HRU', raster='aspect_x', method='sum', column_prefix='aspect_x', flags='c')
grass.run_command('v.rast.stats', map='HRU', raster='aspect_y', method='sum', column_prefix='aspect_y', flags='c')
# Not actually needed, but maybe good to know
#grass.run_command('v.rast.stats', map='HRU', raster='aspect_y', method='number', column_prefix='tmp', flags='c')
#grass.run_command('v.db.renamecolumn', map='HRU', column='tmp_number,ncells_in_hru')
# NO TRIG FUNCTIONS IN SQLITE!
#grass.run_command('v.db.update', map='HRU', column='hru_aspect', query_column='DEGREES(ATN2(aspect_y_sum, aspect_x_sum))') # Getting 0, why?
hru = VectorTopo('HRU')
hru.open('rw')
cur = hru.table.conn.cursor()
cur.execute("SELECT cat,aspect_x_sum,aspect_y_sum FROM %s" %hru.name)
_arr = np.array(cur.fetchall())
_cat = _arr[:,0]
_aspect_x_sum = _arr[:,1]
_aspect_y_sum = _arr[:,2]
aspect_angle = np.arctan2(_aspect_y_sum, _aspect_x_sum) * 180./np.pi
aspect_angle[aspect_angle < 0] += 360 # all positive
aspect_angle_cat = np.vstack((aspect_angle, _cat)).transpose()
cur.executemany("update HRU set hru_aspect=? where cat=?", aspect_angle_cat)
hru.table.conn.commit()
hru.close()

# hru_columns.append('hru_elev double precision') # Mean elevation
grass.run_command('v.rast.stats', map='HRU', raster='srtm', method='average', column='tmp', flags='c')
grass.run_command('v.db.update', map='HRU', column='hru_elev', query_column='tmp_average')
grass.run_command('v.db.dropcolumn', map='HRU', column='tmp_average')

# get x,y of centroid -- but have areas not in database table, that do have
# centroids, and having a hard time finding a good way to get rid of them!
# They have duplicate category values!
# Perhaps these are little dangles on the edges of the vectorization where
# the raster value was the same but pinched out into 1-a few cells?
# From looking at map, lots of extra centroids on area boundaries, and removing
# small areas (though threshold hard to guess) gets rid of these

"""
g.copy vect=HRU,HRUorig # HACK!!!
v.clean in=HRUorig out=HRU tool=rmarea --o thresh=15000
"""

#grass.run_command( 'g.rename', vect='HRU,HRU_too_many_centroids')
#grass.run_command( 'v.clean', input='HRU_too_many_centroids', output='HRU', tool='rmdac')
grass.run_command('v.db.addcolumn', map='HRU', columns='centroid_x double precision, centroid_y double precision')
grass.run_command( 'v.to.db', map='HRU', type='centroid', columns='centroid_x, centroid_y', option='coor', units='meters')

# hru_columns.append('hru_lat double precision') # Latitude of centroid
colNames = np.array(grass.vector_db_select('HRU', layer=1)['columns'])
colValues = np.array(grass.vector_db_select('HRU', layer=1)['values'].values())
xy = colValues[:,(colNames=='centroid_x') + (colNames=='centroid_y')]
np.savetxt('_xy.txt', xy, delimiter='|', fmt='%s')
grass.run_command('m.proj', flags='od', input='_xy.txt', output='_lonlat.txt', overwrite=True)
lonlat = np.genfromtxt('_lonlat.txt', delimiter='|',)[:,:2]
lonlat_cat = np.concatenate((lonlat, np.expand_dims(_cat, 2)), axis=1)

# why not just get lon too?
grass.run_command('v.db.addcolumn', map='HRU', columns='hru_lon double precision')

hru = VectorTopo('HRU')
hru.open('rw')
cur = hru.table.conn.cursor()
cur.executemany("update HRU set hru_lon=?, hru_lat=? where cat=?", lonlat_cat)
hru.table.conn.commit()
hru.close()

# Easting and Northing for other columns
grass.run_command('v.db.update', map='HRU', column='hru_x', query_column='centroid_x')
grass.run_command('v.db.update', map='HRU', column='hru_xlong', query_column='centroid_x*3.28084') # feet
grass.run_command('v.db.update', map='HRU', column='hru_y', query_column='centroid_y')
grass.run_command('v.db.update', map='HRU', column='hru_ylat', query_column='centroid_y*3.28084') # feet


# Streamflow and lake routing
# tosegment
"""
# THIS IS THE NECESSARY PART
# CHANGED (BELOW) TO RE-DEFINE NUMBERS IN SEQUENCE AS HRU'S INSTEAD OF USING
# THE CAT VALUES
# Get the first channels in the segment
tosegment = np.zeros(len(cats)) # default to 0 if they do not flow to another segment
# Loop over all segments
#for i in range(len(cats)):
# From outlet segment
for i in range(len(xy2)):
  # to inlet segment
  inlets = np.prod(xy1 == xy2[i], axis=1)
  # Update inlet segments with ID of outlets
  tosegment[inlets.nonzero()] = cats[i]
tosegment_cat = tosegment.copy()
"""

tosegment_cats = np.zeros(len(cats)).astype(int) # default to 0 if they do not flow to another segment
tosegment = np.zeros(len(cats)).astype(int) # default to 0 if they do not flow to another segment
# From outlet segment
for i in range(len(xy2)):
  # to outlet segment
  outlets = np.prod(xy2 == xy1[i], axis=1)
  # Update outlet segments with ID of inlets
  tosegment[outlets.nonzero()] = nhru[i]
  tosegment_cats[outlets.nonzero()] = cats[i]

"""
  # BACKWARDS!
  # to inlet segment
  inlets = np.prod(xy1 == xy2[i], axis=1)
  # Update inlet segments with ID of outlets
  tosegment_cats[inlets.nonzero()] = cats[i]
"""

# Now, just update tosegment (segments) and hru_segment (hru's)
# In this case, they are the same.
nsegment = nhru.copy() # ONLY FOR THIS SPECIAL CASE -- will be different in general
nsegmentt = nhrut # ONLY FOR THIS SPECIAL CASE -- will be different in general
# Tuple for upload to SQL
# 0 is the default value if it doesn't go into any other segment (i.e flows
# off-map)
tosegmentt = []
tosegment_cats_t = []
for i in range(len(nsegment)):
  tosegmentt.append( (tosegment[i], nsegment[i]) )
  tosegment_cats_t.append( (tosegment_cats[i], cats[i]) )
# Once again, special case
hru_segmentt = tosegmentt

# Loop check!
# Weak loop checker - will only detect direct ping-pong.
loops = []
tosegmenta = np.array(tosegmentt)
for i in range(len(tosegmenta)):
  for j in range(len(tosegmenta)):
    if (tosegmenta[i] == tosegmenta[j][::-1]).all():
      loops.append(tosegmenta[i])

segment = VectorTopo('segment')
segment.open('rw')
cur = segment.table.conn.cursor()
cur.executemany("update segment set tosegment=? where id=?", tosegmentt)
segment.table.conn.commit()
segment.close()

hru = VectorTopo('HRU')
hru.open('rw')
cur = hru.table.conn.cursor()
cur.executemany("update HRU set hru_segment=? where id=?", hru_segmentt)
hru.table.conn.commit()
hru.close()


#grass.run_command('g.rename', vect='HRU_all_2,HRU', overwrite=True)
#grass.run_command('g.rename', vect='segment_all_2,segment', overwrite=True)

# In study basin?
grass.run_command('v.db.addcolumn', map='segment', columns='in_study_basin int')
grass.run_command('v.db.addcolumn', map='HRU', columns='in_study_basin int')
grass.run_command('v.what.vect', map='segment', column='in_study_basin', query_map='studyBasin', query_column='value')
grass.run_command('v.what.vect', map='HRU', column='in_study_basin', query_map='segment', query_column='in_study_basin')

# Save global segment+HRU
grass.run_command('g.rename', vect='HRU,HRU_all')
grass.run_command('g.rename', vect='segment,segment_all')

# Output HRU -- will need to ensure that this is robust!
grass.run_command('v.extract', input='HRU_all', output='HRU', where='in_study_basin=1', overwrite=True)
grass.run_command('v.extract', input='segment_all', output='segment', where='in_study_basin=1', overwrite=True)


colNames = np.array(grass.vector_db_select('segment')['columns'])
colValues = np.array(grass.vector_db_select('segment')['values'].values())
cats = colValues[:,colNames == 'cat'].astype(int).squeeze()
xy1 = colValues[:,(colNames == 'x1') + (colNames == 'y1')].astype(float)
xy2 = colValues[:,(colNames == 'x2') + (colNames == 'y2')].astype(float)
xy  = np.vstack((xy1, xy2))

# Redo nhru down here
nhru = np.arange(1, xy1.shape[0]+1)
nhrut = []
for i in range(len(nhru)):
  nhrut.append( (nhru[i], cats[i]) )
  """
  n = 1
  if i != 1:
    nhrut.append( (n, cats[i]) )
    n += 1
  """
  
hru = VectorTopo('HRU')
hru.open('rw')
cur = hru.table.conn.cursor()
cur.executemany("update HRU set id=? where cat=?", nhrut)
hru.table.conn.commit()
hru.close()

# if you want to append to table
# cur.executemany("update HRU(id) values(?)", nhrut) # "insert into" will add rows

# Same for segments
nsegment = nhru.copy() # ONLY FOR THIS SPECIAL CASE -- will be different in general
nsegmentt = nhrut # ONLY FOR THIS SPECIAL CASE -- will be different in general

# Somehow only works after I v.clean, not right after v.overlay
segment = VectorTopo('segment')
segment.open('rw')
cur = segment.table.conn.cursor()
cur.executemany("update segment set id=? where cat=?", nsegmentt)
segment.table.conn.commit()
segment.close()


tosegment_cats = np.zeros(len(cats)).astype(int) # default to 0 if they do not flow to another segment
tosegment = np.zeros(len(cats)).astype(int) # default to 0 if they do not flow to another segment
# From outlet segment
for i in range(len(xy2)):
  # to outlet segment
  outlets = np.prod(xy2 == xy1[i], axis=1)
  # Update outlet segments with ID of inlets
  tosegment[outlets.nonzero()] = nhru[i]
  tosegment_cats[outlets.nonzero()] = cats[i]

# Now, just update tosegment (segments) and hru_segment (hru's)
# In this case, they are the same.
nsegment = nhru.copy() # ONLY FOR THIS SPECIAL CASE -- will be different in general
nsegmentt = nhrut # ONLY FOR THIS SPECIAL CASE -- will be different in general
# Tuple for upload to SQL
# 0 is the default value if it doesn't go into any other segment (i.e flows
# off-map)
tosegmentt = []
tosegment_cats_t = []
for i in range(len(nsegment)):
  tosegmentt.append( (tosegment[i], nsegment[i]) )
  tosegment_cats_t.append( (tosegment_cats[i], cats[i]) )
# Once again, special case
hru_segmentt = tosegmentt

# Loop check!
# Weak loop checker - will only detect direct ping-pong.
loops = []
tosegmenta = np.array(tosegmentt)
for i in range(len(tosegmenta)):
  for j in range(len(tosegmenta)):
    if (tosegmenta[i] == tosegmenta[j][::-1]).all():
      loops.append(tosegmenta[i])


segment = VectorTopo('segment')
segment.open('rw')
cur = segment.table.conn.cursor()
cur.executemany("update segment set tosegment=? where id=?", tosegmentt)
segment.table.conn.commit()
segment.close()

hru = VectorTopo('HRU')
hru.open('rw')
cur = hru.table.conn.cursor()
cur.executemany("update HRU set hru_segment=? where id=?", hru_segmentt)
hru.table.conn.commit()
hru.close()

v.db.select segment sep=comma > segment.csv
v.db.select HRU sep=comma > HRU.csv
# and then sort by id, manually
# And then manually change the last segment's "tosegment" to 0.
# Except in this case, it was 0!
# Maybe I managed to do this automatically above... but tired and late, 
# so will check later
# but hoping I did something right by re-doing all of the above before
# saving (and doing so inside this smaller basin)

print "COMPLETE."
















































"""
# Extract only that region that is within the study area
grass.run_command('g.rename', vect='HRU,HRU_all')
grass.run_command('g.rename', vect='segment,segment_all')
#grass.run_command('v.extract', input=

# Extract only that region that is within the study area
grass.run_command('g.rename', vect='HRU,HRU_all_2')
grass.run_command('g.rename', vect='segment,segment_all_2')
# 
"""














# query number of topological features
areas   = zipcodes.number_of("areas")
islands = zipcodes.number_of("islands")
print 'Map: <' + vectmap + '> with %d areas and %d islands' % (areas, islands)
 

# http://www.ing.unitn.it/~zambelli/projects/pygrass/attributes.html

  
# This is maybe not true -- but intention is to send the HRU's flow to its
# enclosed river segment for river segment ID === HRU ID. But need consecutive
# numbers [1, 2, ..., n]
# hru_segment = tosegment.copy()

                                            
# Parameters for cascading-flow simulation
# Not doing cascade yet, but preparing for it!
#hru_columns.append('cascade_flg integer') # Cascade type: 0: many to many; 1: one to one
#hru_columns.append('cascade_tol double precision') # Cascade area below which cascade links are ignored
#hru_columns.append('circle_switch integer') # Switch to check for circles (0, don't = efficient; 1, do = safe, good at first!)

#hru_columns.append('gw_down_id integer') # START HERE!!!!!!!!!!
#hru_columns.append('cascade_flg double precision') # 
#hru_columns.append('cascade_flg double precision') # 
#hru_columns.append('cascade_flg double precision') # 
#hru_columns.append('cascade_flg double precision') # 
#hru_columns.append('cascade_flg integer') # 

# hru_up_id --> ncascade
# hru_down_id --> ncascade
# hru_strmseg_down_id --> ncascade
#hru_columns.append('ncascade integer') # CHECK WITH CRYSTAL WHAT THIS SHOULD BE -- THINK IT IS THE SURFACE-WATER ONE
#hru_columns.append('ncascade_gw integer') # CHECK WITH CRYSTAL WHAT THIS SHOULD BE -- THINK IT IS GW_DOWN_ID OR SOMETHING LIKE THIS



"""
# NOT SURE IF THIS IS NECESSARY ANYMORE!

# First, get accumulation at points -- see which is upstream and which is downstream.
# And get elevation while we're at it
xystr = []
for x, y in xy:
  #print x, y
  xystr.append(str(x)+','+str(y))

# Really inefficient, but first step is just to make it work...
z = []
accum = []  
for xyi in xystr:
  z.append(float(grass.parse_command('r.what', map='srtm,flowAccum', coordinates=xyi).keys()[0].split('|')[-2]))
  accum.append(float(grass.parse_command('r.what', map='srtm,flowAccum', coordinates=xyi).keys()[0].split('|')[-1]))
z = np.array(z)
accum_raw = np.array(accum)
accum = np.abs(accum) # This will include some negative values of basin mouths that enter rivers with offmap inputs
                      # Making positive because everything 

# break big lists into those corresponding to x1,y1; x2, y2
L = len(xy1)
accum1 = accum[:L]
accum2 = accum[L:]
z1 = z[:L]
z2 = z[L:]

# And likewise turn small lists into big ones
cat2 = np.vstack((cat, cat))

streamSegments = 
"""


# No longer needed
#incomplete_basin_cats = list(set(list(cat2[accum <= 0].squeeze())))

"""
from contextlib import redirect_stdout
# check: http://eli.thegreenplace.net/2015/redirecting-all-kinds-of-stdout-in-python/
xystr = ''
for x, y in xy:
  #print x, y
  xystr += str(x)+' '+str(y)+'\n'
grass.write_command('r.what', map='srtm,flowAccum', output='-', stdin=xystr)
"""


nsegs_at_point = []
for row in xy:
  nsegs_at_point.append(np.sum( np.prod(xy == row, axis=1)))
nsegs_at_point = np.array(nsegs_at_point)

# Then, see which upstream ones are the upstream-most (i.e., no flow entering)

# Then, go to the ends of these and see what they flow into.

# Fill an array of the numbers of the in_vectors and out_vectors



# Each row in xy1 corresponds to each row in xy2, so these are the endpoints of this stream segment


# CHECK FOR REPEATS -- SHOULD CREATE A LIST OF MINIMUM POINTS
nsegs_at_point = []
for row in xy1:
  nsegs_at_point.append(np.sum( np.prod(xy1 == row, axis=1) + np.prod(xy2 == row, axis=1) ))
for row in xy2:
  nsegs_at_point.append(np.sum((xy1 == row) + (xy2 == row))/2  - 1)
nsegs_at_point = np.array(nsegs_at_point)
(nsegs_at_point == 5).nonzero()

all_points = np.vstack((xy1, xy2))
all_unique_points

# Have to query using coordinates, later
#grass.run_command('v.what.rast', map=

"""
# Get the order of basins
grass.run_command('v.extract', input='streams', output='streamHeads', type='point')
grass.run_command('v.db.addcolumn', map='streamHeads', columns='streamSegmentsHere integer')
grass.run_command('v.what.vect', map='streamHeads', column='streamSegmentsHere', query_map='streams', query_column='cat', dmax=60)

grass.run_command('v.extract', input='streams', output='streamSegments', type='line')
"""

# CASCADE: SKIP THIS IF YOU WANT TO GO DIRECTLY TO BASIN OUTLET (AND IF YOU USE SUBBASINS)









###########
# MODFLOW #
###########

# Generate coarse box for MODFLOW (ADW, 4 September, 2016)

grass.run_command('g.region', rast='srtm')
grass.run_command('g.region', n=7350000, s=7200000, w=170000, e=260000)
reg = grass.region()
MODFLOWres = 2000.
grass.run_command('v.to.rast', input='HRU', output='allHRUs', use='val', val=1.0, overwrite=True)
grass.run_command('r.null', map='allHRUs', null='0')
grass.run_command('r.colors', map='allHRUs', color='grey', flags='n')
grass.run_command('g.region', res=MODFLOWres)
grass.run_command('r.resamp.stats', method='average', input='allHRUs', output='fraction_of_HRU_in_MODFLOW_cell', overwrite=True)
grass.run_command('r.colors', map='fraction_of_HRU_in_MODFLOW_cell', color='grey', flags='n')


