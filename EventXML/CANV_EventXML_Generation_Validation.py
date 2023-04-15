# gcmt_search.py
# Search GCMT CSV file for events

import pandas as pd
from obspy.core import UTCDateTime
from obspy.clients.fdsn import Client
from obspy.core.event import Catalog
from obspy.core.event.event import Event
import sys

from obspy.core.event import Origin, Magnitude
from obspy.core.event.source import FocalMechanism, MomentTensor, Tensor
from math import sqrt, log10

def mw_from_m0(m0):
    """
    Compute mw (moment magnitude) from m0 (seismic moment)
    Assumes m0 is in MKS units, Newton-meters
    """
    return (log10(m0) - 9.1)/1.5

def get_gcmt(gcmt):
    """ 
    Outputs gcmt origin and focal mechanism objects to be added the event metadata
    Based on obspy.io.cmtsolution.core
    Centroid origin
    """
    origin = Origin(time = UTCDateTime(gcmt["time"]),
    latitude = gcmt["latitude"],
    longitude = gcmt["longitude"],
    depth = gcmt["depth_in_meters"],
    origin_type = "centroid")
    # Tensor elements
    tensor = Tensor(m_rr = gcmt["mrr"],
                  m_tt = gcmt["mtt"],
                  m_pp = gcmt["mpp"],
                  m_rt = gcmt["mrt"],
                  m_rp = gcmt["mrp"],
                  m_tp = gcmt["mtp"])
    # Scalar seismic moment
    m0 = 1.0 / sqrt(2.0) * sqrt(
       gcmt["mrr"] ** 2 + 
       gcmt["mtt"] ** 2 + 
       gcmt["mpp"] ** 2 +
       2.0 * gcmt["mrt"] ** 2 + 
       2.0 * gcmt["mrp"] ** 2 + 
       2.0 * gcmt["mtp"] ** 2)

    mt = MomentTensor(tensor=tensor,scalar_moment=m0)
    foc = FocalMechanism()
    foc.moment_tensor = mt

    return (origin, foc)


# +
# First search GCMT catalog, based on simple table from Mondaic with my addition of UTCDateTime time (origin time) 
gcmtfile = '/Users/claired/Desktop/salvus_data_sets/gcmt_catalog.csv'
df = pd.read_csv(gcmtfile, index_col=0, sep=',')

for index, row in df.iterrows():
    row['origin_time'] = UTCDateTime(row['origin_time'])
    
df = df.rename(columns = {'origin_time':'time'})
# -

# CANV* Region
minlat = 31.5
maxlat = 43
minlon = -125
maxlon = -114
mindepth_m = 0
maxdepth_m = 50e3
# magnitudes
minmag = 4.5
maxmag = 6.5

# Narrow search for Napa eq
mintime = UTCDateTime('2020-11-01')
maxtime = UTCDateTime('2022-03-31')

dfs = df [ (df.latitude > minlat) & 
           (df.latitude < maxlat ) &
           (df.longitude > minlon ) &
           (df.longitude < maxlon ) &
           (df.depth_in_meters > mindepth_m) &
           (df.depth_in_meters < maxdepth_m) &
           (df.moment_magnitude > minmag) &
           (df.moment_magnitude < maxmag) &
           (df.time > mintime) &
           (df.time < maxtime) ]

dfs

# ### Grab events from Salvu

# +
from salvus import namespace as sn

p = sn.Project(path="Validation")

events = p.events.list()

import numpy as np

full_utc = []

for i in np.arange(0, len(events)):
    date = events[i][-16:]
    utcdate = str(UTCDateTime(date))
    utcdate = utcdate[:-11]
    full_utc.append(utcdate)

print(full_utc)
# -

# Update event metadata to include gcmt solution
# gcmt = dfs.to_dict('records')[0]
# origin, foc = get_gcmt(gcmt)
# evt = cat[0]
# evt.origins.clear()
# evt.origins.append(origin)
# evt.focal_mechanisms.append(foc)
# evt.preferred_focal_mechanism_id = foc.resource_id.id

# +
catalog = Catalog()
for index, row in dfs.iterrows():
    origin, foc = get_gcmt(row.to_dict())
    mag = Magnitude(resource_id=origin.resource_id,
                    magnitude_type='Mw',
                    mag=row.moment_magnitude)
    event = Event(resource_id=origin.resource_id,
                  event_type='earthquake',
                  origins=[origin],
                  focal_mechanisms=[foc],
                  magnitudes=[mag])
    #event.preferred_origin_id = event.origins[0].resource_id
    catalog.append(event)
    
print('Search returned ',len(catalog),' events')
# -

print('Mw from table: ', dfs.iloc[0].moment_magnitude)
print('Mw from scalar_moment: ', mw_from_m0(foc.moment_tensor.scalar_moment))

catalog.write('CANV_catalog_Validation.xml',format="QUAKEML")








