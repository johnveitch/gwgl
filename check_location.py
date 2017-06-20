#!/usr/bin/env python

# check_location.py: Check for lenses in a GW sky map
# Copyright (C) 2017  John Veitch <john.veitch@ligo.org>

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from __future__ import division
import healpy as hp
import numpy as np
import astropy
from astropy.coordinates import SkyCoord
from astropy import units as u
from optparse import OptionParser
import sys

def load_lenses(file):
    """
    Read in the list of lenses from a file formatted like this:
    ABELL2744            00 14 18.9 -30 23 22
    """
    clusters=[]
    for line in file:
        name, ra_h,ra_m,ra_s, dec_deg,dec_m,dec_s = line.split()
        c = SkyCoord('{0}h{1}m{2}s'.format(ra_h,ra_m,ra_s),
                     '{0}d{1}m{2}s'.format(dec_deg,dec_m,dec_s),
                     frame='icrs')
        clusters.append((name,c))
    return clusters


def check_ci(skymap,clusters,ci=0.9,verbose=False):
    """
    Look for clusters within the given credible interval
    """
    nside=hp.get_nside(skymap)
    thetas,ras = hp.pixelfunc.pix2ang(nside,range(len(skymap)))
    decs = np.pi/2 - thetas

    # Sort the map
    imap=np.argsort(skymap)[::-1]
    idxmap=np.array(range(len(imap)))
    idxlist = idxmap[imap].tolist()
    Ps=np.cumsum(skymap[imap])
 
    skycoords = SkyCoord(ra=ras*u.radian,dec=decs*u.radian,frame='icrs')
    
    found = []
 
    for name,c in clusters:
        idx, d2d, d3d = c.match_to_catalog_sky(skycoords)
    
        n=idxlist.index(idx)
        P=Ps[n]
        if P<ci: found.append((name,c,skycoords[idx],P))
        if verbose: print('Cluster {0} found at P {1}'.format(name,P))
        
    return found

usage=""" %prog [options] skymap.fits.gz
Find lenses within the sky map"""

if __name__=='__main__':
    
    parser=OptionParser(usage)
    parser.add_option('-l','--lens-file',default='strong-lensing-clusters-20170414.txt',help='File with clusters')
    parser.add_option('-P','--credible-interval',default=0.9,type=float,metavar='P',help='Credible interval to search within')
    parser.add_option('-v','--verbose',default=False,action='store_true')
    (opts,args)=parser.parse_args()
    if(len(args)!=1):
        parser.print_help()
        sys.exit(1)
    fitsfile = args[0]
    
    clusters = load_lenses(open(opts.lens_file))

    map=hp.read_map(fitsfile,h=False,verbose=opts.verbose)
    found=check_ci(map,clusters,ci=opts.credible_interval,verbose=opts.verbose)
    
    print('The following clusters were found inside the 90% c.i.')
    for n,c,_,P in found:
        print('{0}:\t{1}\tp={2}'.format(n,c.to_string('hmsdms'),P))
    

