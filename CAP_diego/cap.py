#! /usr/bin/python2
import sys
sys.path.append('../')
from pynoramix import *



prot   = molecule('cap-apo-2wc2-fit.pdb')
calist = [ ii.index  for ii in prot.atom if ii.name == 'CA' and 10<ii.resid.pdb_index<128 ]
calist += [ ii.index  for ii in prot.atom if ii.name == 'C6' ]
CAs    = prot.selection(calist)
nm    = anm( CAs, cutoff=8 )


build_fluct_anm( prot, nm, mode=7,amplitude=8.0,steps=60)
