#! /usr/bin/python2
import sys
sys.path.append('../Beta_Pyno')
import pynoramix as pyno

if 1:
    #apo=pyno.cl_set('2WC2.pdb')
    #prot=pyno.cl_set('1CGP.pdb')
    #apo.delete_coors( frame = range(1,apo.num_frames) )

    #prot=pyno.cl_set('pdb/1g6n.pdb')
    prot=pyno.molecule('pdb/1g6n.pdb')

##    # selections
##    lidx=[]
##    for ii in prot.atom:
##        if ii.name == 'CA' and 10<ii.resid.pdb_index <120:
##            print ii.index
##            lidx.append(ii.index)
##
##    ca=prot.selection(lidx)
##    #apoca=apo.selection( 'atom_name (CA)')
##    #apoca=apo.selection( 'backbone')
##
##    #apoca.write_pdb('cap-ca-10-120.pdb')
##    #apoca.write_pdb('apo-ca-10-120.pdb')
##
##if 1:
##    prot=pyno.cl_set('apo-ca-10-120.pdb')
##    nm=pyno.anm( prot, cutoff=8 )
##    pyno.build_fluct_anm( prot , nm , mode=7, output='tmp' )
##
