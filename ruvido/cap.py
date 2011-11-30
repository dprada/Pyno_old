#! /usr/bin/python
import sys
sys.path.append('../')
import pynoramix as pyno
import cPickle as pickle

db = []


#---------[ CAP - 2CAMP ]-----------------
if 1:
    pdb='cap-camp-1g6n'
    prot   = pyno.molecule('pdb/%s.pdb'%pdb)
    calist = [ ii.index  for ii in prot.atom if ii.name == 'CA' and 10<ii.resid.pdb_index<128 ]
    ca     = prot.selection(calist)
    nm     = pyno.anm( ca, cutoff=8 )

    calist += [ ii.index  for ii in prot.atom if ii.name == 'C6' ]
    ca2    = prot.selection(calist)
    nm2    = pyno.anm( ca2, cutoff=8 )

    data   = [ prot , ca , nm, ca2, nm2 ]
    db.append( data )

#---------[ CAP - APO   ]-----------------

    pdb='cap-apo-2wc2-fit'
    prot   = pyno.molecule('pdb/%s.pdb'%pdb)
    calist = [ ii.index  for ii in prot.atom if ii.name == 'CA' and 10<ii.resid.pdb_index<128 ]
    ca     = prot.selection(calist)
    nm     = pyno.anm( ca, cutoff=8 )

    calist += [ ii.index  for ii in prot.atom if ii.name == 'C6' ]
    ca2    = prot.selection(calist)
    nm2    = pyno.anm( ca2, cutoff=8 )

    data   = [ prot , ca , nm, ca2, nm2 ]
    db.append( data )

#-----------------------------------------
    pickle.dump(db,open('cap.db','w'))

else: db=pickle.load(open('cap.db'))



#--[ ANALYSIS ]---------------------------

# FLUCTUATIONS
if 0:
    for jj in [2,4]:
        for ii in [0,1]:
            #bf=db[ii][jj].bfacts
            nm=db[ii][jj]
            nm.build_correl(modes=7)
            bf=[ nm.correl[kk][kk] for kk in range(len(nm.correl))  ]
            for kk in range(len(bf)):
                print "%10d %10.4f"%(kk,bf[kk])
            print "\n\n"

# DISPLACEMENT VECTORS
if 1:
    if 0:

        prot= db[0][0]
        ca  = db[0][1]
        nm  = db[0][2]
        ref = db[1][1]

        ca.rms_fit(set_reference=ref )
        ca.displ_vector( ref )
        nm.involv_coefficient(modes=range(200),vect=ca.d_vector)

        print '\n\n'
        print '## cap -> apo (no camp molecules)'
        icsum=0.0
        for ii in nm.ic:
            icsum+=nm.ic[ii]*nm.ic[ii]
            print ii, abs(nm.ic[ii]), icsum

    if 0:

        prot= db[1][0]
        ca  = db[1][1]
        nm  = db[1][2]
        ref = db[0][1]

        ca.rms_fit(set_reference=ref )
        ca.displ_vector( ref )
        nm.involv_coefficient(modes=range(200),vect=ca.d_vector)

        print '\n\n'
        print '## apo -> camp (no camp molecules)'
        icsum=0.0
        for ii in nm.ic:
            icsum+=nm.ic[ii]*nm.ic[ii]
            print ii, abs(nm.ic[ii]), icsum


    # ---- CAMP bound -------

    if 0:

        prot= db[0][0]
        ca  = db[0][3] ### !!
        nm  = db[0][4] ### !!
        ref = db[1][3] ### !!

        ca.rms_fit(set_reference=ref )
        ca.displ_vector( ref )
        nm.involv_coefficient(modes=range(200),vect=ca.d_vector)

        print '\n\n'
        print '## cap -> apo (camp)'
        icsum=0.0
        for ii in nm.ic:
            icsum+=nm.ic[ii]*nm.ic[ii]
            print ii, abs(nm.ic[ii]), icsum

    if 1:

        prot= db[1][0]
        ca  = db[1][3] ### !!
        nm  = db[1][4] ### !!
        ref = db[0][3] ### !!

        ca.rms_fit(set_reference=ref )
        ca.displ_vector( ref )
        nm.involv_coefficient(modes=range(200),vect=ca.d_vector)

        print '\n\n'
        print '## apo -> camp (camp)'
        icsum=0.0
        for ii in nm.ic:
            icsum+=nm.ic[ii]*nm.ic[ii]
            print ii, abs(nm.ic[ii]), icsum

        pyno.build_fluct_anm( prot, nm, mode=7,amplitude=8.0,steps=60)
