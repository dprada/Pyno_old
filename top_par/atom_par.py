#################################
######## Charge

charge={}

charge={
'itNA': +1.0,
'itK': +1.0,
'itCL': -1.0
}


#################################
######## Donor - Acceptor

### Donors

donors=[
'atN'   ,   # N for the main chain
'atNE'  ,   # N for ARG
'atNH1' ,   # N for ARG
'atNH2' ,   # N for ARG
'atNE1' ,   # N for TRP
'atNE2' ,   # N for HIS and GLN (presence of H determines if donor)
'atND1' ,   # N for HIS (presence of H determines if donor)
'atND2' ,   # N for ASN
'atNZ'  ,   # N for LYS
'atOG'  ,   # O for SER
'atOH'  ,   # O for TYR
'atOW'  ,   # O for WAT
'atOG1' ,   # O for THR
'atSG'      # S for CYS
]

donors_exception=[
['HIS'  , 'atNE2'  , 'Not Hbonded'  , False],
['HIS'  , 'atND1'  , 'Not Hbonded'  , False]
]

### Acceptors

acceptors=[
'atND1' ,   # N for HIS (presence of H determines if donor-acceptor)
'atNE2' ,   # N for HIS (presence of H determines if donor-acceptor) (Not in GLN)
'atO'   ,   # O for the main chain
'atOD1' ,   # O for ASN and ASP
'atOD2' ,   # O for ASP
'atOE1' ,   # O for GLN and GLU
'atOE2' ,   # O for GLU
'atOG'  ,   # O for SER
'atOG1' ,   # O for THR
'atOH'  ,   # O for TYR
'atOW'  ,   # O for WAT
'atSD'      # S for MET
]

acceptors_exception=[
['GLN'  , 'atNE2'  , 'Always'       , False],
['CYH'  , 'atSG'   , 'Always'       , True ],
['HIS'  , 'atNE2'  , 'Hbonded'      , False],
['HIS'  , 'atND1'  , 'Hbonded'      , False]
]
