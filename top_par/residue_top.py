#residue_atoms={'amino':[A,B,C],...}
#covalent_bonds={'amino':[[A,B],[A,C]],...}

residue_atoms={}
covalent_bonds={}

## ACE:
residue_atoms['ACE']=[
'atCH3',
'atHH31',
'atHH32',
'atHH33',
'atC',
'atO'
]

covalent_bonds['ACE']=[
['atCH3' , 'atHH31' ],
['atCH3' , 'atHH32' ],
['atCH3' , 'atHH33' ],
['atCH3' , 'atC'    ],
['atC'   , 'atO'    ] 
]

## ALA:
residue_atoms['ALA']=[
'atN',
'atH',
'atCA',
'atHA',
'atCB',
'atHB1',
'atHB2',
'atHB3',
'atC',
'atO'
]
 
covalent_bonds['ALA']=[
['atN'   ,'atH'    ],
['atN'   ,'atCA'   ],
['atCA'  ,'atHA'   ],
['atCA'  ,'atCB'   ],
['atCA'  ,'atC'    ],
['atCB'  ,'atHB1'  ],
['atCB'  ,'atHB2'  ],
['atCB'  ,'atHB3'  ],
['atC'   ,'atO'    ] 
]

## ARG:
residue_atoms['ARG']=[
'atN',
'atH',
'atCA',
'atHA',
'atCB',
'atHB1', 
'atHB2', 
'atCG',
'atHG1', 
'atHG2', 
'atCD',
'atHD1', 
'atHD2', 
'atNE',
'atHE',
'atCZ',
'atNH1',  
'atHH11', 
'atHH12', 
'atNH2',  
'atHH21', 
'atHH22', 
'atC',
'atO'
]

covalent_bonds['ARG']=[
['atN'   ,'atH'    ], 
['atN'   ,'atCA'   ], 
['atCA'  ,'atHA'   ], 
['atCA'  ,'atCB'   ], 
['atCA'  ,'atC'    ], 
['atCB'  ,'atHB1'  ], 
['atCB'  ,'atHB2'  ], 
['atCB'  ,'atCG'   ], 
['atCG'  ,'atHG1'  ], 
['atCG'  ,'atHG2'  ], 
['atCG'  ,'atCD'   ], 
['atCD'  ,'atHD1'  ], 
['atCD'  ,'atHD2'  ], 
['atCD'  ,'atNE'   ], 
['atNE'  ,'atHE'   ], 
['atNE'  ,'atCZ'   ], 
['atCZ'  ,'atNH1'  ], 
['atCZ'  ,'atNH2'  ], 
['atNH1' ,'atHH11' ], 
['atNH1' ,'atHH12' ], 
['atNH2' ,'atHH21' ], 
['atNH2' ,'atHH22' ], 
['atC'   ,'atO'    ] 
]

## ASN:

residue_atoms['ASN']=[
'atN',
'atH',
'atCA',
'atHA',
'atCB',
'atHB1', 
'atHB2', 
'atCG',
'atOD1',  
'atND2',  
'atHD21', 
'atHD22', 
'atC',
'atO'
]

covalent_bonds['ASN']=[
['atN'   ,'atH'     ], 
['atN'   ,'atCA'    ], 
['atCA'  ,'atHA'    ], 
['atCA'  ,'atCB'    ], 
['atCA'  ,'atC'     ], 
['atCB'  ,'atHB1'   ], 
['atCB'  ,'atHB2'   ], 
['atCB'  ,'atCG'    ], 
['atCG'  ,'atOD1'   ], 
['atCG'  ,'atND2'   ], 
['atND2' ,'atHD21'  ], 
['atND2' ,'atHD22'  ], 
['atC'   ,'atO'     ] 
]

## ASP:

residue_atoms['ASP']=[
'atN',
'atH',
'atCA',
'atHA',
'atCB',
'atHB1', 
'atHB2', 
'atCG',
'atOD1', 
'atOD2', 
'atC',
'atO'
]

covalent_bonds['ASP']=[
['atN'   ,'atH'    ] 
['atN'   ,'atCA'   ] 
['atCA'  ,'atHA'   ] 
['atCA'  ,'atCB'   ] 
['atCA'  ,'atC'    ] 
['atCB'  ,'atHB1'  ] 
['atCB'  ,'atHB2'  ] 
['atCB'  ,'atCG'   ] 
['atCG'  ,'atOD1'  ] 
['atCG'  ,'atOD2'  ] 
['atC'   ,'atO'    ] 
]

## GLN:
residue_atoms['GLN']=[
'atN',
'atH',
'atCA',
'atHA',
'atCB',
'atHB1', 
'atHB2', 
'atCG',
'atHG1', 
'atHG2', 
'atCD',
'atOE1', 
'atNE2', 
'atHE21',
'atHE22',
'atC',
'atO'
]

covalent_bonds['GLN']=[
['atN'   , 'atH'   ], 
['atN'   , 'atCA'  ], 
['atCA'  , 'atHA'  ], 
['atCA'  , 'atCB'  ], 
['atCA'  , 'atC'   ], 
['atCB'  , 'atHB1' ], 
['atCB'  , 'atHB2' ], 
['atCB'  , 'atCG'  ], 
['atCG'  , 'atHG1' ], 
['atCG'  , 'atHG2' ], 
['atCG'  , 'atCD'  ], 
['atCD'  , 'atOE1' ], 
['atCD'  , 'atNE2' ], 
['atNE2' , 'atHE21'], 
['atNE2' , 'atHE22'], 
['atC'   , 'atO'   ] 
]

## GLY:
residue_atoms['GLY']=[
'atN',
'atH',
'atCA',
'atHA1',
'atHA2',
'atC',
'atO'
]

covalent_bonds['GLY']=[
['atN'   , 'atH'   ],
['atN'   , 'atCA'  ],
['atCA'  , 'atHA1' ],
['atCA'  , 'atHA2' ],
['atCA'  , 'atC'   ],
['atC'   , 'atO'   ] 
]

## NME:
residue_atoms['NME']=[
'atN',
'atH',
'atCH3',
'atHH31',
'atHH32',
'atHH33'
]

covalent_bonds['NME']=[
['atN'   ,'atH'    ],
['atN'   ,'atCH3'  ],
['atCH3' ,'atHH31' ],
['atCH3' ,'atHH32' ],
['atCH3' ,'atHH33' ] 
]

## NAC:
residue_atoms['NAC']=[
'atN',
'atH',
'atCH3',
'atHH31',
'atHH32',
'atHH33'
]

covalent_bonds['NAC']=[
['atN'   ,'atH'    ],
['atN'   ,'atCH3'  ],
['atCH3' ,'atHH31' ],
['atCH3' ,'atHH32' ],
['atCH3' ,'atHH33' ] 
]


##### WATER:

residue_atoms['SOL']=[
'atOW',
'atHW1',
'atHW2'
]

covalent_bonds['SOL']=[
['atOW'  ,'atHW1' ],
['atOW'  ,'atHW2' ]
]

residue_atoms['HO4']=[
'atOW',
'atHW1',
'atHW2'
]

covalent_bonds['HO4']=[
['atOW'  ,'atHW1' ],
['atOW'  ,'atHW2' ]
]




##### IONS:

## NA:
residue_atoms['NA']=[
'itNA'
]

covalent_bonds['NA']=[
]

## NA:
residue_atoms['K']=[
'itK'
]

covalent_bonds['K']=[
]


##### LINKERS DIFFERENT RESIDUES:
inter_residues={
'C':'N'
}



