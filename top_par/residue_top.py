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

## NAC:
residue_atoms['SOL']=[
'atOW',
'atHW1',
'atHW2'
]

covalent_bonds['SOL']=[
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



