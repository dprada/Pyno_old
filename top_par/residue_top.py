#residue_atoms={'amino':[A,B,C],...}
#covalent_bonds={'amino':[[A,B],[A,C]],...}

residue_atoms={}
covalent_bonds={}

### ACE:

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

### ALA:

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

### ARG:

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

### ASN:

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

### ASP:

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
['atN'   ,'atH'    ], 
['atN'   ,'atCA'   ], 
['atCA'  ,'atHA'   ], 
['atCA'  ,'atCB'   ], 
['atCA'  ,'atC'    ], 
['atCB'  ,'atHB1'  ], 
['atCB'  ,'atHB2'  ], 
['atCB'  ,'atCG'   ], 
['atCG'  ,'atOD1'  ], 
['atCG'  ,'atOD2'  ], 
['atC'   ,'atO'    ] 
]

### GLU:

residue_atoms['GLU']=[
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
'atOE2', 
'atC',
'atO'
]

covalent_bonds['GLU']=[
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
['atCD'  , 'atOE2' ], 
['atC'   , 'atO'   ] 
]


### GLN:

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

### GLY:

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

### HISE (ND1 no H, NE2 with H),
### HISD (ND1 with H, NE2 no H),
### HISH (ND1 with H, NE2 with H),
### All included in HIS:

residue_atoms['HIS']=[
'atN',
'atH',
'atCA',
'atHA',
'atCB',
'atHB1', 
'atHB2',
'atCG',
'atND1', 
'atHD1',
'atCD2', 
'atHD2', 
'atCE1', 
'atHE1', 
'atNE2', 
'atHE2', 
'atC',
'atO'
]

covalent_bonds['HIS']=[
['atN'   ,'atH'    ], 
['atN'   ,'atCA'   ], 
['atCA'  ,'atHA'   ], 
['atCA'  ,'atCB'   ], 
['atCA'  ,'atC'    ], 
['atCB'  ,'atHB1'  ], 
['atCB'  ,'atHB2'  ], 
['atCB'  ,'atCG'   ], 
['atCG'  ,'atND1'  ], 
['atCG'  ,'atCD2'  ], 
['atND1' ,'atHD1'  ], 
['atND1' ,'atCE1'  ], 
['atCD2' ,'atHD2'  ], 
['atCD2' ,'atNE2'  ], 
['atCE1' ,'atHE1'  ], 
['atCE1' ,'atNE2'  ], 
['atNE2' ,'atHE2'  ], 
['atC'   ,'atO'    ] 
]


### ILE:

residue_atoms['ILE']=[
'atN', 
'atH', 
'atCA', 
'atHA', 
'atCB', 
'atHB', 
'atCG1', 
'atHG11', 
'atHG12', 
'atCG2', 
'atHG21', 
'atHG22', 
'atHG23', 
'atCD', 
'atHD1', 
'atHD2', 
'atHD3', 
'atC', 
'atO' 
]

covalent_bonds['ILE']=[
['atN',   'atH'    ], 
['atN',   'atCA'   ], 
['atCA',  'atHA'   ], 
['atCA',  'atCB'   ], 
['atCA',  'atC'    ], 
['atCB',  'atHB'   ], 
['atCB',  'atCG1'  ], 
['atCB',  'atCG2'  ], 
['atCG1', 'atHG11' ], 
['atCG1', 'atHG12' ], 
['atCG1', 'atCD'   ], 
['atCG2', 'atHG21' ], 
['atCG2', 'atHG22' ], 
['atCG2', 'atHG23' ], 
['atCD',  'atHD1'  ], 
['atCD',  'atHD2'  ], 
['atCD',  'atHD3'  ], 
['atC',   'atO'    ] 
]

### LEU:

residue_atoms['LEU']=[
'atN',
'atH',
'atCA',
'atHA',
'atCB',
'atHB1',
'atHB2',
'atCG',
'atHG',
'atCD1',
'atHD11', 
'atHD12', 
'atHD13', 
'atCD2',
'atHD21', 
'atHD22', 
'atHD23', 
'atC',
'atO'
]

covalent_bonds['LEU']=[
['atN',   'atH'    ], 
['atN',   'atCA'   ], 
['atCA',  'atHA'   ], 
['atCA',  'atCB'   ], 
['atCA',  'atC'    ], 
['atCB',  'atHB1'  ], 
['atCB',  'atHB2'  ], 
['atCB',  'atCG'   ], 
['atCG',  'atHG'   ], 
['atCG',  'atCD1'  ], 
['atCG',  'atCD2'  ], 
['atCD1', 'atHD11' ], 
['atCD1', 'atHD12' ], 
['atCD1', 'atHD13' ], 
['atCD2', 'atHD21' ], 
['atCD2', 'atHD22' ], 
['atCD2', 'atHD23' ], 
['atC',   'atO'    ] 
]

### LYS:

residue_atoms['LYS']=[
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
'atCE', 
'atHE1',
'atHE2',
'atNZ', 
'atHZ1',
'atHZ2',
'atC',
'atO'
]

covalent_bonds['LYS']=[
['atN',   'atH'   ], 
['atN',   'atCA'  ], 
['atCA',  'atHA'  ], 
['atCA',  'atCB'  ], 
['atCA',  'atC'   ], 
['atCB',  'atHB1' ], 
['atCB',  'atHB2' ], 
['atCB',  'atCG'  ], 
['atCG',  'atHG1' ], 
['atCG',  'atHG2' ], 
['atCG',  'atCD'  ], 
['atCD',  'atHD1' ], 
['atCD',  'atHD2' ], 
['atCD',  'atCE'  ], 
['atCE',  'atHE1' ], 
['atCE',  'atHE2' ], 
['atCE',  'atNZ'  ], 
['atNZ',  'atHZ1' ], 
['atNZ',  'atHZ2' ], 
['atC',   'atO'   ] 
]

### MET:

residue_atoms['MET']=[
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
'atSD', 
'atCE', 
'atHE1',
'atHE2',
'atHE3',
'atC',
'atO'
]

covalent_bonds['MET']=[
['atN',   'atH'    ], 
['atN',   'atCA'   ], 
['atCA',  'atHA'   ], 
['atCA',  'atCB'   ], 
['atCA',  'atC'    ], 
['atCB',  'atHB1'  ], 
['atCB',  'atHB2'  ], 
['atCB',  'atCG'   ], 
['atCG',  'atHG1'  ], 
['atCG',  'atHG2'  ], 
['atCG',  'atSD'   ], 
['atSD',  'atCE'   ], 
['atCE',  'atHE1'  ], 
['atCE',  'atHE2'  ], 
['atCE',  'atHE3'  ], 
['atC',   'atO'    ] 
]

### NME:

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

### NAC:

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

### NHE:

residue_atoms['NHE']=[
'atN',
'atH1',
'atH2'
]

covalent_bonds['NHE']=[
['atN'   ,'atH1'   ],
['atN'   ,'atH2'   ]
]

### PHE:

residue_atoms['PHE']=[
'atN',
'atH',
'atCA',
'atHA',
'atCB',
'atHB1', 
'atHB2', 
'atCG',
'atCD1', 
'atHD1', 
'atCD2', 
'atHD2', 
'atCE1', 
'atHE1', 
'atCE2', 
'atHE2', 
'atCZ',
'atHZ',
'atC',
'atO'
]

covalent_bonds['PHE']=[
['atN',   'atH'    ], 
['atN',   'atCA'   ], 
['atCA',  'atHA'   ], 
['atCA',  'atCB'   ], 
['atCA',  'atC'    ], 
['atCB',  'atHB1'  ], 
['atCB',  'atHB2'  ], 
['atCB',  'atCG'   ], 
['atCG',  'atCD1'  ], 
['atCG',  'atCD2'  ], 
['atCD1', 'atHD1'  ], 
['atCD1', 'atCE1'  ], 
['atCD2', 'atHD2'  ], 
['atCD2', 'atCE2'  ], 
['atCE1', 'atHE1'  ], 
['atCE1', 'atCZ'   ], 
['atCE2', 'atHE2'  ], 
['atCE2', 'atCZ'   ], 
['atCZ',  'atHZ'   ], 
['atC',   'atO'    ] 
]

### PRO:

residue_atoms['PRO']=[
'atN',
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
'atC', 
'atO'
]

covalent_bonds['PRO']=[
['atN',   'atCA'   ], 
['atCA',  'atHA'   ], 
['atCA',  'atCB'   ], 
['atCA',  'atC'    ], 
['atCB',  'atHB1'  ], 
['atCB',  'atHB2'  ], 
['atCB',  'atCG'   ], 
['atCG',  'atHG1'  ], 
['atCG',  'atHG2'  ], 
['atCG',  'atCD'   ], 
['atCD',  'atHD1'  ], 
['atCD',  'atHD2'  ], 
['atCD',  'atN'    ], 
['atC',   'atO'    ] 
]

### SER:

residue_atoms['SER']=[
'atN',
'atH',
'atCA', 
'atHA', 
'atCB', 
'atHB1',
'atHB2',
'atOG', 
'atHG', 
'atC',
'atO'
]

covalent_bonds['SER']=[
['atN',  'atH'    ], 
['atN',  'atCA'   ], 
['atCA', 'atHA'   ], 
['atCA', 'atCB'   ], 
['atCA', 'atC'    ], 
['atCB', 'atHB1'  ], 
['atCB', 'atHB2'  ], 
['atCB', 'atOG'   ], 
['atOG', 'atHG'   ], 
['atC',  'atO'    ] 
]

### THR:

residue_atoms['THR']=[
'atN',
'atH',
'atCA',
'atHA',
'atCB',
'atHB',
'atOG1',
'atHG1',
'atCG2',
'atHG21',
'atHG22', 
'atHG23', 
'atC',
'atO'
]

covalent_bonds['THR']=[
['atN',   'atH'     ], 
['atN',   'atCA'    ], 
['atCA',  'atHA'    ], 
['atCA',  'atCB'    ], 
['atCA',  'atC'     ], 
['atCB',  'atHB'    ], 
['atCB',  'atOG1'   ], 
['atCB',  'atCG2'   ], 
['atOG1', 'atHG1'   ], 
['atCG2', 'atHG21'  ], 
['atCG2', 'atHG22'  ], 
['atCG2', 'atHG23'  ], 
['atC',   'atO'     ] 
]

### TRP:

residue_atoms['TRP']=[
'atN',
'atH',
'atCA',
'atHA',
'atCB',
'atHB1', 
'atHB2', 
'atCG',
'atCD1', 
'atHD1', 
'atCD2', 
'atNE1', 
'atHE1', 
'atCE2', 
'atCE3', 
'atHE3', 
'atCZ2', 
'atHZ2', 
'atCZ3', 
'atHZ3', 
'atCH2', 
'atHH2', 
'atC',
'atO'
]

covalent_bonds['TRP']=[
['atN',   'atH'   ], 
['atN',   'atCA'  ], 
['atCA',  'atHA'  ], 
['atCA',  'atCB'  ], 
['atCA',  'atC'   ], 
['atCB',  'atHB1' ], 
['atCB',  'atHB2' ], 
['atCB',  'atCG'  ], 
['atCG',  'atCD1' ], 
['atCG',  'atCD2' ], 
['atCD1', 'atHD1' ], 
['atCD1', 'atNE1' ], 
['atCD2', 'atCE2' ], 
['atCD2', 'atCE3' ], 
['atNE1', 'atHE1' ], 
['atNE1', 'atCE2' ], 
['atCE2', 'atCZ2' ], 
['atCE3', 'atHE3' ], 
['atCE3', 'atCZ3' ], 
['atCZ2', 'atHZ2' ], 
['atCZ2', 'atCH2' ], 
['atCZ3', 'atHZ3' ], 
['atCZ3', 'atCH2' ], 
['atCH2', 'atHH2' ], 
['atC',   'atO'   ] 
]

### TYR:

residue_atoms['TYR']=[
'atN',
'atH',
'atCA', 
'atHA', 
'atCB', 
'atHB1',
'atHB2',
'atCG', 
'atCD1',
'atHD1', 
'atCD2', 
'atHD2', 
'atCE1', 
'atHE1', 
'atCE2', 
'atHE2', 
'atCZ', 
'atOH', 
'atHH', 
'atC',
'atO'
]

covalent_bonds['TYR']=[
['atN',   'atH'    ], 
['atN',   'atCA'   ], 
['atCA',  'atHA'   ], 
['atCA',  'atCB'   ], 
['atCA',  'atC'    ], 
['atCB',  'atHB1'  ], 
['atCB',  'atHB2'  ], 
['atCB',  'atCG'   ], 
['atCG',  'atCD1'  ], 
['atCG',  'atCD2'  ], 
['atCD1', 'atHD1'  ], 
['atCD1', 'atCE1'  ], 
['atCD2', 'atHD2'  ], 
['atCD2', 'atCE2'  ], 
['atCE1', 'atHE1'  ], 
['atCE1', 'atCZ'   ], 
['atCE2', 'atHE2'  ], 
['atCE2', 'atCZ'   ], 
['atCZ',  'atOH'   ], 
['atOH',  'atHH'   ], 
['atC',   'atO'    ] 
]

### VAL:

residue_atoms['VAL']=[
'atN',
'atH', 
'atCA',
'atHA',
'atCB',
'atHB',
'atCG1', 
'atHG11',
'atHG12',
'atHG13',
'atCG2', 
'atHG21',
'atHG22',
'atHG23',
'atC',
'atO'
]

covalent_bonds['VAL']=[
['atN',   'atH'     ], 
['atN',   'atCA'    ], 
['atCA',  'atHA'    ], 
['atCA',  'atCB'    ], 
['atCA',  'atC'     ], 
['atCB',  'atHB'    ], 
['atCB',  'atCG1'   ], 
['atCB',  'atCG2'   ], 
['atCG1', 'atHG11'  ], 
['atCG1', 'atHG12'  ], 
['atCG1', 'atHG13'  ], 
['atCG2', 'atHG21'  ], 
['atCG2', 'atHG22'  ], 
['atCG2', 'atHG23'  ], 
['atC',   'atO'     ] 
]

##### WATER:

### SOL:

residue_atoms['SOL']=[
'atOW',
'atHW1',
'atHW2'
]

covalent_bonds['SOL']=[
['atOW'  ,'atHW1' ],
['atOW'  ,'atHW2' ]
]

### HOH:

residue_atoms['HOH']=[
'atOW',
'atHW1',
'atHW2'
]

covalent_bonds['HOH']=[
['atOW'  ,'atHW1' ],
['atOW'  ,'atHW2' ]
]

### HO4:

residue_atoms['HO4']=[
'atOW',
'atHW1',
'atHW2',
'atvir'
]

covalent_bonds['HO4']=[
['atOW'  ,'atHW1' ],
['atOW'  ,'atHW2' ]
]

### HO5:

residue_atoms['HO5']=[
'atOW',
'atHW1',
'atHW2',
'atvir',
'atvir'
]

covalent_bonds['HO5']=[
['atOW'  ,'atHW1' ],
['atOW'  ,'atHW2' ]
]


### SWM:

residue_atoms['SWM']=[
'atOW',
'atHW1',
'atHW2',
'atvir',
'atvir'
]

covalent_bonds['SWM']=[
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

## K:
residue_atoms['K']=[
'itK'
]

covalent_bonds['K']=[
]

## LI:
residue_atoms['LI']=[
'itLI'
]

covalent_bonds['LI']=[
]

## CL:
residue_atoms['CL']=[
'itCL'
]

covalent_bonds['CL']=[
]




##### LINKERS DIFFERENT RESIDUES:
inter_residues={
'C':'N'
}
