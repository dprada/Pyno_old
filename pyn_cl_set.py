############ GENERAL COMMENTS
#
# This module requires:
#
#from pyn_cl_unit import *
#from pyn_cl_coors import *
#from os import system
#from os import path
#from numpy import *
#import pylab
#import pyn_fort as f
#import pickle as pic
#
# Structure of the file:
#
# -Class Set
#            - Instantiation
#            - Functions
# -Class Residue
# -Class Water
# -External functions

class cl_set:

    
    ###
    ### Class Instantiation
    ###

    def __init__(self,input_file=None,input_selection=None,download=None,verbose=True):

        # > Instantation options:
        self.file=input_file
        self.select=input_selection
        self.file_hbonds=''
        self.file_mss=''
        self.file_shell=''
        

        # > Topological properties
        self.num_atoms=0
        self.name=''
        self.index=0
        self.pdb_index=0
        self.atom=[]       #list of objects unit
        self.list_atoms=[]
        self.acceptors=[]
        self.donors=[]
        self.dimensionality=0
        self.ss_pdb=[]
        self.chains=[]
        self.water_model=''

        # > Coordinates
        self.recent_frame=0
        self.coors=[]      # list of objects coor
        self.dist_matrix=[]

        ##################################

        # BUILDING THE TOPOLOGY OF THE SET FROM A FILE

        ### Download option:

        if download:
            if not path.exists(download):
                temp='wget -nv http://www.rcsb.org/pdb/files/'+download
                system(temp)

            input_file=download
            self.file=input_file

        ### Reading the file and attaching the atoms to the set

        if self.file:

            if self.file.endswith('pdb'):
                self.read_pdb(self.file)

            if self.file.endswith('gro'):
                self.read_gro(self.file)

            ### Setting up other atoms' attributes



            ### Setting up the chains
            for aa in self.atom[:]:
                if aa.type_pdb in ['ATOM']:
                    if aa.chain not in self.chains:
                        self.chains.append(aa.chain)


            ### Auxiliary dictionary
            
            before=-99
            ii=-1
            for aa in self.atom:
                if aa.resid_pdb_index != before :
                    before=aa.resid_pdb_index
                    ii+=1
                aa.resid_index=ii


            aux={}

            for aa in self.atom[:]:
                ii=aa.resid_index
                try: 
                    aux[ii][aa.name]=aa.index
                except:
                    aux[ii]={}
                    aux[ii][aa.name]=aa.index

            for aa in self.atom[:]:
                if aa.name in ['OW']:
                    aa.acceptor=True
                    aa.polar_class='acceptor'
                    aa.polarizability=True
                    aa.covalent_bond.append(aux[aa.resid_index]['HW1'])
                    aa.covalent_bond.append(aux[aa.resid_index]['HW2'])
                if aa.name in ['O'] and aa.resid_name in ['HOH','SOL','HO4','water']:
                    aa.acceptor=True
                    aa.polar_class='acceptor'
                    aa.polarizability=True
                if aa.name in ['HW1','HW2']:
                    aa.donor=True
                    aa.polar_class='donor'
                    aa.polarizability=True
                    aa.covalent_bond.append(aux[aa.resid_index]['OW'])
        
            ### Setting up the residues

            self.residue=[]
            for aa in aux.keys():
                temp_residue=cl_residue()
                temp_residue.index=aa
                temp_residue.list_atoms=aux[aa].values()
                ii=temp_residue.list_atoms[0]
                temp_residue.pdb_index=self.atom[ii].resid_pdb_index
                temp_residue.name=self.atom[ii].resid_name
                self.residue.append(temp_residue)


            ### Setting up the waters

            self.water=[]
            for aa in self.residue[:]:
                if aa.name in ['HOH','SOL','HO4']:
                    temp_water=cl_water()
                    temp_water.name=aa.name
                    temp_water.index=aa.index
                    temp_water.pdb_index=aa.pdb_index
                    temp_water.list_atoms=aa.list_atoms
                    if self.water_model==None:
                        self.water_model='tip'+str(len(aa.list_atoms))+'p'
                    if 'HW1' in aux[aa.index].keys():
                        temp_water.H1=aux[aa.index]['HW1']
                        temp_water.H2=aux[aa.index]['HW2']
                        temp_water.O=aux[aa.index]['OW']
                    if 'O' in aux[aa.index].keys():
                        temp_water.O=aux[aa.index]['O']
                    self.water.append(temp_water)

            ### Setting up global attributes

            ii=self.file[::-1].find('.')
            self.name=self.file[:-ii]
            self.num_atoms=len(self.atom)
            self.dimensionality=self.num_atoms*3
            for aa in self.atom[:]:
                if aa.acceptor: self.acceptors.append(aa.index)
                if aa.donor: self.donors.append(aa.index)
            self.acceptors=array(self.acceptors,order='Fortran')
            self.donors=array(self.donors,order='Fortran')
            self.num_residues=len(self.residue)
            self.num_waters=len(self.water)


            ### Print info:
            if self.verbose:
                self.info()



    ###
    ### Functions
    ###

    # Info function

    def info(self):

        print '#','System created from the file ',self.file,':'
        print '#',self.num_atoms,' atoms'
        print '#',self.num_residues,' residues'
        print '#',len(self.chains),' chains'
        print '#',self.num_waters,' waters'


    # To handle files

    def read_pdb (self,name_file):

        for line in open(name_file,'r'):
            ss=line.split()

            if ss[0] in ['HELIX','SHEET','TURN']:
                self.ss_pdb.append(line)
            if ss[0] in ['ATOM','HETATM']:

                temp_atom=cl_unit()
                temp_atom.type_pdb=line[0:6].replace(' ', '')
                temp_atom.pdb_index=int(line[6:11])
                temp_atom.name=(line[12:16].split())[0]
                temp_atom.alt_loc=line[16]
                temp_atom.resid_name=(line[17:20])
                temp_atom.chain=line[21]
                temp_atom.resid_pdb_index=int(line[22:26])
                temp_atom.code_ins_res=line[26]
                temp_atom.occup=float(line[54:60])
                temp_atom.bfactor=float(line[60:66])
                temp_atom.seg_ident=line[72:76].replace(' ', '')
                temp_atom.elem_symb=line[76:78].replace(' ', '')
                temp_atom.charge=line[78:80].replace(' ', '')

                temp_atom.index=len(self.atom)
                self.atom.append(temp_atom)


    def read_gro (self,name_file):

        f=open(name_file,'r')

        line=f.readline()                                          # Header of the gro file

        line=f.readline()                                        
        self.num_atoms=int(line)

        for i in range(self.num_atoms):           
            
            temp_atom=cl_unit()

            line=f.readline().split()
            temp_atom.pdb_index=int(line[2])
            temp_atom.name=line[1]
            temp_atom.resid_name=line[0][-3:]
            temp_atom.resid_pdb_index=int(line[0][:-3]) 

            temp_atom.index=i           

            self.atom.append(temp_atom)


    def write_pdb (self,filename=None):
        
        if filename==None:
            print 'Enter filename: '
            print '      foo.write_pdb("foo.pdb")'
        else:
            file=open(filename,'w')
            for ii in self.ss_pdb:
                file.write(str(ii))

            for ii in range(self.num_atoms):
                a='ATOM  '                                  # 1-6
                a+="%5d" % (ii+1)                           # 7-11
                #a+="%5d" % self.atom[ii].pdb_index          # 7-11
                a+=' '                                      # 12
                a+="%-4s" % self.atom[ii].name           # 13-16
                a+=' '                                      # 17
                a+="%3s" % self.atom[ii].resid_name          # 18-20
                a+=' '                                      # 21
                a+="%1s" % self.atom[ii].chain               # 22
                a+="%4d" % self.atom[ii].resid_pdb_index     # 23-26
                a+=' '                                      # 27
                a+='   '                                    # 28-30
                a+="%8.3f" % float(self.coors[0].xyz[ii][0]) # 31-38
                a+="%8.3f" % float(self.coors[0].xyz[ii][1]) # 39-46
                a+="%8.3f" % float(self.coors[0].xyz[ii][2]) # 47-54
                a+="%6.2f" % self.atom[ii].occup             # 55-60
                a+="%6.2f" % self.atom[ii].bfactor           # 61-66
                a+='          '                             # 67-76
                a+="%2s" % self.atom[ii].elem_symb           # 77-78
                a+="%2s" % self.atom[ii].charge              # 79-80
                a+='\n' 
                file.write(str(a))         
            file.close()
        return None

    def write_set_to_file(self,name_of_file):
        file=open(name_of_file,'w')
        pic.dump(self,file)
        file.close()

    def read_set_from_file(self,name_of_file):
        file=open(name_of_file,'r')
        A=pic.load(file)
        file.close()
        return A


    # To handle coordinates


    def load_coors (self,input_file,frame=None,begin=None,end=None):

        self.coors_file=input_file
        self.traj_mark=0

        if self.coors_file.endswith('pdb'):
            temp_frame=cl_coors(self.file)
            self.coors.append(temp_frame)

        elif self.coors_file.endswith('gro'):
            temp_frame=cl_coors(self.file)
            self.coors.append(temp_frame)

        elif self.coors_file.endswith('bin'):
            if begin==None and frame==None and end==None:
                temp_frame=cl_coors(self.coors_file,self.recent_frame)
                self.coors.append(temp_frame)
                self.recent_frame+=1
            elif begin==None and end==None and frame!=None:
                temp_frame=cl_coors(self.coors_file,frame)
                self.coors.append(temp_frame)
                self.recent_frame=frame
            elif begin!=None and end!=None:
                for ii in range(begin,end):
                    temp_frame=cl_coors(self.coors_file,ii)
                    self.coors.append(temp_frame)
                self.recent_frame=end

    def delete_coors (self,begin=None,end=None,frame=None):
        
        if frame==begin==end==None :
            del self.coors[:]


    # The functions, really...

    def distance(self,pbc=False):
                
        for frame in self.coors:

            dist_frame=f.aux_funcs_general.dist(pbc,frame.xyz,frame.box,frame.xyz,self.num_atoms,self.num_atoms)
            
        self.dist_matrix.append(dist_frame)

    def neighbs(self,system2=None,limit=0,dist=0.0,pbc=False):


        if system2==None:
            system2=self
            ident=True
        else:
            ident=False

        if limit != 0:
            neighbs=f.aux_funcs_general.neighbs_limit(pbc,ident,limit,self.coors[0].xyz,self.coors[0].box,system2.coors[0].xyz,self.num_atoms,system2.num_atoms)
        
        else:
            print type(self.coors[0].xyz[1][:]),self.coors[0].xyz[1][:],system2.num_atoms
            neighbs=f.aux_funcs_general.neighbs_dist(pbc,ident,dist,self.coors[0].xyz,self.coors[0].box,system2.coors[0].xyz,system2.num_atoms)


        return neighbs

        
    def plot_contact_map(contact_map):
        
        pylab.gray()
        pylab.imshow(contact_map==False,origin='lower',interpolation=None) # Would be better not to interpolate
        #pylab.matshow(contact_map==False)
        return pylab.show()






#######################################################
#### EXTERNAL OBJECTS AND FUNCTIONS
#######################################################


    
###
### Classes
###


class cl_water(cl_set):
    '''Water specific attributes and methods'''

    def __init__( self ,water_model=None):

        self.name=''
        self.index=''
        self.pdb_index=''
        self.list_atoms=''
        self.water_model=water_model
                
        self.O=''
        self.H1=''
        self.H2=''
        self.uvect_norm=[]
        #def get_uvect_norm(self):
        
    pass


class cl_residue(cl_set):

    def __init__( self ):

        self.index=''
        self.list_atoms=''
        self.pdb_index=''
        self.name=''

    pass


###
### Functions
###



def xtc2bin(xtc_name,bin_name):
    command=home_path+'xtc2bin %s %s'%(xtc_name,bin_name)

    if path.exists(bin_name):
        command2='mv %s %s#'%(bin_name,bin_name)
        print 'file',bin_name,' was moved to ', bin_name+'#'
        system(command2)

    system(command)



    #######################################################
    #### Selection algorithm:    
    #######################################################

def make_selection(system,condition):                 #####system - your system
						       #####condition - keywords to select,as example: 	OO=make_selection2(wbox,'atom_name (OW or MW) and resid_name (HO4)')
    if type(condition)==str:                           #####also can work with list of indexes	

        condition=prolong_string(condition)

    	AA=condition.split(' and ')                        #####finding keyword 'and'
    
    	aux_lists={}
    	Num_con=len(AA)
    	for jj in range(Num_con):
		aux_lists[jj]=good_select(system,AA[jj])      ######making selection for splitted keywords and putting the indexes in the dictionary aux_lists

    	last_list=[]					
    	ll2=[]
    	for jj in range(Num_con):			       ###### Comparing lists from different keywords,and creating final list of indexed ll2
		if jj==Num_con-1:
		 	for ii in range(len(aux_lists[jj])):
				last_list.append(aux_lists[jj][ii])
			
				if last_list.count(aux_lists[jj][ii])==Num_con:
    					ll2.append(ii)

		else:
			 for ii in range(len(aux_lists[jj])):
				last_list.append(aux_lists[jj][ii])

    	sux=appending_sel(system,ll2)

    elif type(condition) in [list,ndarray,tuple]:      ##########if you want to select from list with indexes of atoms
        LIST=condition
        sux=appending_sel(system,LIST)
    else:
        print "ERROR sel01"

    return sux



def good_select(system,condition):
    list_of_ind=[]

    AA=condition.split('(')[1][:-2]		      ##condition to illustrate "atom_name (C or OW or CA or CB)
    param=condition.split('(')[0].replace(' ','')     #param - key word like atom_name
    possib=AA.split(' or ')			      #possib - array of possibilies, like [C,OW,CA,CB] if 
    

				
    for jj in range(len(possib)):		      #making selection using auxiliary functions g_parameter; if you want to select by new parameter - you need to add your own function
	    if param=='atom_name':  
		
                list_of_ind=g_atom_name(system,possib[jj].replace(' ',''),list_of_ind) 
              
            elif param=='resid_name':
                list_of_ind=g_resid_name(system,possib[jj].replace(' ',''),list_of_ind)
            elif param=='donors':
                possib='donor'
                list_of_ind=g_donors(system,possib[jj].replace(' ',''),list_of_ind)
            elif param=='acceptors':
                possib='acceptor'
                list_of_ind=g_acceptors(system,possib[jj].replace(' ',''),list_of_ind)
            elif param=='atom_index':
                list_of_ind=g_atom_index(system,possib[jj].replace(' ',''),list_of_ind)
            else:
                print 'ERROR sel02: unknown parameter ',param

    return(list_of_ind)  #######after all -  output is  list of index ##############


#############Subfunction for make parenthesis in any way

def prolong_string(condition):
    b=' '
    for ii in range(len(condition)):
        if condition[ii]=='(' or condition[ii]==')':
           b+=' '
           b+=condition[ii]
           b+=' '
        else:
            b+=condition[ii]
    return b

### additional func for all parameters ###
def g_atom_name(system,possib,list_of_ind):    
        for ii in range(system.num_atoms):
            
            if system.atom[ii].name==possib:
                
                list_of_ind.append(ii)
                
        return list_of_ind

def g_resid_name(system,possib,list_of_ind):    
        for ii in range(system.num_atoms):

            if system.atom[ii].resid_name==possib:
               
                list_of_ind.append(ii)
             
        return list_of_ind

def g_donors(system,possib,list_of_ind):    
        for ii in range(system.num_atoms):
            if system.atom[ii].donor==True:
                list_of_ind.append(ii)
             
        return list_of_ind

def g_acceptors(system,possib,list_of_ind):    
        for ii in range(system.num_atoms):

            if system.atom[ii].acceptor==True:
               
                list_of_ind.append(ii)
             
        return list_of_ind

def g_atom_index(system,possib,list_of_ind):
    
    if ':' in possib:
        a=possib
        i = a[a.index('[')+1:a.index(':')]
        j = a[a.index(':')+1:a.index(']')]
       
        if int(j)>system.num_atoms:
            j=system.num_atoms
            print 'Total number of atoms in your system is ', system.num_atoms
        for kk in range(int(i),int(j)):
            list_of_ind.append(kk)
        
    return list_of_ind

##############appending from list of atom indexes###########################

def appending_sel(syst,list_of_ind):   
    temp_set=cl_set()

    ############Appending atoms to selection#################
    for ii in range(syst.num_atoms):
        if ii in list_of_ind:
           temp_set.atom.append(syst.atom[ii])
      

    # Building global attributes:
    temp_set.num_atoms=len(temp_set.atom)

 
    ##########Appending coordinates to selection#########################
    for jj in syst.coors:
        temp_frame=cl_coors()
        temp_frame.box=jj.box
  
        for ii in range(syst.num_atoms):         
                    if ii in list_of_ind:
                     
                        temp_frame.xyz.append(jj.xyz[ii])
        temp_set.coors.append(temp_frame)

   
         
    return temp_set














