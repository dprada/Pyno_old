####################################
# GENERAL COMMENTS
#
# This module requires:
#

from os import system
from os import path
from os import sys
from pyn_cl_coors import *
sys.path.append('/home/diego/Projects/Pynoramix/backdoor_Water3/top_par/')
import top_par as tp
from numpy import *
import copy
import pylab
import pyn_fort_general as f
import pickle as pic
import datetime as datetime

#
# Structure of the file:
#
# -Class Set
#            - Instantiation
#            - Functions
# -Class Unit
# -Class Residue
# -Class Water
# -External functions
#

#
# END GENERAL COMMENTS
####################################




#######################################################
#######################################################
#### CLASSES
    
####
#### Labels or common attributes to be inherited
####

class labels_unit():                           # Every unit (atom) has at least these attributes
    def __init__(self):
        self.name=None
        self.index=None
        self.pdb_index=None
        self.covalent_bonds=[]                 # List of atoms covalently bonded.
        self.hbonds=[]                         # Atom h-bonded: [[atom_index][strength]]

class labels_set(labels_unit):                 # Every set of units (chain, residue, molecule) has at least these attributes
    def __init__(self):
        self.num_atoms=0
        self.list_atoms=[]

class labels_parent(labels_unit):               # A selection always keep memory of the previous location
    def __init__(self,parent,argument=None,unit=False):
        if unit:
            self.name=parent.name
            self.index=parent.index
            self.pdb_index=parent.pdb_index
        else:
            self.name=parent.name
            self.condition=argument

####
#### Class Unit (atom)
####


class cl_unit(labels_unit):                     # Attributes of an atom

    def __init__(self):
        '''Initialize an atom object'''

        # From labels_unit: .name, .index, .pdb_index, .covalent_bonds    

        # > Topological properties

        self.resid=labels_unit()        # Residue which this atom belongs to. 
        self.chain=labels_unit()        # Chain which this atom belongs to.
        self.parent=None                # Selection which this atom comes from.

        self.alt_loc=0                  # Alternate location (PDB)
        self.code_ins_res=0             # Code of insertion of residues (PDB)
        self.seg_ident=''               # Index segment (PDB)
        self.elem_symb=''               # Element symbol (PDB)
        self.type_pdb=''                # Type of atom for the PDB (ATOM,HETATM)

        self.covalent_bonds=[]          # esto deberia estar heredado de labels_unit @!!!!!!@

        # > Physical and Chemical properties

        self.mass=0.0                   # Mass
        self.charge=0.0                 # Charge
        self.vdw=0.0                    # VdW radius
        self.occup=0.0                  # Occupation (PDB)
        self.bfactor=0.0                # B-Factor
        self.acceptor=False             # True or false 
        self.donor=False                # True or false
        self.polar_class='Nothing'      # Acceptor or Donnor or nothing
        self.polarizability=False       # True of falsel

####
#### Class residue (set of atoms)
####


class cl_residue(labels_set):           # Attributes of a residue (inherits from label_set)

    def __init__( self ):

        # From labels_set: .name, .index, .pdb_index, .num_atoms, .list_atoms

        self.chain=labels_unit()         # Chain which this residue belongs to.

        pass


####
#### Class water (set of atoms)
####

class cl_water(labels_set):             # Attributes of a water molecule
    '''Water specific attributes and methods'''

    def __init__( self ,water_model=None):

        # From labels_set: .name, .index, .pdb_index, .num_atoms, .list_atoms

        self.model=water_model
                
        self.O=labels_unit()
        self.H1=labels_unit()
        self.H2=labels_unit()
        self.uvect_norm=[]
        #def get_uvect_norm(self):
        
        pass

####
#### Class set (set of atoms: Molecule)
####


class molecule(labels_set):               # The suptra-estructure: System (waters+cofactors+proteins...)

    
    def __init__(self,input_file=None,download=None,coors=True,verbose=True,with_bonds=False):

        # From labels_set: .name, .index, .pdb_index, .num_atoms, .list_atoms

        # > Instantation options:
        self.file=input_file            # pdb,gro...
        self.file_hbonds=''             # still not useful -do not remove-
        self.file_mss=''                # still not useful -do not remove-
        self.file_shell=''              # still not useful -do not remove-
        #self.selection=False            # is the set a selection (a mirror -pointer- of a parent set)
        
        # > Topological properties
        self.atom=[]                    # list of atoms    (objects: cl_unit)
        self.resid=[]                   # list of residues (objects: molecule)
        self.chain=[]                   # list of chains   (objects: molecule)
        self.chains=[]                  # list of chain names (strings)
        self.ion=[]                     # list of ions (objects: molecule)
        self.water=[]                   # list of waters   (objects: cl_water)
        self.water_model=None           # water model
        self.parent=None                # Parent set if it comes from selection (object: labels_parent)

        # > Physical and Chemical properties
        self.acceptors=[]               # list of acceptors
        self.donors=[]                  # list of donnors
        self.donors_hydrogen={}         # dict. (index donor atom : index H)
        self.dimensionality=0           # dimensionality (num_atoms*3)

        # > Coordinates and Trajectory
        self.file_traj=labels_file_traj()
        self.frame=[]                   # list of coordinates (objects: cl_coors)
        self.num_frames=0               # number of frames or models
        self.last_frame=0               # auxilary index of the last frame analysed
        self.dist_matrix=[]             # distance matrix (This should go to cl_coors?)

        # > Info PDB
        self.pdb_header=[]              # PDB Header (HEAD + TITLE)
        self.pdb_ss=[]                  # PDB Secondary structure
        

        ##################################

        # A SET CAN BE BUILT FROM A FILE OR FROM A SELECTION

        # IF IT COMES FROM A FILE, DOES IT NEED TO BE DOWNLOADED?

        if download:
            if not download.endswith('.pdb'):
                    download=download+'.pdb'
            input_file=download
            if not path.exists(input_file):
                temp='wget -nv http://www.rcsb.org/pdb/files/'+input_file
                system(temp)

            else:
                print '# The file '+input_file+' exists in the local folder.'
                print '# Loading it...'

            self.file=input_file

        # BUILDING THE SET FROM A FILE

        if input_file:

            # The file does not exist?
            if not download and not path.exists(input_file):
                print "# The file "+input_file+" does not exist."
                return

            # Reading the file and
            # attaching the atoms: self.atom[] (cl_unit)

            if self.file.endswith('pdb'):
                self.read_pdb(self.file)

            if self.file.endswith('gro'):
                self.read_gro(self.file)

            # Finnal set up of the attributes of the cl_units in self.atom[]:
            
            before_resid=None
            before_chain=None
            jj=-1
            ii=-1
            kk=-1
            for atom in self.atom:
                if atom.resid.pdb_index!=before_resid :
                    before_resid=atom.resid.pdb_index
                    jj+=1   
                if atom.chain.name!=before_chain :
                    before_chain=atom.chain.name
                    kk+=1
                ii+=1
                atom.index=ii                   # atom index for pynoramix
                atom.resid.index=jj             # resid index for pynoramix
                atom.chain.index=kk
                atom.hbonds=[]
            ### Setting up the subsets: residues, waters, chains.

            # Auxiliary dictionary to build resids and waters.

            aux={}

            for atom in self.atom[:]:
                ii=atom.resid.index
                try: 
                    aux[ii][atom.name]=atom.index
                except:
                    aux[ii]={}
                    aux[ii][atom.name]=atom.index

            # Residues:

            for ii in aux.keys():
                temp_residue=cl_residue()
                temp_residue.index=ii
                temp_residue.list_atoms=aux[ii].values()
                jj=temp_residue.list_atoms[0]
                temp_residue.pdb_index=self.atom[jj].resid.pdb_index
                temp_residue.name=self.atom[jj].resid.name
                self.resid.append(temp_residue)


            # Waters

            for residue in self.resid[:]:
                if tp.residue_type[residue.name]=='Water':
                    if self.water_model==None:
                        self.water_model='tip'+str(len(residue.list_atoms))+'p'
                    temp_water=cl_water()
                    temp_water.name=residue.name
                    temp_water.index=residue.index
                    temp_water.pdb_index=residue.pdb_index
                    temp_water.list_atoms=residue.list_atoms
                    temp_water.model=self.water_model

                    for ii in residue.list_atoms:
                        atom=self.atom[ii]
                        if atom.name in ['OW']:
                            atom.polar_class='acceptor'
                            atom.polarizability=True
                            temp_water.O.index=atom.index
                            temp_water.O.pdb_index=atom.pdb_index
                            temp_water.O.name='OW'
                        if atom.name in ['O']:
                            atom.polar_class='acceptor'
                            atom.polarizability=True
                            temp_water.O.index=atom.index
                            temp_water.O.pdb_index=atom.pdb_index
                            temp_water.O.name='O'
                            temp_water.H1=None
                            temp_water.H2=None
                        if atom.name in ['HW1','HW2']:
                            atom.polar_class='donor'
                            atom.polarizability=True
                            if atom.name in ['HW1']:
                                temp_water.H1.index=atom.index
                                temp_water.H1.pdb_index=atom.pdb_index
                                temp_water.H1.name='HW1'
                            elif atom.name in ['HW2']:
                                temp_water.H2.index=atom.index
                                temp_water.H2.pdb_index=atom.pdb_index
                                temp_water.H2.name='HW2'

                    self.water.append(temp_water)

            # Chains
            ii=-1
            for atom in self.atom[:]:
                if atom.type_pdb in ['ATOM']:
                    if atom.chain.name not in self.chains:
                        ii+=1
                        self.chains.append(atom.chain.name)
                        temp_chain=labels_set()
                        temp_chain.name=atom.chain.name
                        temp_chain.index=ii
                        self.chain.append(temp_chain)
                    self.chain[ii].list_atoms.append(atom.index)
            for residue in self.resid[:]:
                ii=residue.list_atoms[0]
                residue.chain.name=self.atom[ii].chain.name
                residue.chain.index=self.atom[ii].chain.index

            # Ions

            for atom in self.atom[:]:
                if tp.residue_type[atom.resid.name]=='Ion':
                    temp_residue=cl_residue()
                    temp_residue.list_atoms=[atom.index]
                    temp_residue.index=atom.resid.index
                    temp_residue.pdb_index=atom.resid.pdb_index
                    temp_residue.name=atom.resid.name
                    self.ion.append(temp_residue)

            # Deleting the auxiliary dictionary:
            del(aux)

            ### Setting up the local attributes

            # Covalent bonds

            if with_bonds:
                for residue in self.resid[:]:
                    aux_name={}
                    for ii in residue.list_atoms:
                        aux_name[tp.atom[self.atom[ii].name]]=ii
                    for ii in tp.covalent_bonds[residue.name]:
                        aa=aux_name[ii[0]]
                        bb=aux_name[ii[1]]
                        self.atom[aa].covalent_bonds.append(bb)
                        self.atom[bb].covalent_bonds.append(aa)

            # Charge

            for atom in self.atom[:]:
                if tp.atom[atom.name] in tp.charge:
                    atom.charge=tp.charge[tp.atom[atom.name]]


            # Acceptors-Donors

                # Default:
            for atom in self.atom[:]:
                if tp.atom[atom.name] in tp.donors: atom.donor=True
                if tp.atom[atom.name] in tp.acceptors: atom.acceptor=True

                # Exceptions: (This needs to be polished)
            exc_res_don=[ii[0] for ii in tp.donors_exception]
            exc_res_acc=[ii[0] for ii in tp.acceptors_exception]
            for residue in self.resid[:]:
                if residue.name in exc_res_don:
                    for exception in tp.donors_exception:
                        if residue.name == exception[0]:
                            for ii in residue.list_atoms:
                                if tp.atom[self.atom[ii].name]==exception[1]:
                                    cov=0
                                    for jj in self.atom[ii].covalent_bonds:
                                        if tp.atom_nature[self.atom[jj].name]=='H': 
                                            cov=1
                                            break
                                        if exception[2]=='Always': self.atom[ii].donor=exception[2]
                                        elif exception[2]=='Hbonded' and cov==1: self.atom[ii].donor=exception[2]
                                        elif exception[2]=='Not Hbonded' and cov==0: self.atom[ii].donor=exception[2]
                if residue.name in exc_res_acc:
                    for exception in tp.acceptors_exception:
                        if residue.name == exception[0]:
                            for ii in residue.list_atoms:
                                if tp.atom[self.atom[ii].name]==exception[1]:
                                    cov=0
                                    for jj in self.atom[ii].covalent_bonds:
                                        if tp.atom_nature[self.atom[jj].name]=='H': 
                                            cov=1
                                            break
                                        if exception[2]=='Always': self.atom[ii].acceptor=exception[2]
                                        elif exception[2]=='Hbonded' and cov==1: self.atom[ii].acceptor=exception[2]
                                        elif exception[2]=='Not Hbonded' and cov==0: self.atom[ii].acceptor=exception[2]

            for atom in self.atom[:]:
                if atom.acceptor: self.acceptors.append(atom.index)
                if atom.donor:
                    self.donors.append(atom.index)
                    for ind_cov in atom.covalent_bonds:
                        if tp.atom_nature[tp.atom[self.atom[ind_cov].name]]=='H':
                            try: 
                                self.donors_hydrogen[atom.index].append(ind_cov)
                            except:
                                self.donors_hydrogen[atom.index]=[ind_cov]

            self.acceptors=array(self.acceptors,order='Fortran')
            self.donors=array(self.donors,order='Fortran')

            ### Setting up the global attributes

            self.name=self.file[:-self.file[::-1].find('.')-1]       # file=FOO.N.pdb -> name=FOO.N
            self.num_atoms=len(self.atom)
            self.dimensionality=self.num_atoms*3
            self.num_residues=len(self.resid)
            self.num_waters=len(self.water)
            self.num_chains=len(self.chains)
            self.num_ions=len(self.ion)
            self.list_atoms=[ii for ii in range(self.num_atoms)]

            ### Loading coordinates
            if coors:
                self.load_coors(self.file)
 
            if verbose:
                self.info()

        ## END of IF input_file


    # END OF INSTANTATION

    ###
    ### INTERNAL FUNCTIONS OF A SET
    ###

    # Info function

    def info(self):

        self.num_atoms=len(self.atom)
        self.num_residues=len(self.resid)
        self.num_chains=len(self.chain)
        self.num_waters=len(self.water)
        self.num_frames=len(self.frame)
        self.num_ions=len(self.ion)
        print '#','System created from the file ',self.file,':'
        print '#',self.num_atoms,' atoms'
        print '#',self.num_residues,' residues'
        print '#',self.num_chains,' chains'
        print '#',self.num_waters,' waters'
        print '#',self.num_ions,' ions'
        print '#',self.num_frames,' frames/models'

    # To handle files

    def read_pdb (self,name_file):

        for line in open(name_file,'r'):

            ss=line.split()

            if ss[0] in ['HEADER','TITLE','CRYST1']: self.pdb_header.append(line)

            if ss[0].startswith('END'): break  # To read only the 1st model

            if ss[0] in ['HELIX','SHEET','TURN']: self.pdb_ss.append(line)

            if ss[0] in ['ATOM','HETATM']:

                temp_atom=cl_unit()
                temp_atom.type_pdb=line[0:6].replace(' ', '')
                temp_atom.pdb_index=int(line[6:11])
                temp_atom.name=(line[12:16].split())[0]
                temp_atom.alt_loc=line[16]
                temp_atom.resid.name=(line[17:20]).replace(' ', '')
                temp_atom.chain.name=line[21]
                temp_atom.resid.pdb_index=int(line[22:26])
                temp_atom.code_ins_res=line[26]
                temp_atom.occup=float(line[54:60])
                temp_atom.bfactor=float(line[60:66])
                temp_atom.seg_ident=line[72:76].replace(' ', '')
                temp_atom.elem_symb=line[76:78].replace(' ', '')
                temp_atom.charge=line[78:80].replace(' ', '')
                
                temp_atom.index=len(self.atom)
                self.atom.append(temp_atom)


    def read_gro (self,name_file):

        ## Fixed format taken from http://manual.gromacs.org/online/gro.html
        # C:   "%5d%5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f"
        # F90: "(i5,2a5,i5,3f8.3,3f8.4)"

        ff=open(name_file,'r')

        line=ff.readline()                                          # Header of the gro file

        line=ff.readline()                                        
        self.num_atoms=int(line)

        for i in range(self.num_atoms):           

            temp_atom=cl_unit()

            line=ff.readline().split()

            temp_atom.pdb_index=int(line[15:20])
            temp_atom.name=line[10:15].replace(" ", "")
            temp_atom.resid.name=line[5:10].replace(" ", "")
            temp_atom.resid.pdb_index=int(line[0:5]) 

            temp_atom.index=i           

            self.atom.append(temp_atom)


    def write_pdb (self,filename=None):
        
        if filename==None:
            print 'Enter filename: '
            print '      foo.write_pdb("foo.pdb")'
        else:
            if not filename.endswith('.pdb'): filename+='.pdb'
            if path.exists(filename): 
                print '# The file '+filename+' already exists.'
                return

            file=open(filename,'w')

            a='HEADER    '+'> CREATED BY PYNORAMIX '+datetime.datetime.now().strftime("%Y-%m-%d %H:%M")+' <\n'
            file.write(str(a))

            for ii in self.pdb_header: file.write(str(ii))
            
            for ii in self.pdb_ss:     file.write(str(ii))

            dct_aux={'ATOM': 'ATOM  ', 'HETATM': 'HETATM'}
            
            new_index=0
            for ii in range(self.num_atoms):
                new_index+=1
                a=dct_aux[self.atom[ii].type_pdb]              # 1-6
                a+="%5d" % (new_index)                         # 7-11
                #a+="%5d" % self.atom[ii].pdb_index            # 7-11
                a+=' '                                         # 12
                a+=' '+"%-3s" % self.atom[ii].name             # 13-16
                a+=' '                                         # 17
                a+="%3s" % self.atom[ii].resid.name            # 18-20
                a+=' '                                         # 21
                a+="%1s" % self.atom[ii].chain.name            # 22
                a+="%4d" % self.atom[ii].resid.pdb_index       # 23-26
                a+=' '                                         # 27
                a+='   '                                       # 28-30
                a+="%8.3f" % float(self.frame[0].coors[ii][0])   # 31-38
                a+="%8.3f" % float(self.frame[0].coors[ii][1])   # 39-46
                a+="%8.3f" % float(self.frame[0].coors[ii][2])   # 47-54
                a+="%6.2f" % self.atom[ii].occup               # 55-60
                a+="%6.2f" % self.atom[ii].bfactor             # 61-66
                a+='          '                                # 67-76
                a+="%2s" % self.atom[ii].elem_symb             # 77-78
                a+="%2s" % self.atom[ii].charge                # 79-80
                a+='\n' 
                file.write(str(a))         
                if ii<(self.num_atoms-1):
                    if self.atom[ii].type_pdb!=self.atom[ii+1].type_pdb or self.atom[ii].chain.name!=self.atom[ii+1].chain.name :
                        new_index+=1
                        a="TER   "
                        a+="%5d" % (new_index)
                        a+=' '
                        a+='  '
                        a+=' '                                         
                        a+="%3s" % self.atom[ii].resid.name            
                        a+='\n' 
                        file.write(str(a))
            a='END   '+'\n'
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


    def load_coors (self,input_file,frame=None,begin=None,end=None):  # This function needs to be polished.

        self.coors_file=input_file
        self.traj_mark=0

        if self.coors_file.endswith('pdb'):

            if frame==None:
                if begin==None:
                    if end==None:
                        frame=[]
                        for line in open(self.coors_file,'r'): 
                            if line.startswith('MODEL'): 
                                frame.append(int(line.split()[1]))
                        if len(frame)==0:
                            begin=0
                            end=0
                        else:
                            begin=frame[0]
                            end=frame[-1]
                    else:
                        for line in open(self.coors_file,'r'): 
                            if line.startswith('MODEL'): 
                                begin=frame.append(int(line.split()[1]))
                                break
                else:
                    if end==None:
                        frame=[]
                        for line in open(self.coors_file,'r'): 
                            if line.startswith('MODEL'): 
                                frame.append(int(line.split()[1]))
                        end=frame[-1]
                    
                frame=[ii for ii in range(begin,end+1)]

            if type(frame) in [int]: frame=[frame]

            for ii in frame:
                temp_frame=cl_coors(self.file,ii)
                self.frame.append(temp_frame)

        elif self.coors_file.endswith('gro'):
            temp_frame=cl_coors(self.file)
            self.frame.append(temp_frame)

        elif self.coors_file.endswith('bin'):
            if begin==None and frame==None and end==None:
                temp_frame=cl_coors(self.coors_file,self.last_frame)
                self.frame.append(temp_frame)
                self.last_frame+=1
            elif begin==None and end==None and frame!=None:
                temp_frame=cl_coors(self.coors_file,frame)
                self.frame.append(temp_frame)
                self.last_frame=frame
            elif begin!=None and end!=None:
                for ii in range(begin,end):
                    temp_frame=cl_coors(self.coors_file,ii)
                    self.frame.append(temp_frame)
                self.last_frame=end

        elif self.coors_file.endswith('xtc'):
            if begin==None and frame==None and end==None:
                if self.file_traj.status==None or self.file_traj.name!=self.coors_file:
                    self.file_traj=xdrfile(self.coors_file)
                temp_frame=cl_coors(self.coors_file,self.last_frame,self.file_traj)
                if self.file_traj.status=='END': return
                self.frame.append(temp_frame)
                self.last_frame+=1

        elif self.coors_file.endswith('trr'):
             if begin==None and frame==None and end==None:
                if self.file_traj==None:
                    self.file_traj=xdrfile(self.coors_file)
                elif self.file_traj.name!=self.coors_file:
                    self.file_traj=xdrfile(self.coors_file)
                temp_frame=cl_coors(self.coors_file,self.last_frame,self.file_traj)
                self.frame.append(temp_frame)
                self.last_frame+=1

        self.num_frames=len(self.frame)

    def delete_coors (self,frame=None,begin=None,end=None):  # This function needs to be polished
        
        if frame==begin==end==None :
            del self.frame[:]
            self.num_frames=0
            return

        if frame==None:

            if end==None:
                end=self.num_frames
            if begin==None:
                begin=1
                    
            frame=[ii for ii in range(begin,end)]

        if type(frame) in [int]: frame=[frame]
        frame.sort(); frame.reverse()

        for ii in frame:
            self.frame.__delitem__(ii)
            self.num_frames-=1


    # To handle the set

    def selection(self,condition=None):
        return make_selection(self,condition)

    def distance(self,pbc=False):
                
        for frame in self.frame:

            dist_frame=f.aux_funcs_general.dist(pbc,frame.coors,frame.box,frame.coors,self.num_atoms,self.num_atoms)
            
        self.dist_matrix.append(dist_frame)

    def neighbs(self,system2=None,limit=0,dist=0.0,pbc=False):


        if system2==None:
            system2=self
            ident=True
        else:
            ident=False

        if limit != 0:
            neighbs=f.aux_funcs_general.neighbs_limit(pbc,ident,limit,self.frame[0].coors,self.frame[0].box,system2.frame[0].coors,self.num_atoms,system2.num_atoms)
        
        else:
            print type(self.frame[0].coors[1][:]),self.frame[0].coors[1][:],system2.num_atoms
            neighbs=f.aux_funcs_general.neighbs_dist(pbc,ident,dist,self.frame[0].coors,self.frame[0].box,system2.frame[0].coors,system2.num_atoms)


        return neighbs

        
    def plot_contact_map(contact_map):
        
        pylab.gray()
        pylab.imshow(contact_map==False,origin='lower',interpolation=None) # Would be better not to interpolate
        #pylab.matshow(contact_map==False)
        return pylab.show()

    def rms_fit(self,set_reference=None,selection='all',new=False):
        
        coors_original=make_selection(self,selection).frame[0].coors
        coors_reference=make_selection(set_reference,selection).frame[0].coors

        if len(coors_original)!=len(coors_reference):
            print '# Error: Different number of atoms'
            return
        
        rot,center_ref,center_orig,rmsd,g=f.aux_funcs_general.min_rmsd(coors_reference,coors_original,len(coors_original))
        
        coors_original=self.frame[0].coors
        coors_new=f.aux_funcs_general.rot_trans(coors_original,rot,center_orig,center_ref,len(coors_original))

        if new:
            fitted_set=copy.deepcopy(self)
            fitted_set.frame[0].coors=coors_new
#            fitted_set.pdb_header="Mensaje del fitteo"
            fitted_set.rmsd=rmsd

            return fitted_set
        else:
            self.frame[0].coors=copy.deepcopy(coors_new)
            self.rmsd=rmsd

        print '# RMSD fitting:',rmsd
            # Use coors_new

        return

    def displ_vector(self,set_reference=None):

        self.d_vector=set_reference.frame[0].coors - self.frame[0].coors





#### END CLASSES
#######################################################
#######################################################



#######################################################
#######################################################
#### EXTERNAL OBJECTS AND FUNCTIONS
    
####
#### Functions
####

def min_distance(system,set_a,set_b=None,pbc=True,type_b='atoms'):

    if set_b==None:

        ind_a1,ind_a2,min_dist = f.aux_funcs_general.min_dist_atoms(pbc,True,system.frame[0].coors,system.frame[0].box,set_a,set_a,system.num_atoms,len(set_a),len(set_a))

    else:
        if type_b=='atoms':
            ind_a1,ind_a2,min_dist = f.aux_funcs_general.min_dist_atoms(pbc,False,system.frame[0].coors,system.frame[0].box,set_a,set_b,system.num_atoms,len(set_a),len(set_b))
        elif type_b=='vectors':
            l_vects=shape(set_b)
            if len(l_vects)==1:
                set_b=[set_b]
                l_vects=shape(set_b)
            array(set_b,order='Fortran')
            ind_a1,ind_a2,min_dist = f.aux_funcs_general.min_dist_atoms_ref(pbc,system.frame[0].coors,system.frame[0].box,set_a,set_b,system.num_atoms,len(set_a),l_vects[0])
    
    return ind_a1,ind_a2,min_dist


def xtc2bin(xtc_name,bin_name):
    command=home_path+'xtc2bin %s %s'%(xtc_name,bin_name)

    if path.exists(bin_name):
        command2='mv %s %s#'%(bin_name,bin_name)
        print 'file',bin_name,' was moved to ', bin_name+'#'
        system(command2)

    system(command)


def dot_product_3d(vect1,vect2):

    return f.aux_funcs_general.proj3d(vect1,vect2,len(vect1))

def isothermal_compressibility(system,Temp,input_file,frame=None,begin=None,end=None):

    frame=[ii for ii in range(begin,end+1)]

    V2a=0.0
    Va=0.0
    Kt=0.0
    for ii in frame :
        system.delete_coors()
        system.load_coors (input_file,frame=ii)
        xx=0.0
        xx=system.frame[0].box[0][0]*system.frame[0].box[1][1]*system.frame[0].box[2][2]
        Va+=xx
        V2a+=xx**2
    
    Va=Va/(len(frame)*1.0)
    V2a=V2a/(len(frame)*1.0)
    Kt=(V2a-Va**2)/(Temp*Va)
    
    return Kt*0.10,'(nm/Kb)'

#######################################################
####### Selection algorithm:    (Any suggestion to Roman)


def make_selection(system,condition):                 #####system - your system
						       #####condition - keywords to select,as example: 	OO=make_selection2(wbox,'atom_name (OW or MW) and resid.name (HO4)')

    # Some Keywords:
    if condition=='backbone':
        condition='atom_name (N or CA or CB)'
    if condition=='all':
        condition=range(system.num_atoms)

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
    					ll2.append(aux_lists[jj][ii])

		else:
			 for ii in range(len(aux_lists[jj])):
				last_list.append(aux_lists[jj][ii])

        ll2.sort()
    	sux=extracting_sel(system,ll2)


    elif type(condition) in [list,ndarray,tuple]:      ##########if you want to select from list with indexes of atoms
        LIST=condition
        sux=extracting_sel(system,LIST)


    else:
        print "ERROR sel01"
        return
    

    sux.num_frames=len(sux.frame)   ### Global variables
    sux.chains=list(set([ii.chain.name for ii in sux.atom]))
    sux.num_chains=len(sux.chains)
    #sux.parent=labels_parent(system,condition)
    #sux.selection=True
    sux.pdb_ss=system.pdb_ss

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
            elif param=='chain':
                list_of_ind=g_chain(system,possib[jj].replace(' ',''),list_of_ind)
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

            if system.atom[ii].resid.name==possib:
               
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

def g_chain(system,possib,list_of_ind):    
        for ii in range(system.num_atoms):
            
            if system.atom[ii].chain==possib:
                
                list_of_ind.append(ii)
                
        return list_of_ind

##############extracting from list of atom indexes###########################

def extracting_sel(syst,list_of_ind):   
    temp_set=molecule()

    ############Appending atoms to selection#################


    for ii in list_of_ind:
        temp_set.atom.append(syst.atom[ii])


    #########################Another attrs#####3

    # Building global attributes:

    temp_set.num_atoms=len(temp_set.atom)
    temp_set.num_waters=''
    temp_set.list_atoms=list(list_of_ind)
    temp_set.list_atoms.sort()
    temp_set.num_residues=''
    

    ##########Appending coordinates to selection#########################
    for jj in syst.frame:
        temp_frame=cl_coors()
        temp_frame.box=jj.box
  
        for ii in range(syst.num_atoms):         
                    if ii in list_of_ind:
                     
                        temp_frame.coors.append(jj.coors[ii])
        temp_set.frame.append(temp_frame)

    for ii in temp_set.frame :
        ii.coors=npy.array(ii.coors,order='Fortran')
         
    return temp_set



####### END Selection algorithm
#######################################################


#### END EXTERNAL FUNCTIONS
#######################################################
#######################################################









