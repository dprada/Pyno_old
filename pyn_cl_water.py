from numpy import *
import pyn_fort_general as f
import pyn_fort_water as f_water
from pyn_cl_set import *
from pyn_cl_net import *
import copy
import pickle as pic

#####################################################################################
##### Water tools
#####################################################################################



def hbonds_water(definition=None,system1=None,system2=None,frame=None,sk_param=0.00850,optimize=False):

    # Set up parameters as the Skinner Parameter:

    f_water.wat.sk_param=sk_param

    # Reset of previous hbonds

    for ii in system1.atom :
        ii.hbonds=[]

    for ii in system1.water :
        ii.O.hbonds=[]
        ii.H1.hbonds=[]
        ii.H2.hbonds=[]

    # Skinner

    if definition=='Skinner':

        if system2==None and optimize==False:
            if frame==None:
                frame=system1.last_frame

        f_water.wat.switch=0     # Optimization for hbonds=False 
            
        f_water.wat.xarr=zeros(shape=(system1.num_waters,3,3),order='Fortran')
        f_water.wat.iarr=zeros(shape=(system1.num_waters,2,f_water.wat.nparts),order='Fortran')
        f_water.wat.darr=zeros(shape=(system1.num_waters,2,f_water.wat.nparts),order='Fortran')

        for jj in range(system1.num_waters):
            f_water.wat.xarr[jj,0,:]=system1.frame[0].coors[system1.water[jj].O.index,:]
            f_water.wat.xarr[jj,1,:]=system1.frame[0].coors[system1.water[jj].H1.index,:]
            f_water.wat.xarr[jj,2,:]=system1.frame[0].coors[system1.water[jj].H2.index,:]

        f_water.wat.hbonds_skinner(system1.num_waters,system1.frame[0].box[0,0])

        # hbonds already in fortran variables... :
        ## f_water.wat.num_o2h[i]                 number of hbonds of the atom O of water i
        ## f_water.wat.o2h[i,dim=6]               list of water index bonded to water i because of the O
        ## f_water.wat.o2which[i,dim=6]           index of hydrogen -1 or 2- corresponding to f_water.wat.o2h
        ## f_water.wat.strength_o2h[i,dim=6]      skinner parameter corresponding to f_water.wat.o2h
        ## f_water.wat.num_h2i[i,j-1]             number of hbonds of the atom Hj of water i
        ## f_water.wat.h2o[i,j-1,dim=6]           list of water index bonded to water i because of the atom Hj
        ## f_water.wat.strength_h2o[i,j-1,dim=6]  skinner parameter corresponding to f_water.wat.h2o

        # Reformatting data:
        list_hbonds={}
        for ii in range(system1.num_waters):
            index_water_o=ii
            index_o=system1.water[ii].O.index
            for jj in range(f_water.wat.num_o2h[ii]):
                index_water_h=f_water.wat.o2h[ii,jj]-1
                if f_water.wat.o2which[ii,jj]==1:
                    index_h=system1.water[index_water_h].H1.index
                    aux_link_h=system1.water[index_water_h].H1.hbonds
                else:
                    index_h=system1.water[index_water_h].H2.index
                    aux_link_h=system1.water[index_water_h].H2.hbonds
                list_hbonds[str(index_o)+'-'+str(index_h)]=f_water.wat.strength_o2h[ii,jj]
                system1.water[index_water_o].O.hbonds.append([index_h,f_water.wat.strength_o2h[ii,jj]])
                system1.atom[index_o].hbonds.append([index_h,f_water.wat.strength_o2h[ii,jj]])
                aux_link_h.append([index_o,f_water.wat.strength_o2h[ii,jj]])
                system1.atom[index_h].hbonds.append([index_o,f_water.wat.strength_o2h[ii,jj]])

        # Free memory:
        f_water.wat.free_hbonds()

        return list_hbonds   # dict: 'index_O'-'indexH'=Skinner_parameter

def skinner_parameter(system=None,index_wat_o=None,index_wat_h=None,index_h=None,frame=None):

    if frame==None:
        frame=system.last_frame

    f_water.wat.switch=0     # Optimization for hbonds=False 
            
    f_water.wat.xarr=zeros(shape=(system.num_waters,3,3),order='Fortran')
    f_water.wat.iarr=zeros(shape=(system.num_waters,2,f_water.wat.nparts),order='Fortran')
    f_water.wat.darr=zeros(shape=(system.num_waters,2,f_water.wat.nparts),order='Fortran')

    for jj in range(system.num_waters):
        f_water.wat.xarr[jj,0,:]=system.frame[0].coors[system.water[jj].O.index,:]
        f_water.wat.xarr[jj,1,:]=system.frame[0].coors[system.water[jj].H1.index,:]
        f_water.wat.xarr[jj,2,:]=system.frame[0].coors[system.water[jj].H2.index,:]

    sk=f_water.wat.skinner_parameter(index_wat_o,index_wat_h,index_h,system.frame[0].box[0,0])

    return sk   # dict: 'index_O'-'indexH'=Skinner_parameter



def mss_water (system=None,ind_waters=False,hb_definition='Skinner',sk_param=0.00850):
    
     # Set up parameters as the Skinner Parameter:

    f_water.wat.sk_param=sk_param

    if system==None:
        print 'Error: input variables needed'
        print 'mss_water(system=None)'
        return None

        ####### INITIALIZE FORTRAN OBJECTS NEEDED#####
    f_water.wat.switch=0          ## Optimization for hbonds=False in first frame
    f_water.wat.xarr=zeros(shape=(system.num_waters,3,3),order='Fortran')
    f_water.wat.iarr=zeros(shape=(system.num_waters,2,f_water.wat.nparts),order='Fortran')
    f_water.wat.darr=zeros(shape=(system.num_waters,2,f_water.wat.nparts),order='Fortran')


    for jj in range(system.num_waters):
        f_water.wat.xarr[jj,0,:]=system.frame[0].coors[system.water[jj].O.index,:]
        f_water.wat.xarr[jj,1,:]=system.frame[0].coors[system.water[jj].H1.index,:]
        f_water.wat.xarr[jj,2,:]=system.frame[0].coors[system.water[jj].H2.index,:]

    if ind_waters==False:
        mss=f_water.wat.microstates(system.num_waters,system.frame[0].box[0,0])
    elif ind_waters==True:
        mss=f_water.wat.microstates_ind_wat(system.num_waters,system.frame[0].box[0,0])

    return mss
        


class kinetic_network(cl_net):
    
    def __init__(self,system=None,file_traj=None,begin=None,end=None,hb_definition='Skinner',sk_param=0.00850,verbose=True):

        self.init_net()
        self.file_traj=file_traj

        if system==None or file_traj==None or begin==None or end==None:
            print 'Error: input variables needed'
            print 'kinetic_network(system=None,file_traj=None,begin=None,end=None)'
            return None

        
        # Set up parameters as the Skinner Parameter:
        
        f_water.wat.sk_param=sk_param

        ####### INITIALIZE FORTRAN OBJECTS NEEDED#####
        system.last_frame=begin
        f_water.wat.switch=0          ## Optimization for hbonds=False in first frame
        f_water.wat.xarr=zeros(shape=(system.num_waters,3,3),order='Fortran')
        f_water.wat.iarr=zeros(shape=(system.num_waters,2,f_water.wat.nparts),order='Fortran')
        f_water.wat.darr=zeros(shape=(system.num_waters,2,f_water.wat.nparts),order='Fortran')

        ####### INITIALIZE NET#####
        nodes_ant=[0 for ii in range(system.num_waters)]
        nodes_post=[0 for ii in range(system.num_waters)]

        num_nodes=-1
        ####### INITIALIZE NET#####



        ###################################### first frame
        system.delete_coors()
        system.load_coors(file_traj)

        for jj in range(system.num_waters):
            f_water.wat.xarr[jj,0,:]=system.frame[0].coors[system.water[jj].O.index,:]
            f_water.wat.xarr[jj,1,:]=system.frame[0].coors[system.water[jj].H1.index,:]
            f_water.wat.xarr[jj,2,:]=system.frame[0].coors[system.water[jj].H2.index,:]

        mss=f_water.wat.microstates(system.num_waters,system.frame[0].box[0,0])
        


        ###### NET: 1ST FRAME NODES ########
        for jj in range(system.num_waters):
            aa=str(mss[jj])
            try:
                nodes_ant[jj]=self.keys[aa]
            except:
                num_nodes+=1
                self.keys[aa]=num_nodes
                nodes_ant[jj]=num_nodes
                self.links.append({})
        ###### NET: 1ST FRAME NODES ########

        system.delete_coors()
        
        ###################################### Remaining frames

        for ii in range(begin+1,end+1):

            system.load_coors(file_traj)

            for jj in range(system.num_waters):
                f_water.wat.xarr[jj,0,:]=system.frame[0].coors[system.water[jj].O.index,:]
                f_water.wat.xarr[jj,1,:]=system.frame[0].coors[system.water[jj].H1.index,:]
                f_water.wat.xarr[jj,2,:]=system.frame[0].coors[system.water[jj].H2.index,:]
            
            mss=f_water.wat.microstates(system.num_waters,system.frame[0].box[0,0])


            ###### NET: FRAME NODES ########
            for jj in range(system.num_waters):
                bb=str(mss[jj])
                try:
                    nodes_post[jj]=self.keys[bb]
                except:
                    num_nodes+=1
                    self.keys[bb]=num_nodes
                    nodes_post[jj]=num_nodes
                    self.links.append({})
            ###### NET: FRAME NODES ########

            system.delete_coors()
            

            ###### NET: LINKS ##############
            for jj in range(system.num_waters):

                aa=nodes_ant[jj]
                bb=nodes_post[jj]

                try:
                    self.links[aa][bb]+=1
                except:
                    self.links[aa][bb]=1
            ###### NET: LINKS ##############
            
            ###### NET: UPDATE FRAME #######
            nodes_ant=copy.deepcopy(nodes_post)
            ###### NET: UPDATE FRAME #######
            
        ################################################# END

        self.keys_inv=dict((v,k) for k, v in self.keys.iteritems())
        self.num_nodes=len(self.keys)
        for ii in range(self.num_nodes):
            self.k_out_node.append(len(self.links[ii]))
            self.weight_node.append(sum(self.links[ii].values()))
        self.num_links=sum(self.k_out_node)
        self.k_total=self.num_links
        self.k_max=max(self.k_out_node)
        self.weight_total=sum(self.weight_node)

        self.T_ind=zeros(shape=(self.k_total),dtype=int,order='Fortran')
        self.T_start=zeros(shape=(self.num_nodes+1),dtype=int,order='Fortran')
        self.T_weight=zeros(shape=(self.k_total),dtype=int,order='Fortran')

        kk=0
        for ii in range(self.num_nodes):
            self.T_start[ii]=kk
            aux_links=self.links[ii].keys()
            if ii in aux_links:
                self.T_ind[kk]=ii+1
                self.T_weight[kk]=self.links[ii][ii]
                aux_links.remove(ii)
                kk+=1
            for jj in aux_links:
                self.T_ind[kk]=jj+1
                self.T_weight[kk]=self.links[ii][jj]
                kk+=1
        self.T_start[self.num_nodes]=kk
        self.Ts=True

        if verbose:
            self.info()

        return 
        

'''
    ## Provisional auxiliary functions
    def reformatting(self):
    
        for jj in range(system.num_waters):
            f_water.wat.xarr[jj,0,:]=system.frame[0].coors[system.water[jj].O.index,:]
            f_water.wat.xarr[jj,1,:]=system.frame[0].coors[system.water[jj].H1.index,:]
            f_water.wat.xarr[jj,2,:]=system.frame[0].coors[system.water[jj].H2.index,:]
'''

            


