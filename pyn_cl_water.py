from numpy import *
#import pyn_fort_water as f_water
import pyn_fort_general as f
import pyn_fort_water as f_water
from pyn_cl_set import *
import copy
import pickle as pic

#####################################################################################
##### Gaussian Network Model
#####################################################################################



def hbonds_water(definition=None,system1=None,system2=None,frame=None,optimize=False):

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
                else:
                    index_h=system1.water[index_water_h].H2.index
                list_hbonds[str(index_o)+'-'+str(index_h)]=f_water.wat.strength_o2h[ii,jj]
        
        # Free memory:
        f_water.wat.free_hbonds()

        return list_hbonds   # dict: 'index_O'-'indexH'=Skinner_parameter

class kinetic_network():
    
    def __init__(self,system=None,file_traj=None,begin=None,end=None):

        self.file=file_traj

        if system==None or file_traj==None or begin==None or end==None:
            print 'Error: input variables needed'
            print 'kinetic_network(system=None,file_traj=None,begin=None,end=None)'
            return None

        
        ####### INITIALIZE FORTRAN OBJECTS NEEDED#####
        system.last_frame=begin
        f_water.wat.switch=0          ## Optimization for hbonds=False in first frame
        f_water.wat.xarr=zeros(shape=(system.num_waters,3,3),order='Fortran')
        f_water.wat.iarr=zeros(shape=(system.num_waters,2,f_water.wat.nparts),order='Fortran')
        f_water.wat.darr=zeros(shape=(system.num_waters,2,f_water.wat.nparts),order='Fortran')

        ####### INITIALIZE NET#####
        nodes_ant=[0 for ii in range(system.num_waters)]
        nodes_post=[0 for ii in range(system.num_waters)]
        net=[]
        clave={}
        num_nodes=-1
        ####### INITIALIZE NET#####
        


        ###################################### first frame

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
                nodes_ant[jj]=clave[aa]
            except:
                num_nodes+=1
                clave[aa]=num_nodes
                nodes_ant[jj]=num_nodes
                net.append({})
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
                    nodes_post[jj]=clave[bb]
                except:
                    num_nodes+=1
                    clave[bb]=num_nodes
                    nodes_post[jj]=num_nodes
                    net.append({})
            ###### NET: FRAME NODES ########

            system.delete_coors()
            

            ###### NET: LINKS ##############
            for jj in range(system.num_waters):

                aa=nodes_ant[jj]
                bb=nodes_post[jj]

                try:
                    net[aa][bb]+=1
                except:
                    net[aa][bb]=1
            ###### NET: LINKS ##############
            
            ###### NET: UPDATE FRAME #######
            nodes_ant=copy.deepcopy(nodes_post)
            ###### NET: UPDATE FRAME #######
            
        ################################################# END

        self.net=net
        self.keys=clave
        self.keys_inv=dict((v,k) for k, v in clave.iteritems())

        return 

'''
    ## Provisional auxiliary functions
    def reformatting(self):
    
        for jj in range(system.num_waters):
            f_water.wat.xarr[jj,0,:]=system.frame[0].coors[system.water[jj].O.index,:]
            f_water.wat.xarr[jj,1,:]=system.frame[0].coors[system.water[jj].H1.index,:]
            f_water.wat.xarr[jj,2,:]=system.frame[0].coors[system.water[jj].H2.index,:]
'''

            



