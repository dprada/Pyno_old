from numpy import *
#import pyn_fort_water as f_water
import pyn_fort_general as f
import pyn_fort_water as f_water
from pyn_cl_set import *


#####################################################################################
##### Gaussian Network Model
#####################################################################################


class kinetic_network():
    
    def __init__(self,system=None,file_traj=None,begin=None,end=None):

        self.file=file_traj

        if system==None or file_traj==None or begin==None or end==None:
            print 'Error: input variables needed'
            print 'kinetic_network(system=None,file_traj=None,begin=None,end=None)'
            return None

        #coors_wat=zeros(shape=(system.num_waters,3,3),order='Fortran')
        
        system.last_frame=begin
        f_water.wat.switch=1
        f_water.wat.xarr=zeros(shape=(system.num_waters,3,3),order='Fortran')
        f_water.wat.iarr=zeros(shape=(system.num_waters,2,f_water.wat.nparts),order='Fortran')
        f_water.wat.darr=zeros(shape=(system.num_waters,2,f_water.wat.nparts),order='Fortran')

        for ii in range(begin,end+1):

            system.load_coors(file_traj)

            for jj in range(system.num_waters):
                f_water.wat.xarr[jj,0,:]=system.frame[0].coors[system.water[jj].O.index,:]
                f_water.wat.xarr[jj,1,:]=system.frame[0].coors[system.water[jj].H1.index,:]
                f_water.wat.xarr[jj,2,:]=system.frame[0].coors[system.water[jj].H2.index,:]

            mss=f_water.wat.microstates(system.num_waters,system.frame[0].box[0,0])

            system.delete_coors()

        return None


        net={}
        net['node1']={'node3':1}
        net['node1']['node5']=1


        '''
        def ff(aa,bb):
    if aa in net.keys() and bb in net[aa].keys():
        net[aa][bb]+=1
    elif aa in net.keys() and bb not in net[aa].keys():
        net[aa][bb]=1
    elif aa not in net.keys():
        net[aa]={bb:1}
        

net={}

def fff(aa,bb):

    try: 
        net[aa][bb]+=1
    except KeyError,e:
        if e[0]==aa:
            net[aa]={bb:1}
        elif e[0]==bb:
            net[aa][bb]=1
for ii in range(10000000):
    fff('node1','node1')
        '''



            



