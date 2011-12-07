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

        coors_wat=zeros(shape=(system.num_waters,3,3),order='Fortran')
        
        system.last_frame=begin
        switch=1

        for ii in range(begin,end+1):

            system.load_coors(file_traj)

            for jj in range(system.num_waters):
                coors_wat[jj,0,:]=system.frame[0].coors[system.water[jj].O.index,:]
                coors_wat[jj,1,:]=system.frame[0].coors[system.water[jj].H1.index,:]
                coors_wat[jj,2,:]=system.frame[0].coors[system.water[jj].H2.index,:]

            mss=f_water.wat.microstates(switch,coors_wat,system.frame[0].box[0,0],system.num_waters)
            print '              ',switch

#            print ii,system.frame[0].step, mss[0]


            system.delete_coors()

        return None


            



