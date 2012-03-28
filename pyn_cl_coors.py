from numpy import *
import struct as stc #To read binary files

class cl_coors:


    def __init__(self,file_input=None,frame=None,f_traj=None):

        self.file=file_input
        self.time=0
        self.step=0
        self.precision=0
        self.model=''
        self.coors=[]
        self.box=zeros(shape=(3,3),order='Fortran')

        if file_input!=None:

            if self.file.endswith('inp'):
                self.read_inp(self.file)
                
            if self.file.endswith('pdb'):
                self.read_pdb(self.file,frame)
            
            if self.file.endswith('gro'):
                self.read_gro(self.file)

            if self.file.endswith('bin'):
                self.read_bin(self.file,frame)

            if self.file.endswith('xtc'):
                self.read_xtc(f_traj,frame)

            if self.file.endswith('trr'):
                self.read_trr(f_traj,frame)

            self.coors=array(self.coors,order='Fortran')
        
            if self.file.endswith('bin') or self.file.endswith('gro') or self.file.endswith('xtc'):
                self.coors=10.0*self.coors
                self.box=10.0*self.box


#>>>>>>>>>>#>>>>>>>>>>
#>>>>>>>>>>#>>>>>>>>>>


    def read_pdb (self,name_file,frame):

        model=0
        for line in open(name_file,'r'):
            ii=line.split()

            if ii[0]=='CRYST1':                               # Reading the size of the box:
                self.box[0][0]=float(ii[1])                   # Using the global variable (pyn_var_glob.py)
                self.box[1][1]=float(ii[2])                   # for the size of the box vg.box
                self.box[2][2]=float(ii[3])               


            if (ii[0] in ['ATOM','HETATM']) and model==frame :

                aux=(float(line[30:38]),float(line[38:46]),float(line[46:54]))
                self.coors.append(aux)
                
            if ii[0]=='MODEL': 
                model=int(ii[1])
                if model>frame: 
                    self.model=frame
                    break

#>>>>>>>>>> FILE.GRO

    def read_gro (self,name_file):


        f=open(name_file,'r')

        line=f.readline()                                          # Header of the gro file

        line=f.readline()                         
        num_atoms=int(line)                                        # Number of atoms

        for i in range(num_atoms): 
            line=f.readline().split()
            self.coors.append(map(float,line[3:6]))

        line=f.readline().split()              # Reading the size of the box
        
        self.box[0][0]=float(line[0])          # Using the global variable (pyn_var_glob.py) 
        self.box[1][1]=float(line[1])          # for the size of the box vg.box      
        self.box[2][2]=float(line[2])       

        
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    READING THE TRAJECTORY   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



    def update_coors ( self, name_file,frame=None):     
        '''Updating coordinates'''

        self.traj_file=name_file                     # Name of the trajectory file
                                                     # self.frame was defined at the beginning
                                                     # If frame==None, the next frame will be read

        if self.traj_file.endswith('bin'):            # Reading binary traj.
            self.read_bin(name_file,frame)
            
        #if self.name_aux.endswith('xtc'):            # Readint binary traj.
        #    self.read_traj_bin(name_file,frame)      <--- To be inplemented


#>>>>>>>>>> TRAJ.BIN (Binary)
    def read_bin (self, name_file,frame):
        self.traj_file=name_file  
        FF=file(self.traj_file)

      
       # if FF.closed==False:                     # Checking if the file was not opened before
        self.f_traj=open(name_file,'rb')      
           
           
            
        N_A=stc.unpack('i', self.f_traj.read(4))[0]
	self.box[0][:]=stc.unpack('3f', self.f_traj.read(3*4))[0:3]    # Using the global variable (pyn_var_glob.py)
        self.box[1][:]=stc.unpack('3f', self.f_traj.read(3*4))[0:3]    # for the size of the box vg.box        
        self.box[2][:]=stc.unpack('3f', self.f_traj.read(3*4))[0:3]         
        
        if frame!= None: 
            self.frame=frame            # Going to the choosen frame
            bytes_frame=(13+N_A*3)*4       # Bytes per frame
            self.f_traj.seek(bytes_frame*self.frame,1)     # Jumping to the choosen frame
            
        else:
            self.frame+=1
     
        # print self.frame
        # temp_coors.num_atoms=stc.unpack('i', self.f_traj.read(4))[0]            # Reading Attributes
        N_A=stc.unpack('i', self.f_traj.read(4))[0]                              # N_A  = number of atoms

        self.step=stc.unpack('i', self.f_traj.read(4))[0]
        self.time=stc.unpack('f', self.f_traj.read(4))[0]
        self.box[0][:]=stc.unpack('3f', self.f_traj.read(3*4))[0:3]    # Using the global variable (pyn_var_glob.py)
        self.box[1][:]=stc.unpack('3f', self.f_traj.read(3*4))[0:3]    # for the size of the box vg.box        
        self.box[2][:]=stc.unpack('3f', self.f_traj.read(3*4))[0:3]                                                 
        
        format=str(3*N_A)+'f'                                   # Format of system's coordinates
        temp=stc.unpack(format, self.f_traj.read(N_A*4*3))[0:N_A*3]
        
        for ii in range(0,3*N_A,3):
            aux=temp[ii:ii+3]
            
            self.coors.append(aux)
        self.precision=stc.unpack('f',self.f_traj.read(4))[0]             # Precision of the trajectory
       
        self.f_traj.close()
       # return temp_coors


#>>>>>>>>>> TRAJ.XTC
    def read_xtc (self, f_traj,frame):

        f=f_traj.upload_frame()

        self.coors=f.x
        self.box=f.box
        self.step=f.step
        self.time=f.time
        self.precision=f.prec

#>>>>>>>>>> TRAJ.TRR
    def read_trr (self, f_traj,frame):

        f=f_traj.upload_frame()

        self.coors=f.x
        self.box=f.box
        self.step=f.step
        self.time=f.time
        self.precision=f.lam



      
###################################################
######## XDRFILES (XTC and TRR)
###################################################
        
from ctypes import *
import os.path
from numpy.ctypeslib import ndpointer

class frame:
    #variables
    #x: rvec*natoms / numpy array if installed
    #box DIM*DIM
    #step 
    #time 
    #prec 
    #lam: lambda

    def __init__(self,n):
        #create vector for x
        self.x=empty((n,3),dtype=float32)
        self.box = empty((3,3),float32)


class xdrfile:
    exdrOK, exdrHEADER, exdrSTRING, exdrDOUBLE, exdrINT, exdrFLOAT, exdrUINT, exdr3DX, exdrCLOSE, exdrMAGIC, exdrNOMEM, exdrENDOFFILE, exdrNR = range(13)

    #
    def __init__(self,fn,ft="Auto"):

        self.name=fn

        if ft=="Auto":
          ft = os.path.splitext(fn)[1][1:]
        self.mode=ft

        #load libxdrfil
        try: 
          self.xdr=cdll.LoadLibrary("libxdrfile.so")
        except:
          raise IOError("libxdrfile.so can't be loaded")
          
        #open file
        self.xd = self.xdr.xdrfile_open(fn,"r")
        if not self.xd: raise IOError("Cannot open file: '%s'"%fn)
        
        #read natoms
        natoms=c_int()
        if self.mode == 'trr':
            r=self.xdr.read_trr_natoms(fn,byref(natoms))
        else:
            r=self.xdr.read_xtc_natoms(fn,byref(natoms))
        if r!=self.exdrOK: raise IOError("Error reading: '%s'"%fn)
        self.natoms=natoms.value
        
        #for NumPy define argtypes - ndpointer is not automatically converted to POINTER(c_float)
        #alternative of ctypes.data_as(POINTER(c_float)) requires two version for numpy and c_float array
        self.xdr.read_xtc.argtypes=[c_int,c_int,POINTER(c_int),POINTER(c_float),
           ndpointer(ndim=2,dtype=float32),ndpointer(ndim=2,dtype=float32),POINTER(c_float)]
        self.xdr.read_trr.argtypes=[c_int,c_int,POINTER(c_int),POINTER(c_float),POINTER(c_float),
           ndpointer(ndim=2,dtype=float32),ndpointer(ndim=2,dtype=float32),POINTER(c_float),POINTER(c_float)]

    def __iter__(self):
        f = frame(self.natoms)
        #temporary c_type variables (frame variables are python type)
        step = c_int()
        time = c_float()
        prec = c_float()
        lam = c_float()
        while 1:
            #read next frame
            if self.mode=='xtc':
                result = self.xdr.read_xtc(self.xd,self.natoms,byref(step),byref(time),f.box,
                        f.x,byref(prec))
                f.prec=prec.value
            else:
                result = self.xdr.read_trr(self.xd,self.natoms,byref(step),byref(time),byref(lam),
                        f.box,f.x,None,None) #TODO: make v,f possible
                f.lam=lam.value
                
            #check return value
            if result==self.exdrENDOFFILE: break
            if result==self.exdrINT and self.mode=='trr': 
              break  #TODO: dirty hack. read_trr return exdrINT not exdrENDOFFILE
            if result!=self.exdrOK: raise IOError("Error reading xdr file")
            
            #convert c_type to python 
            f.step=step.value
            f.time=time.value
            yield f

    def upload_frame(self):
        f = frame(self.natoms)
        step = c_int()
        time = c_float()
        prec = c_float()
        lam = c_float()
        if self.mode=='xtc':
            result = self.xdr.read_xtc(self.xd,self.natoms,byref(step),byref(time),f.box,
                     f.x,byref(prec))
            f.prec=prec.value
        else:
            result = self.xdr.read_trr(self.xd,self.natoms,byref(step),byref(time),byref(lam),
                     f.box,f.x,None,None) #TODO: make v,f possible
            f.lam=lam.value
        #check return value
        if result==self.exdrENDOFFILE: 
            print 'End of file'
            return
        if result==self.exdrINT and self.mode=='trr': 
            print 'mmm...'
            return
            #TODO: dirty hack. read_trr return exdrINT not exdrENDOFFILE
        if result!=self.exdrOK: raise IOError("Error reading xdr file")
            
        #convert c_type to python 
        f.step=step.value
        f.time=time.value
        return f
