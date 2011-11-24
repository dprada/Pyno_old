##################################

    def hbonds(self,syst=None,opt=None,neighbs=None):

        if syst==None:
            syst=self
            ident=True
        else:
            ident=False

        if opt=='Skinner':
            hbonds_out=f.aux_funcs.hbonds_skinner(self.coors[0].xyz,self.coors[0].box,self.donors,self.acceptors,self.num_atoms,len(self.donors),len(self.acceptors))
            
        if opt=='Roh':
            hbonds_out=f.aux_funcs.hbonds_roh(ident,self.coors[0].xyz,self.coors[0].box,self.donors,self.acceptors,
                                          syst.coors[0].xyz,syst.donors,syst.acceptors,neighbs[0],neighbs[1],neighbs[2],
                                          self.num_atoms,len(self.donors),len(self.acceptors),syst.num_atoms,
                                          len(syst.donors),len(syst.acceptors),len(neighbs[0][0]))
        


        hbonds={}

        if ident==True :
            for ii in self.donors:
                for jj in syst.acceptors:
                    if hbonds_out[0][ii,jj]==-1 and jj not in self.atom[ii].covalent_bond[:]:
                        
                        try: 
                            hbonds[ii][jj]=hbonds_out[1][ii,jj]
                        except:
                            hbonds[ii]={}
                            hbonds[ii][jj]=hbonds_out[1][ii,jj]

                        try: 
                            hbonds[jj][ii]=hbonds_out[1][jj,ii]
                        except:
                            hbonds[jj]={}
                            hbonds[jj][ii]=hbonds_out[1][jj,ii]  



        num_hbonds=0
        for aa in hbonds.values():
            num_hbonds+=len(aa)
        print num_hbonds

        return hbonds
 ###########################################

#######################################################
#######################################################
### Analysis of the entire water trajectory:

    def water_analysis (self,traj_name=None,write_out=None,init_frame=0,last_frame=-1,RAM=0):

        if traj_name == None:
            print 'input the traj_name'
            return None
        if traj_name.endswith('.xtc'):
            if not path.exists(traj_name[:-3]+'bin'):
                xtc2bin(traj_name,traj_name[:-3]+'bin')
            traj_name=traj_name[:-3]+'bin'

        if write_out==None:
            write_out='No'
            
        if RAM==0:

            command=self.path+'./pyn_anw '+traj_name+' '+write_out+' '+self.water_model+' '+str(self.num_residues)+' '+str(init_frame)+' '+str(last_frame)
            system(command)

            self.file_hbonds='aux_hbs.bin'
            self.file_mss='aux_mss.bin'
            self.file_shell='aux_shell.bin'
            self.file_net='aux_net.oup'
            self.file_key_mss='aux_key_mss.oup'

        else:

            command=self.path+'./pyn_anw_long '+traj_name+' '+write_out+' '+self.water_model+' '+str(self.num_residues)+' '+str(init_frame)+' '+str(last_frame)
            system(command)
            system('rm trj_mss.bin')
            self.file_net='aux_net.oup'
            self.file_key_mss='aux_key_mss.oup'
         


    def get_hbonds(self,frame=None,molecule=None):

        L=[]
        HB=open(self.file_hbonds,'rb')
        HB.seek((2*4*4*(frame-1)*self.num_waters)+(4*(frame)))
        for ii in range(self.num_waters*2):    
            B=(stc.unpack('3i',HB.read(12)))
            S=(stc.unpack('1f',HB.read(4)))
            L.append([B,S])
        if molecule!=None:
            
            print L[(molecule-1)*2][0][0],'H1', L[(molecule-1)*2][0][2],'O',L[(molecule-1)*2][1][0]
            print L[(molecule-1)*2+1][0][0],'H2', L[(molecule-1)*2+1][0][2],'O',L[(molecule-1)*2+1][1][0]

            
            for jj in range(self.num_waters*2):
                if molecule==L[jj][0][2]:
                    if L[jj][0][1]==1:
                        print L[jj][0][0],'O', L[jj][0][2],'H1',L[jj][1][0]
                    else:
                        print L[jj][0][0],'O', L[jj][0][2],'H2',L[jj][1][0]
        else:
            print L


    def get_mss(self,frame=None,molecule=None):
        L=[]    
       
        BB=open(self.file_mss,'rb')

        BB.seek((17*4*(frame-1)*self.num_waters)+(4*(frame)))
        for ii in range(self.num_waters):    
            B=(stc.unpack('17i',BB.read(68)))
            L.append(B)
        if molecule!=None:
            print L[molecule-1]
        else:
            print L

    def get_shell(self,frame=None,molecule=None):
        L=[]    
       
        BB=open(self.file_shell,'rb')

        BB.seek((17*4*(frame-1)*self.num_waters)+(4*(frame)))
        for ii in range(self.num_waters):    
            B=(stc.unpack('17i',BB.read(68)))
            L.append(B)
        if molecule!=None:
            print L[molecule-1]
        else:
            print L

    def get_network(self):

        return cl_net(self.file_net,self.file_key_mss)
