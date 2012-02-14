from numpy import *
import pyn_fort_net as f_net
import copy

#####################################################################################
##### Networks
#####################################################################################

class cl_net():

    def __init__(self,file_net=None,file_keys=None,verbose=True):

        self.init_net()

        if file_net!=None:
            self.read_net(file_net)
        if file_keys!=None:
            self.read_keys(file_keys)

        if verbose:
            self.info()

        return

    def init_net(self):

        self.num_nodes=0
        self.links=[]
        self.num_links=0
        self.k_out_node=[]
        self.k_in_node=[]
        self.k_total=0
        self.k_max=0
        self.weight_node=[]
        self.weight_total=0
        self.keys={}
        self.keys_inv={}
        self.file_net=None
        self.file_keys=None

        self.Ts=False
        self.T_ind=[]
        self.T_start=[]
        self.T_weight=[]
        

    def info(self):

        print '# Network:'
        print '#', self.num_nodes, 'nodes'
        print '#', self.k_total, 'links'
        print '#', self.k_max, 'max. conectivity'
        print '#', self.weight_total, 'total weight'


    def read_net(self,name_file):

        self.file_net=name_file


        ff=open(name_file,'r')
        line=ff.readline()
        self.num_nodes=int(line.split()[0])

        k_max=int(line.split()[1])
        k_total=int(line.split()[2])

        self.T_ind=zeros(shape=(k_total),dtype=int,order='Fortran')
        self.T_start=zeros(shape=(self.num_nodes+1),dtype=int,order='Fortran')
        self.T_weight=zeros(shape=(k_total),dtype=int,order='Fortran')
        
        data=ff.read()
        ff.close()
        data2=[]
        for aa in data.split():
            data2.append(int(aa))


        jj=-1
        sumk=0
        for ii in range(self.num_nodes):
            jj+=1
            node_ind=data2[jj]
            k_out=data2[jj+1]
            weight=data2[jj+2]
            self.T_start[ii]=sumk
            self.links.append({})
            jj=jj+2
            for kk in range(k_out):
                jj+=1
                neigh=data2[jj]
                jj+=1
                flux=data2[jj]
                self.T_ind[sumk]=neigh
                self.T_weight[sumk]=flux
                sumk+=1
                self.links[ii][neigh-1]=flux
        self.T_start[self.num_nodes]=sumk

        for ii in range(self.num_nodes):
            self.k_out_node.append(len(self.links[ii]))
            self.weight_node.append(sum(self.links[ii].values()))


        self.k_max=max(self.k_out_node)
        self.num_links=sum(self.k_out_node)
        self.k_total=self.num_links
        self.weight_total=sum(self.weight_node)
        self.Ts=True
        
    def build_Ts(self):

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


    def read_keys(self,name_file):

        self.file_keys=name_file

        ff=open(name_file,'r')

        for ii in range(self.num_nodes):
            line=ff.readline()
            mss=line.split()[1]+' |'
            for jj in range(2,6):
                mss=mss+' '+line.split()[jj]
            mss=mss+' |' 
            for jj in range(6,9):
                mss=mss+' '+line.split()[jj]
            mss=mss+' |'
            for jj in range(9,12):
                mss=mss+' '+line.split()[jj]
            mss=mss+' |'
            for jj in range(12,15):
                mss=mss+' '+line.split()[jj]
            mss=mss+' |'
            for jj in range(15,18):
                mss=mss+' '+line.split()[jj]
            self.keys[mss]=int(line.split()[0])-1

        self.keys_inv=dict((v,k) for k, v in self.keys.iteritems())

    def symmetrize(self,new_net=False,verbose=True):

        temp=cl_net(verbose=False)
        temp.keys=copy.deepcopy(self.keys)
        temp.keys_inv=copy.deepcopy(self.keys_inv)
        temp.file_net=copy.deepcopy(self.file_net)
        temp.file_keys=copy.deepcopy(self.file_keys)
        temp.num_nodes=copy.deepcopy(self.num_nodes)

        aux={}
        for ii in range(self.num_nodes):
            for jj in self.links[ii].keys():
                aux[(ii,jj)]=0
                aux[(jj,ii)]=0
        temp.k_total=len(aux)
        del(aux)        

        if self.Ts==False :

            self.build_Ts()
            

        pfff=f_net.funcs.symmetrize_net(temp.k_total,self.T_ind,self.T_weight,self.T_start,self.num_nodes,self.k_total)
        temp.k_max=pfff[0]
        temp.T_weight=pfff[1]
        temp.T_ind=pfff[2]
        temp.T_start=pfff[3]
        temp.Ts=True
        temp.weight_node=pfff[4]
        temp.weight_total=sum(temp.weight_node)
        temp.k_out_node=zeros(temp.num_nodes)
        temp.links=[]
        for ii in range(temp.num_nodes):
            temp.links.append({})
            for jj in range(temp.T_start[ii],temp.T_start[ii+1]):
                neigh=temp.T_ind[jj]
                temp.links[ii][neigh-1]=temp.T_weight[jj]
                temp.k_out_node[ii]+=1

        temp.k_in=[]
        
        if verbose==True :
            temp.info()

        return temp



    def gradient_clusters(self,verbose=True):

        if self.Ts==False :

            self.build_Ts()

            

        self.num_clusters,self.node_belongs2=f_net.funcs.grad(self.weight_node,self.T_ind,self.T_weight,self.T_start,self.num_nodes,self.k_total)

        Clust={}
        Aux={}
        for ii in range(len(self.node_belongs2)):
            try:
                Clust[self.node_belongs2[ii]].append(ii)
            except:
                Clust[self.node_belongs2[ii]]=[]
                Clust[self.node_belongs2[ii]].append(ii)
        a=0
        repre=[]
        for ii in Clust.keys():
            Aux[a]=Clust[ii]
            repre.append(int(ii))
            a+=1
        del Clust
        self.representants=repre
        self.clust_info=Aux
        Aux1={}
        for ii in range(len(self.representants)):
            Aux1[self.representants[ii]]=ii
        for ii in range(len(self.node_belongs2)):
            self.node_belongs2[ii]=Aux1[self.node_belongs2[ii]]
        self.cluster_weight=zeros(self.num_clusters)
        for ii in range(self.num_clusters):
            for jj in self.clust_info[ii]:
                self.cluster_weight[ii]+=self.weight_node[jj]

        # Output: self.clust_info, self.representants, self.node_belongs2, self.cluster_weight, self.num_clusters
        if verbose:
            print '# Number of clusters: ',self.num_clusters


    def cfep_pfold(self,A=0,B=0,num_bins=1000):

        A=A+1
        B=B+1
        self.cfep,self.key_cfep1,self.key_cfep2=f_net.funcs.cfep_pfold(A,B,self.T_ind,self.T_weight,self.T_start,num_bins,self.num_nodes,self.k_total)



#### External Functions

def traj2net(filename=None,num_particles=0,num_frames=0,output=None):

    if output==None :
        ii=filename.rfind('.')
        if ii>0 :
            output=filename[:ii]+'.pxn'
        else:
            output=filename+'.pxn'

    if filename.endswith('.bin'):
        f_net.funcs.build_net_bin(filename,output,num_particles,num_frames)
    else:
        f_net.funcs.build_net(filename,output,num_particles,num_frames)


    print ' # New network file:', output
    return None


def merge_nets(net1=None,net2=None,verbose=True):

    net_total=cl_net(verbose=False)

    # merging the labels and keys of nodes

    net_total.keys=copy.deepcopy(net1.keys)
    net_total.keys_inv=copy.deepcopy(net1.keys_inv)
    total_num_nodes=net1.num_nodes

    net2_to_total=[0 for x in range(net2.num_nodes)]

    for ii in range(net2.num_nodes):
        try :
            jj=net1.keys[net2.keys_inv[ii]]
            net2_to_total[ii]=jj
        except:
            net_total.keys[net2.keys_inv[ii]]=total_num_nodes
            net_total.keys_inv[total_num_nodes]=net2.keys_inv[ii]
            net2_to_total[ii]=total_num_nodes
            total_num_nodes+=1

    net_total.num_nodes=total_num_nodes

    
    # merging the links

    net_total.links=copy.deepcopy(net1.links)

    for ii in range(net1.num_nodes,total_num_nodes):
        net_total.links.append({})

    for aa in range(net2.num_nodes):
        aaa=net2_to_total[aa]
        for ii,jj in net2.links[aa].iteritems():
            iii=net2_to_total[ii]
            try:
                net_total.links[aaa][iii]+=jj
            except:
                net_total.links[aaa][iii]=jj

    for ii in range(net_total.num_nodes):
       net_total.k_out_node.append(len(net_total.links[ii]))
       net_total.weight_node.append(sum(net_total.links[ii].values()))
       net_total.num_links=sum(net_total.k_out_node)
       net_total.k_total=net_total.num_links
       net_total.k_max=max(net_total.k_out_node)
       net_total.weight_total=sum(net_total.weight_node)

    net_total.build_Ts()

    if verbose:
        net_total.info()

    return net_total
