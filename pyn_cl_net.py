from numpy import *
import pyn_fort_net as f_net
import copy

#####################################################################################
##### Networks
#####################################################################################

class cl_node():

    def __init__(self):
        self.label=''
        self.link={}
        self.k_out=0
        self.k_in=0
        self.weight=0
        self.cluster=0
        self.coors=[]

    def most_weighted_links(self,length=1):
        aux_bak=[[self.link[x],x] for x in self.link.keys()]
        aux_bak.sort(reverse=True)
        most_w_destin=[]
        for ii in range(length):
            most_w_destin.append(aux_bak[ii][1])
        return most_w_destin

class cl_cluster():

    def __init__(self):
        self.label=''
        self.link={}
        self.nodes=[]
        self.num_nodes=0
        self.weight=0
        self.k_out=0
        self.k_in=0

class cl_net():

    def __init__(self,file_net=None,file_keys=None,directed=True,verbose=True):

        self.init_net(directed)

        if file_net!=None:
            self.read_net(file_net)
        if file_keys!=None:
            self.read_keys(file_keys)

        if verbose:
            self.info()

        return

    def init_net(self,directed):

        self.directed=directed
        self.num_nodes=0
        self.num_links=0
        self.num_clusters=0
        self.num_components=0
        self.node=[]
        self.cluster=[]
        self.component=[]
        self.k_total=0
        self.k_max=0
        self.k_out=0
        self.k_in=0
        self.k_out_max=0
        self.k_in_max=0
        self.weight=0
        self.labels={}
        self.clustering_method=' '
        self.directed=directed

        self.file_net=None
        self.file_keys=None

        self.Ts=False
        self.T_ind=[]
        self.T_start=[]
        self.T_weight=[]
        
    def info(self):

        print '# Network:'
        print '#', self.num_nodes, 'nodes'
        print '#', self.k_total, 'links out'
        print '#', self.weight, 'total weight nodes'


    def add_node(self, new_node, weight=0):

        node=cl_node()
        node.label=str(new_node)
        if weight!=0:
            node.weight=weight
            self.weight+=weight
        new_ind=len(self.node)
        self.node.append(node)
        self.labels[node.label]=new_ind
        self.num_nodes=new_ind+1

        return

    def add_link(self,node_origin,node_final,weight=0):

        v=[]
        for aa in [node_origin,node_final]:
            try:
                ind=self.labels[str(aa)]
            except:
                ind=self.num_nodes
                bb=self.add_node(aa)
            v.append(ind)

        if self.directed :
            try:
                self.node[v[0]].link[v[1]]+=weight
            except:
                self.node[v[0]].link[v[1]]=weight
                self.node[v[0]].k_out+=1
                self.node[v[1]].k_in+=1
                self.k_out+=1
                self.k_total+=1
        else:
            try:
                self.node[v[1]].link[v[0]]+=weight
            except:
                self.node[v[1]].link[v[0]]=weight
                self.node[v[1]].k_out+=1
                self.node[v[0]].k_in+=1
                self.k_out+=1
                self.k_total+=1
        pass

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

            node=cl_node()
            node_ind=data2[jj]
            k_out=data2[jj+1]
            weight=data2[jj+2]
            self.T_start[ii]=sumk

            jj=jj+2
            for kk in range(k_out):
                jj+=1
                neigh=data2[jj]
                jj+=1
                flux=data2[jj]
                self.T_ind[sumk]=neigh
                self.T_weight[sumk]=flux
                sumk+=1
                node.link[neigh-1]=flux

            node.k_out=len(node.link)
            node.weight=sum(node.link.values())
            self.node.append(node)

        self.T_start[self.num_nodes]=sumk

        self.k_max=k_max
        self.num_links=0
        self.weight=0
        for ii in self.node:
            self.num_links+=ii.k_out
            self.weight+=ii.weight
        self.k_total=k_total
        self.Ts=True
        
    def build_Ts(self):

        self.T_ind=zeros(shape=(self.k_total),dtype=int,order='Fortran')
        self.T_start=zeros(shape=(self.num_nodes+1),dtype=int,order='Fortran')
        self.T_weight=zeros(shape=(self.k_total),dtype=int,order='Fortran')

        kk=0
        for ii in range(self.num_nodes):
            self.T_start[ii]=kk
            aux_links=self.node[ii].link.keys()
            if ii in aux_links:
                self.T_ind[kk]=ii+1
                self.T_weight[kk]=self.node[ii].link[ii]
                aux_links.remove(ii)
                kk+=1
            for jj in aux_links:
                self.T_ind[kk]=jj+1
                self.T_weight[kk]=self.node[ii].link[jj]
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
            index=int(line.split()[0])-1
            self.labels[mss]=index
            self.node[index].label=mss
            
    def detailed_balance_distance(self,p=1.000):

        if self.Ts==False :

            self.build_Ts()
            
        db_dist=f_net.funcs.detailed_balance_distance(p,self.T_start,self.T_ind,self.T_weight,self.num_nodes,self.k_total)

        return db_dist

    def evolution_step(self,vect_in):

        if self.Ts==False:
            self.build_Ts()

        vect_out=f_net.funcs.evolution_step(self.T_start,self.T_ind,self.T_weight,vect_in,self.num_nodes,self.k_total)

        return vect_out

    def symmetrize(self,new_net=False,verbose=True):

        temp=cl_net(verbose=False)
        temp.labels=copy.deepcopy(self.labels)
        temp.file_net=copy.deepcopy(self.file_net)
        temp.file_keys=copy.deepcopy(self.file_keys)
        temp.num_nodes=copy.deepcopy(self.num_nodes)

        aux={}
        for ii in range(self.num_nodes):
            for jj in self.node[ii].link.keys():
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
        temp.weight=0
        for ii in range(temp.num_nodes):
            node=cl_node()
            node.weight=pfff[4][ii]
            for jj in range(temp.T_start[ii],temp.T_start[ii+1]):
                neigh=temp.T_ind[jj]
                node.link[neigh-1]=temp.T_weight[jj]
            node.k_out=temp.T_start[ii+1]-temp.T_start[ii]
            temp.weight+=node.weight
            node.label=self.node[ii].label
            temp.node.append(node)

        if verbose==True :
            temp.info()

        del(pfff)

        return temp


    def gradient_clusters(self,verbose=True):

        if self.Ts==False :

            self.build_Ts()


        self.num_clusters,pfff=f_net.funcs.grad(self.T_ind,self.T_weight,self.T_start,self.num_nodes,self.k_total)


        Clust={}

        for ii in range(self.num_nodes):
            try:
                Clust[pfff[ii]].append(ii)
            except:
                Clust[pfff[ii]]=[]
                Clust[pfff[ii]].append(ii)


        a=0
        for ii in Clust.keys():
            temp=cl_cluster()
            temp.label=self.node[int(ii)].label
            temp.nodes=Clust[ii]
            temp.num_nodes=len(temp.nodes)
            temp.weight=0
            for jj in temp.nodes:
                self.node[jj].cluster=a
                temp.weight+=self.node[jj].weight
            self.cluster.append(temp)
            a+=1


        # Output: self.clust_info, self.representants, self.node_belongs2, self.cluster_weight, self.num_clusters
        if verbose:
            print '# Number of clusters: ',self.num_clusters

    def clusters_links(self,verbose=True):

        if self.Ts==False:

            self.build_Ts()

        if self.num_clusters < 2:

            print '#Error: Number of clusters lower than 2'
            return

        for ii in self.node:
            c_a=ii.cluster
            for jj,ww in ii.link.items():
                c_b=self.node[jj].cluster
                try:
                    self.cluster[c_a].link[c_b]+=ww
                except:
                    self.cluster[c_a].link[c_b]=ww



    def cfep_pfold(self,A=0,B=0,num_bins=1000):

        A=A+1
        B=B+1
        self.cfep,self.key_cfep1,self.key_cfep2=f_net.funcs.cfep_pfold(A,B,self.T_ind,self.T_weight,self.T_start,num_bins,self.num_nodes,self.k_total)



#### External Functions

#def traj2net(filename=None,num_particles=0,num_frames=0,output=None):
# 
#    if output==None :
#        ii=filename.rfind('.')
#        if ii>0 :
#            output=filename[:ii]+'.pxn'
#        else:
#            output=filename+'.pxn'
# 
#    if filename.endswith('.bin'):
#        f_net.funcs.build_net_bin(filename,output,num_particles,num_frames)
#    else:
#        f_net.funcs.build_net(filename,output,num_particles,num_frames)
# 
# 
#    print ' # New network file:', output
#    return None


def merge_nets(net1=None,net2=None,verbose=True):

    net_total=cl_net(verbose=False)

    # merging the labels and keys of nodes

    net_total.labels=copy.deepcopy(net1.labels)
    total_num_nodes=net1.num_nodes

    net2_to_total=[0 for x in range(net2.num_nodes)]

    for ii in range(net1.num_nodes):
        temp=cl_node()
        temp.label=copy.deepcopy(net1.node[ii].label)
        temp.link=copy.deepcopy(net1.node[ii].link)
        net_total.node.append(temp)


    for ii in range(net2.num_nodes):
        try :
            jj=net1.labels[net2.node[ii].label]
            net2_to_total[ii]=jj
        except:
            net_total.labels[net2.node[ii].label]=total_num_nodes
            temp=cl_node()
            temp.label=net2.node[ii].label
            net_total.node.append(temp)
            net2_to_total[ii]=total_num_nodes
            total_num_nodes+=1

    net_total.num_nodes=total_num_nodes

    
    # merging the links

    for aa in range(net2.num_nodes):
        aaa=net2_to_total[aa]
        for ii,jj in net2.node[aa].link.iteritems():
            iii=net2_to_total[ii]
            try:
                net_total.node[aaa].link[iii]+=jj
            except:
                net_total.node[aaa].link[iii]=jj
    
    net_total.num_links=0
    net_total.weight=0
    for ii in net_total.node:
        ii.k_out=len(ii.link)
        ii.weight=sum(ii.link.values())
        net_total.num_links+=ii.k_out
        net_total.weight+=ii.weight
        if (net_total.k_max<ii.k_out): 
            net_total.k_max=ii.k_out

    net_total.k_total=net_total.num_links

    net_total.build_Ts()

    if verbose:
        net_total.info()

    return net_total
