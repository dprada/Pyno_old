
Tutorials and Examples
*****************************

First of all, lets load Pynoramix in our script or in a ipython session:

.. sourcecode:: ipython

     In [1]: from pynoramix_beta import *


Some basic notions on python will be assumed along this tutorial. If you just landed here without any idea on python, have a look to the section *First steps on python*.

.. todo:: Make a short tutorial on python, enough to run pynoramix.

----------------------
 

Networks
===================

How to create, load and handle a network.

A network from scratch
++++++++++++++++++++++

Lets create as simple example a network of cities:

.. sourcecode:: ipython

   In [2]: cities=network()
   # Network:
   # 0 nodes
   # 0 links out
   # 0 total weight nodes

Nodes can be added in two ways, along or inferred by the links addition:

.. sourcecode:: ipython

   In [3]: cities.add_node('Zaragoza')
   In [4]: cities.add_link('Rome','Brescia',444)
   In [5]: cities.add_link('Rome','Zaragoza',1112)
   In [6]: cities.add_link('Zaragoza','Kiev',2561)
   In [7]: cities.info()
   # Network:
   # 4 nodes
   # 3 links out
   # 0 total weight nodes

The nodes are attached to the network in order of creation:

.. sourcecode:: ipython

   In [7]: cities.node[0].label
   Out[7]: 'Zaragoza'
   In [8]: cities.node[3].label
   Out[8]: 'Kiev'

Links are stored as a dictionary for each node (see ref). Nodes are from now on referred because of their indexes.

.. sourcecode:: ipython

   In [9]: cities.labels['Rome']
   Out[9]: 1
   In [10]: cities.node[1].link.keys()
   Out[10]: [0, 2]
   In [11]: print 'From Rome to Zaragoza: ', cities.node[1].link[0], 'km.'
   From Rome to Zaragoza:  1112 km.

.. seealso:: include here a link to the class definition and attributes.


Loading networks from files
+++++++++++++++++++++++++++++

A network can be loaded together with their labels. Pynoramix uses its
own compact format for the network, while the labels can be readed with many formats.
This way a network can be initialized with the files or a posteriori:

Loading a simple network from a columns file such as
:download:`net_ex1.inp <../tutorials/nets/net_ex1.inp>`,
:download:`net_ex2.inp <../tutorials/nets/net_ex2.inp>` or
:download:`net_ex3.inp <../tutorials/nets/net_ex3.inp>` can be
done as follows:

.. sourcecode:: ipython

   In [2]: net1=network('net_ex1.inp')
   # Network:
   # 5 nodes
   # 7 links out
   # 0 total weight nodes

   In [4]: net2=network(verbose=False)
   In [5]: net2.read_net('net_ex2.inp')
   # Network:
   # 5 nodes
   # 5 links out
   # 6.0 total weight nodes

And using an extra file for the labels, :download:`labels_ex3.inp
<../tutorials/nets/labels_ex3.inp>`, we can also:

.. sourcecode:: ipython

   In [6]: net3=network('net_ex3.inp','labels_ex3.inp')
   # Network:
   # 4 nodes
   # 5 links out
   # 24.0 total weight nodes
    
   In [7]: net3.labels()
   Out[7]: {'Alexandra': 1, 'Bob': 2, 'Liliana': 0, 'Tom': 3}

.. note:: Describe the parameters needed in the text input files, and make a call to the subroutines of this part (like read_labels()).

Merging networks
++++++++++++++++

Two networks can be merged into one. The function updates one of the
networks appending the new nodes and links and adding up the value of
the weights of overlapping nodes and links.

.. sourcecode:: ipython

   In [8]: net_12=network()
   # Network:
   # 0 nodes
   # 0 links out
   # 0 total weight nodes
    
   In [9]: net_12.merge_net(net1,verbose=False)
   In [10]: net_12.merge_net(net2,verbose=False)
   In [11]: net_12.info()
   # Network:
   # 6 nodes
   # 11 links out
   # 6.0 total weight nodes

----------------------

Kinetic Networks
===================

This section is a tutorial on how to analyze kinetic networks. To
illustrate the analysis some test networks are available.

Examples
+++++++++

1D double well
..............

A kinetic network has been obtained for particle in a 1D potential: 

.. math::

   x^4-4x^2+x+sin(10x) 

The files for this network are available as :download:`2w_1D.net
<../tutorials/nets/1D_2well/2w_1D.net>` and :download:`2w_1D.aux
<../tutorials/nets/1D_2well/2w_1D.aux>`. Where the topology file is in
the native format and the labels in text format.


.. sourcecode:: ipython

   In [2]: net_1D=network('2w_1D.net','2w_1D.aux',net_format='native')
   # Network:
   # 970 nodes
   # 195638 links out
   # 1000090000 total weight nodes

Since the label of each node corresponds to the bin of coordinate x, a
single value on the midle of the bin can be given to each node as
coordinate for representations.

.. sourcecode:: ipython

   In [3]: for nn in net_1D.node:
     ....:     aa=nn.label[1:-1].split(',')
     ....:     nn.coors=(float(aa[0])+float(aa[1]))/2.0
     ....: 
    
   In [4]: print net_1D.node[0].label, net_1D.node[0].coors
   [-1.530,-1.525] -1.5275

This way we can plot the stationary probability distribution of the
particle along x:

.. sourcecode:: ipython

   In [5]: xx=[]; yy=[]; delta_x=0.025
    
   In [6]: for nn in net_1D.node:
     ....:         xx.append(nn.coors); yy.append(nn.weight/(net_1D.weight*delta_x))
     ....: 
    
   In [7]: pyl.plot(xx,yy,'bo')
   Out[7]: [<matplotlib.lines.Line2D object at 0x515a110>]
    
   In [8]: pylab.show()

.. plot:: pyplots/2w_1D_fig1.py

cFEPs
++++++

How to build a cFEP from a kinetic network.

Dijkstra
++++++++

MFPT
++++

MCL
+++

Components
++++++++++

Gigant Component
................

Extracting subnetwork
+++++++++++++++++++++

Weight-core
+++++++++++

K-core
++++++


----------------------

Water
=====

How to analize the Conformational Space Network of bulk water. Add references.

The second solvation shell
++++++++++++++++++++++++++

The system is loaded from a pdb or gro file.

.. sourcecode:: ipython

   In [2]: watbox=molecule('tip4p-2005.pdb')
   # System created from the file  tip4p-2005.pdb :
   # 4096  atoms
   # 1024  residues
   # 1  chains
   # 1024  waters
   # 0  ions
   # 1  frames/models

We can already calculate the microstates for the coordinates stored from the pdb:

.. sourcecode:: ipython

   In [3]: mss_water(watbox,definition='Skinner')
   # Water microstates updated
   In [4]: watbox.water[500].microstate
   Out[4]: '1 | 2 3 4 5 | 6 7 8 | 9 10 11 | 12 13 14 | 15 0 17'

The former function can return the microstates of the system as an
array or the indexes of the water molecules behind it.

.. sourcecode:: ipython

   In [5]: mss_frame=mss_water(watbox,definition='Skinner',output_array='microstates',verbose=False)
   In [6]: ind_frame=mss_water(watbox,definition='Skinner',output_array='indexes_waters',verbose=False)
   In [7]: print mss_frame[500]; print ind_frame[500]
   [ 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15  0 17]
   [500 323 670 973 151 566  10 722 942  71 306 777 212  97 865  -1 573]


Notice that a '0' in any position of the microstate corresponds to a
'-1' in the array of water indexes. This is due to the fact that zero
is the first index of water.


The Kinetic network
+++++++++++++++++++

The system is loaded as it was described in the previous section but
only the topology will be used, not the coordinates of the initial
frame. This way these data can be removed:

.. sourcecode:: ipython

   In [2]: watbox=molecule('tip4p-2005.pdb',verbose=False)
   In [3]: watbox.delete_coors()


We can now build the kinetic network reading the frames of a trajectory:

.. sourcecode:: ipython

   In [6]: watnet=kinetic_network(watbox,'md_test.xtc',begin=0,end=100,definition='Skinner')
   # Network:
   # 2597 nodes
   # 11277 links out
   # 102400.0 total weight nodes

A kinetic network 'watnet' has been created analysing the first 100 frames of the trajectory.





