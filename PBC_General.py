import numpy as np

# ABAQUS
from part import *
from material import *
from section import *
from assembly import *
from step import *                                               # Generalized code for applying Periodic Boundary conditions
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *

## Provide the model and instance name from ABAQUS

modelname = 'Model-1'
instancename = 'Part-1-1'



## Define the geometrical properties, as required
l, le = 5.0e-3, 1.0e-3
wf = 0.6e-3                   # 0.6, 0.4, 0.2, 0.8, 1.0
t = 0.163e-3                  # 0.163, 0.147, 0.128, 0.181, 0.2
x = 0.5*(le-wf)*(1+t/l)

LRVE = l + t
BRVE = 2*(wf+x+t)

## Accesing the model and instance

mod = mdb.models[modelname]
a = mod.rootAssembly
p = a.instances[instancename]

## U - Up/Top, D - Down/Bottom, L - Left, R - Right

## For applying PBC, the nodes at all the edges must be obtained. If each edge (U,D,L,R) contains several partitions, all those partitions 
## are found and combined to form a set.

# Left and right edges - 4 partitions
L1 = p.edges.findAt((0.0, round(wf/4, 6), 0.0))
L2 = p.edges.findAt((0.0, round(BRVE/2, 6), 0.0))
L3 = p.edges.findAt((0.0, round(BRVE-wf/4, 6), 0.0))
# L4 = p.edges.findAt((0.0, round(0.1, 6), 0.0))

R1 = p.edges.findAt((round(LRVE, 6), round(wf/4, 6), 0.0))
R3 = p.edges.findAt((round(LRVE, 6), round(BRVE/2, 6), 0.0))
R2 = p.edges.findAt((round(LRVE, 6), round(BRVE-wf/4, 6), 0.0))
# R4 = p.edges.findAt((round(LRVE, 6), round(0.1, 6), 0.0))


# Top and bottom - 3 and no partitions respectively
U1 = p.edges.findAt((round(l/4, 6), round(BRVE, 6), 0.0))
U2 = p.edges.findAt((round(LRVE/2, 6), round(BRVE, 6), 0.0))
U3 = p.edges.findAt((round(LRVE-l/4, 6), round(BRVE, 6), 0.0))

D1 = p.edges.findAt((round(l/4, 6), 0.0, 0.0))
D2 = p.edges.findAt((round(LRVE/2, 6), 0.0, 0.0))
D3 = p.edges.findAt((round(LRVE-l/4, 6), 0.0, 0.0))


## For each partition, the index of the particular partition is selected and stored. Index are used to easily access the edges, considered as unique identifier. 
qU1 = U1.index
qU2 = U2.index
qU3 = U3.index

qD1 = D1.index
q22 = D2.index
q23 = D3.index

qR1 = R1.index
qR2 = R2.index
qR3 = R3.index
# qR4 = R4.index

qL1 = L1.index
qL2 = L2.index
qL3 = L3.index
# qL4 = L4.index


# Here, the edges are found using the indices and stored. This is done to facilitate the combination of all the edges to form single edge.
# qU1:qU1+1 means selecting only one edge of the particular index
EdUe1 = p.edges[qU1:qU1+1]
EdUe2 = p.edges[qU2:qU2+1]
EdUe3 = p.edges[qU3:qU3+1]  

EdDe1 = p.edges[qD1:qD1+1]
EdDe2 = p.edges[q22:q22+1]
EdDe3 = p.edges[q23:q23+1]

EdRe1 = p.edges[qR1:qR1+1]
EdRe2 = p.edges[qR2:qR2+1]
EdRe3 = p.edges[qR3:qR3+1]
# EdRe4 = p.edges[qR4:qR4+1]

EdLe1 = p.edges[qL1:qL1+1]
EdLe2 = p.edges[qL2:qL2+1]
EdLe3 = p.edges[qL3:qL3+1]
# EdLe4 = p.edges[qL4:qL4+1]


# Combining all the partitions, forming a single edge and collecting the nodes. This can be done with the help of above indices code.
combinedEdgesL = EdLe1 + EdLe2 + EdLe3 
a.Set(edges=combinedEdgesL, name='Combined_LeftEdges')
Leftnodes = a.sets['Combined_LeftEdges'].nodes

combinedEdgesR = EdRe1 + EdRe2 + EdRe3 
a.Set(edges=combinedEdgesR, name='Combined_RightEdges')
Rightnodes = a.sets['Combined_RightEdges'].nodes

combinedEdgesU = EdUe1 + EdUe2 + EdUe3 
a.Set(edges=combinedEdgesU, name='Combined_UpEdges')
Upnodes = a.sets['Combined_UpEdges'].nodes

combinedEdgesD = EdDe1 + EdDe2 + EdDe3
a.Set(edges=combinedEdgesD, name='Combined_DownEdges')
Downnodes = a.sets['Combined_DownEdges'].nodes

## Creating null arrays for putting all nodal coordinates and labels
Upcoord, Downcoord, Leftcoord, Rightcoord = [], [], [], []


## Putting all the nodal information into the empty array. Final array consists of 3 columns
for node in Upnodes:
    Upcoord = Upcoord + [[node.coordinates[0], node.coordinates[1], node.label]]    # node.coordinates[0] - x coordinate of a ptcular. node
                                                                                    # node.coordinates[1] - y coordinate of a ptcular. node
                                                                                    # node.label - label of a ptcular. node
for node in Downnodes:
    Downcoord = Downcoord + [[node.coordinates[0], node.coordinates[1], node.label]]

for node in Leftnodes:
    Leftcoord = Leftcoord + [[node.coordinates[0], node.coordinates[1], node.label]]
for node in Rightnodes:
    Rightcoord = Rightcoord + [[node.coordinates[0], node.coordinates[1], node.label]]

for i in range(0,len(Leftcoord)):
    if (Leftcoord[i][0]<10E-6):           ## To prevent the problem of coordinates collection in wavy cases
        Leftcoord[i][0] = 0
    else:
        continue


## Sorting of the nodes based on the labels
Upcoord.sort()
Downcoord.sort()
Leftcoord.sort()
Rightcoord.sort()


NodeTol = 0.05/20   ## A tolerance value which will be used in selection of nodal pairs from opposite edges while applying the constraints

## Creating nodal pairs from opposite edges
Num = len(Upcoord)
for i in range(0, Num):
    if (abs(Upcoord[i][0]-Downcoord[i][0])<NodeTol):                       # Checking the x - coord. of i-th node of Up and down edge, if the values are almost same
        Nlabel = Upcoord[i][2]                                             # collecting the node label
        a.Set(nodes=p.nodes[Nlabel-1:Nlabel], name='UpNode-'+str(i))       # Selecting that ptcular node and creating a set for up nodes 
        Nlabel = Downcoord[i][2]
        a.Set(nodes=p.nodes[Nlabel-1:Nlabel], name='DownNode-'+str(i))     # for down nodes


Num = len(Leftcoord)
for i in range(0, Num):
    if (abs(Leftcoord[i][1]-Rightcoord[i][1])<NodeTol):                          # Similar process for nodes from left and right edges
        Nlabel = Rightcoord[i][2]
        a.Set(nodes=p.nodes[Nlabel-1:Nlabel], name='RightNode-'+str(i))
        Nlabel = Leftcoord[i][2]
        a.Set(nodes=p.nodes[Nlabel-1:Nlabel], name='LeftNode-'+str(i))

## UpNode-0 => Top left corner node (Y-roller - x-disp =0)
## DownNode-0, LeftNode-0 => Bottom Left corner node (Fix)
## RightNode-0 => Bottom right corner node  (X-roller - y-disp = 0)

### Applying constraints
        
        ## Constraint: 
        # Along X-dirn: UR_i - UL_i = UR_0 - UL_0   (Here, as Leftnode-0 is fixed, UL_0 is 0. So only 3 terms are present)
        #               UD_i - Uu_i = UD_0 - Uu_0   (Here, as DownNode-0 is fixed and UpNode-0 has roller support which cannot move in x-dirn, both displacements are equal to zero and only 2 terms are present in the equation)
        
        # Along Y-dirn: VU_i - VD_i = VU_0 - VD_0 (DownNode-0 - fixed => Only 3 terms)
        #               VL_i - VR_i = VL_0 - VR_0 (LeftNode-0 fixed, RightNode-0 roller support - cannot move in y-dirn => Only 2 terms in equation)  

        # The coefficients above are given according to the code below. Coefficients can be changed accordingly also, as long as the above constraints hold true.

 ## Left and Right edges       
for i in range(1, len(Rightcoord)-1):                    # Loop starts from 1 to facilitate the node numbering, and doesnot include the zeroth node (Corners)
    mod.Equation(name='Eqn-LR-X-'+str(i), terms=((-1.0, 'LeftNode-'+str(i), 1), (1.0, 'RightNode-'+str(i), 1), (-1.0, 'RightNode-0', 1)))      #X-dirn
for i in range(1, len(Rightcoord)-1):
    mod.Equation(name='Eqn-LR-Y-'+str(i), terms=((1.0, 'LeftNode-'+str(i), 2), (-1.0, 'RightNode-'+str(i), 2)))                # Y-dirn

for i in range(1, len(Upcoord)-1):
    mod.Equation(name='Eqn-UD-Y-'+str(i), terms=((-1.0, 'DownNode-'+str(i), 2), (1.0, 'UpNode-'+str(i), 2), (-1.0, 'UpNode-0', 2)))         # Y-dirn
for i in range(1, len(Upcoord)-1):
    mod.Equation(name='Eqn-UD-X-'+str(i), terms=((1.0, 'DownNode-'+str(i), 1), (-1.0, 'UpNode-'+str(i), 1)))                                # X_dirn


## Applying constraints to top right corner in both directions w.r.t the bottom right corner and top leftcorner as the corner nodes were neglected in the above code
mod.Equation(name='Eqn-TR-X', terms=((1.0, 'UpNode-'+str(len(Upcoord)-1), 1), (-1.0, 'DownNode-'+str(len(Upcoord)-1), 1)))  # X_dirn
mod.Equation(name='Eqn-TR-Y', terms=((1.0, 'UpNode-'+str(len(Upcoord)-1), 2), (-1.0, 'UpNode-0', 2)))   # Y_dirn



### Applying BCs

# Fixing Bottom left corner
v = a.instances[instancename].vertices
BL_ver = v.findAt((0.0, 0.0, 0.0))                        # Selecting the node/vertex
qbl = BL_ver.index                                        # Finding the index of the node
FixVer1 = v[qbl:qbl+1]                                    # Selecting the particular vertex (Only 1)
region1 = a.Set(vertices=FixVer1, name='Set-Fix-1')
mod.EncastreBC(name='Fix1', createStepName='Initial', region=region1, localCsys=None)            # Encastre BC
 
# Fixing Upper Left Node along X direction
UL_ver = v.findAt((0.0, round(BRVE, 6), 0.0))
qul = UL_ver.index
RollVer1 = v[qul:qul+1]
region2 = a.Set(vertices=RollVer1, name='Y_roll') 
mod.DisplacementBC(name='Y_roll', createStepName="Step-1", region=region2, u1=0.0, u2=UNSET, ur3=UNSET, amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', localCsys=None)

# Fix Bottom Right Node along Y direction
BR_ver = v.findAt((round(LRVE, 6), 0.0, 0.0))
qbr = BR_ver.index
RollVer2 = v[qbr:qbr+1]
region3 = a.Set(vertices=RollVer2, name='X_roll')
mod.DisplacementBC(name='X_roll', createStepName="Step-1", region=region3, u1=1, u2=0.0, ur3=UNSET, amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', localCsys=None) 