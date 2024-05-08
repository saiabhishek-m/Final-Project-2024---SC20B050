import numpy as np

# ABAQUS directories
from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *                                 ##################
from mesh import *                     
from optimization import *                         ## Code for creating the RVE with varying waviness level and thickness (Symmetric overlap)
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *


# Model details
modelname = 'Model-1'
partname = 'Part-1'
sheetsize = 5000
gridspace = 0.5


# Geometric details of RVE
l, we = 5.0e-3, 1.0e-3                           # length and width of platelet
wf = 0.6e-3                                   # width of the platelet at the centre
t = 0.163e-3                                  # thickness of matrix
x = 0.5*(we-wf)*(1+t/l)                    # a length factor used for modelling (necessary)

LRVE = l + t                               # length of RVE
BRVE = 2*(wf+x+t)                          # breadth of RVE


mod = mdb.models[modelname]
mod.ConstrainedSketch(name='__profile__', sheetSize=sheetsize, ).setPrimaryObject(option=STANDALONE)
s = mod.sketches['__profile__']
s.rectangle(point1=(0.0, 0.0), point2=(LRVE, BRVE))                                                   # Creating basic rectangular geometry (LRVE x BRVE)
session.viewports['Viewport: 1'].view.fitView()
mod.Part(name=partname, dimensionality=TWO_D_PLANAR, type=DEFORMABLE_BODY)
mod.parts[partname].BaseShell(sketch=mod.sketches['__profile__'])
mod.ConstrainedSketch(name='__profile__', sheetSize=sheetsize).unsetPrimaryObject()                   # Creating a part

# Creating separate partitions to make the microstructure 
session.viewports['Viewport: 1'].setValues(displayedObject=mod.parts[partname])
p = mod.parts[partname]
f = p.faces
RVEcent = (0.0 ,0.0, 0.0)                                                                       # Temporary origin - may or may not be required. all zero means no change. Can be useful sometimes for modelling
tr = p.MakeSketchTransform(sketchPlane=f[0], sketchPlaneSide=SIDE1, origin=RVEcent)             # Transforming the temporary origin as main origin to create partitions
                                                                                                # Temporary origin is defined with respect to the origin used for part creation

s = mod.ConstrainedSketch(name='__profile__', sheetSize=sheetsize, gridSpacing=gridspace, transform=tr)
s.setPrimaryObject(option=SUPERIMPOSE)
p = mod.parts[partname]
p.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)

# Creating partitions using lines and coordinates. 
s.Line(point1=(0-RVEcent[0],wf/2-RVEcent[1]),point2=(l/2-RVEcent[0],we/2-RVEcent[1]))                  #HO      (HO - name of line used in modelling)
s.Line(point1=(l/2-RVEcent[0],we/2-RVEcent[1]),point2=(l/2-RVEcent[0],0-RVEcent[1]))             #OA
s.Line(point1=(t+l/2-RVEcent[0],0-RVEcent[1]),point2=(l/2+t-RVEcent[0],0.5*we-RVEcent[1]))      #BP
s.Line(point1=(l/2+t-RVEcent[0],0.5*we-RVEcent[1]),point2=(LRVE-RVEcent[0],0.5*wf-RVEcent[1]))         #PC
s.Line(point1=(l/2+t-RVEcent[0],BRVE-0.5*we-RVEcent[1]),point2=(LRVE-RVEcent[0],BRVE-wf/2-RVEcent[1]))                   #QD
s.Line(point1=(l/2+t-RVEcent[0],BRVE-RVEcent[1]),point2=(l/2+t-RVEcent[0],BRVE-we/2-RVEcent[1]))            #EQ
s.Line(point1=(l/2-RVEcent[0],BRVE-RVEcent[1]),point2=(l/2-RVEcent[0],BRVE-we/2-RVEcent[1]))             #FR
s.Line(point1=(l/2-RVEcent[0],BRVE-we/2-RVEcent[1]),point2=(0-RVEcent[0],BRVE-wf/2-RVEcent[1]))              #RG
s.Line(point1=(t/2-RVEcent[0],0.5*(BRVE+we)-RVEcent[1]),point2=(LRVE/2-RVEcent[0],0.5*(BRVE+wf)-RVEcent[1]))            #NM
s.Line(point1=(LRVE/2-RVEcent[0],0.5*(BRVE+wf)-RVEcent[1]),point2=(LRVE-t/2-RVEcent[0],0.5*(BRVE+we)-RVEcent[1]))          #ML 
s.Line(point1=(LRVE-t/2-RVEcent[0],0.5*(BRVE+we)-RVEcent[1]),point2=(LRVE-t/2-RVEcent[0],0.5*(BRVE-we)-RVEcent[1]))                   #LK
s.Line(point1=(LRVE-t/2-RVEcent[0],BRVE/2-0.5*we-RVEcent[1]),point2=(LRVE/2-RVEcent[0],BRVE/2-0.5*wf-RVEcent[1]))                   #KJ
s.Line(point1=(LRVE/2-RVEcent[0],BRVE/2-0.5*wf-RVEcent[1]),point2=(t/2-RVEcent[0],BRVE/2-0.5*we-RVEcent[1]))            #JI
s.Line(point1=(t/2-RVEcent[0],BRVE/2-0.5*we-RVEcent[1]),point2=(t/2-RVEcent[0],BRVE/2+we/2-RVEcent[1]))       #IN

s.Line(point1=(t/2-RVEcent[0],(BRVE-we)/2-RVEcent[1]),point2=(t/2-RVEcent[0],(BRVE-we)/2-t-RVEcent[1]))                 #IS
s.Line(point1=(l/2-RVEcent[0],we/2+t-RVEcent[1]),point2=(l/2-RVEcent[0],0.5*we-RVEcent[1]))                 #WO
s.Line(point1=(l/2+t-RVEcent[0],we/2+t-RVEcent[1]),point2=(l/2+t-RVEcent[0],we/2-RVEcent[1]))             #XP
s.Line(point1=(LRVE-t/2-RVEcent[0],0.5*(BRVE-we)-RVEcent[1]),point2=(LRVE-t/2-RVEcent[0],0.5*(BRVE-we)-t-RVEcent[1]))            #KT
s.Line(point1=(LRVE-t/2-RVEcent[0],BRVE/2+0.5*we+t-RVEcent[1]),point2=(LRVE-t/2-RVEcent[0],BRVE/2+0.5*we-RVEcent[1]))        #UL
s.Line(point1=(l/2+t-RVEcent[0],BRVE-we/2-RVEcent[1]),point2=(l/2+t-RVEcent[0],BRVE-we/2-t-RVEcent[1]))             #QY
s.Line(point1=(l/2-RVEcent[0],BRVE-we/2-RVEcent[1]),point2=(l/2-RVEcent[0],BRVE-we/2-t-RVEcent[1]))          #RZ
s.Line(point1=(t/2-RVEcent[0],0.5*(BRVE+we)+t-RVEcent[1]),point2=(t/2-RVEcent[0],0.5*(BRVE+we)-RVEcent[1]))         #VN


p = mod.parts[partname]

# Creating faces for the generated partitions (Points defined from the real origin and not from the transformed origin)
pt1 = (l/4,0.25*wf,0)
f1 = f.findAt(pt1)
pt2 = (LRVE-l/4,0.25*wf,0)                  
f2 = f.findAt(pt2)                                            
pt3 = (LRVE/2,0.5*BRVE,0)
f3 = f.findAt(pt3)
pt4 = (LRVE-l/4,BRVE-0.25*wf,0)
f4 = f.findAt(pt4) 
pt5 = (0.25*l,BRVE-0.25*wf,0)
f5 = f.findAt(pt5)  
pt6 = (t/4,BRVE/2,0)
f6 = f.findAt(pt6)
pt7 = (l/4,0.25*(we+wf+t),0)
f7 = f.findAt(pt7)
pt8 = (LRVE/2,0.25*we,0)
f8 = f.findAt(pt8)
pt9 = (LRVE-l/4,0.25*(we+wf+t),0)
f9 = f.findAt(pt9)
pt10 = (LRVE-t/4,BRVE/2,0)
f10 = f.findAt(pt10)
pt11 = (LRVE-0.25*l,BRVE-0.25*(we+wf+t),0)
f11 = f.findAt(pt11)
pt12 = (LRVE/2,BRVE-we/4,0)
f12 = f.findAt(pt12)
pt13 = (we/4,BRVE-0.25*(we+wf+t),0)
f13 = f.findAt(pt13)


pickedfaces = (f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13 )
p.PartitionFaceBySketch(faces=pickedfaces, sketch=s)
s.unsetPrimaryObject()                                                             # Final creation of part
