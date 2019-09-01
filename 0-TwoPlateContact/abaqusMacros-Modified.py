#coding:utf-8
from abaqus import *
import testUtils
testUtils.setBackwardCompatibility()
from abaqusConstants import *
import __main__
import section
import regionToolset
import displayGroupMdbToolset as dgm
import part
import material
import assembly
import step
import interaction
import load
import mesh
import optimization
import job
import sketch
import visualization
import xyPlot
import displayGroupOdbToolset as dgo
import connectorBehavior

### Create a new model named 'Model-1'
myModel = mdb.Model(name='Model-1')

### Create a constrained sketch for 'Part-PlateUpper' 
s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=200.0)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=STANDALONE)
s.rectangle(point1=(-10.0, 1.0), point2=(10.0, 0.0))
# Create part 'Part-PlateUpper'
p = mdb.models['Model-1'].Part(name='Part-PlateUpper', dimensionality=TWO_D_PLANAR, type=DEFORMABLE_BODY)
p = mdb.models['Model-1'].parts['Part-PlateUpper']
p.BaseShell(sketch=s)
s.unsetPrimaryObject()
del mdb.models['Model-1'].sketches['__profile__']


### Create a constrained sketch for 'Part-PlateLower' ###
s1 = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=200.0)
g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
s1.setPrimaryObject(option=STANDALONE)
s1.rectangle(point1=(-10.0, 0.0), point2=(10.0, -1.0))
p = mdb.models['Model-1'].Part(name='Part-PlateLower', dimensionality=TWO_D_PLANAR, type=DEFORMABLE_BODY)
p = mdb.models['Model-1'].parts['Part-PlateLower']
p.BaseShell(sketch=s1)
s1.unsetPrimaryObject()
del mdb.models['Model-1'].sketches['__profile__']


### Create material named 'Material-PlateUpper'
mdb.models['Model-1'].Material(name='Material-PlateUpper')
mdb.models['Model-1'].materials['Material-PlateUpper'].Elastic(table=((210000.0, 0.3), ))
mdb.models['Model-1'].materials['Material-PlateUpper'].Expansion(table=((0.3, ), ))
### Create material named 'Material-PlateLower'
mdb.models['Model-1'].Material(name='Material-PlateLower')
mdb.models['Model-1'].materials['Material-PlateLower'].Elastic(table=((210000.0, 0.3), ))
mdb.models['Model-1'].materials['Material-PlateLower'].Expansion(table=((0.1, ), ))
### Create HomogeneousSolidSection objects and assign material property to them ###
mdb.models['Model-1'].HomogeneousSolidSection(name='Section-PlateUpper', material='Material-PlateUpper', 
    thickness=None)
mdb.models['Model-1'].HomogeneousSolidSection(name='Section-PlateLower', material='Material-PlateLower', 
    thickness=None)


### Assign 'Section-PlateUpper' to 'Part-PlateUpper' 
p = mdb.models['Model-1'].parts['Part-PlateUpper']
f = p.faces
faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
# Create a region that represents 'Part-PlateUpper'
region = regionToolset.Region(faces=faces)
p.SectionAssignment(region=region, sectionName='Section-PlateUpper', 
    offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)
### Assign 'Section-PlateLower' to 'Part-PlateLower' ###
p = mdb.models['Model-1'].parts['Part-PlateLower']
f = p.faces
faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
region = regionToolset.Region(faces=faces)
p.SectionAssignment(region=region, sectionName='Section-PlateLower', 
    offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)
mdb.save()


### Create assembly consisting of two plate instances
a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByDefault(CARTESIAN)
p = mdb.models['Model-1'].parts['Part-PlateLower']
# Create an instance named 'Part-PlateLower-1'
a.Instance(name='Part-PlateLower-1', part=p, dependent=ON)
p = mdb.models['Model-1'].parts['Part-PlateUpper']
# Create an instance named 'Part-PlateUpper-1'
a.Instance(name='Part-PlateUpper-1', part=p, dependent=ON)
# Make both instances independent
a.makeIndependent(instances=(a.instances['Part-PlateLower-1'], ))
a.makeIndependent(instances=(a.instances['Part-PlateUpper-1'], ))


### Generate mesh for the assembly
# Create tuple of part instances
partInstances =(a.instances['Part-PlateLower-1'], a.instances['Part-PlateUpper-1'], )
# Seed part instances
a.seedPartInstance(regions=partInstances, size=0.4, deviationFactor=0.1, minSizeFactor=0.1)
# Generate mesh for partInstances based on previous seeding
a.generateMesh(regions=partInstances)
# Create a new step named 'Step-ThermalExpansion' after initial step
mdb.models['Model-1'].StaticStep(name='Step-ThermalExpansion', previous='Initial')


### Create interaction edges and define interaction relations 
s1 = a.instances['Part-PlateLower-1'].edges
side1Edges1 = s1.getSequenceFromMask(mask=('[#1 ]', ), )
# Create a surface named 'Surf-PlateLower'
a.Surface(side1Edges=side1Edges1, name='Surf-PlateLower')
s1 = a.instances['Part-PlateUpper-1'].edges
# Define number of the sideedge
side1Edges1 = s1.getSequenceFromMask(mask=('[#4 ]', ), )
# Create a surface named 'Surf-PlateUpper'
a.Surface(side1Edges=side1Edges1, name='Surf-PlateUpper')
region1=a.surfaces['Surf-PlateUpper']
region2=a.surfaces['Surf-PlateLower']
# Create tie surface contact between upper and lower surface of PlateLower and PlateUpper
mdb.models['Model-1'].Tie(name='Constraint-PlateContact', master=region1, 
    slave=region2, positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, thickness=ON)


### Create boundary conditions for the assembly
a = mdb.models['Model-1'].rootAssembly
e1 = a.instances['Part-PlateUpper-1'].edges
edges1 = e1.getSequenceFromMask(mask=('[#2 ]', ), )
# Create a set named 'Set-PlateUpperFix'
a.Set(edges=edges1, name='Set-PlateUpperFix')
e1 = a.instances['Part-PlateLower-1'].edges
edges1 = e1.getSequenceFromMask(mask=('[#2 ]', ), )
# Create a set named 'Set-PlateLowerFix'
a.Set(edges=edges1, name='Set-PlateLowerFix')
# Set boundary conditions
region = a.sets['Set-PlateUpperFix']
mdb.models['Model-1'].DisplacementBC(name='BC-PlateUpperFix', 
    createStepName='Initial', region=region, u1=SET, u2=UNSET, ur3=SET, 
    amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
region = a.sets['Set-PlateLowerFix']
mdb.models['Model-1'].DisplacementBC(name='BC-PlateLowerFix', 
    createStepName='Initial', region=region, u1=SET, u2=UNSET, ur3=SET, 
    amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)


### Define ambient temperature field
a = mdb.models['Model-1'].rootAssembly
f1 = a.instances['Part-PlateLower-1'].faces
# Get the face of 'Part-PlateLower-1'
faces1 = f1.getSequenceFromMask(mask=('[#1 ]', ), )
f2 = a.instances['Part-PlateUpper-1'].faces
# Get the face of 'Part-PlateUpper-1'
faces2 = f2.getSequenceFromMask(mask=('[#1 ]', ), )
# Create a set consisting of all faces
a.Set(faces=faces1+faces2, name='Set-PlateWhole')
region = a.sets['Set-PlateWhole']
mdb.models['Model-1'].Temperature(name='Predefined Field-Temperature', 
    createStepName='Initial', region=region, distributionType=UNIFORM, 
    crossSectionDistribution=CONSTANT_THROUGH_THICKNESS, magnitudes=(20.0, ))
mdb.models['Model-1'].predefinedFields['Predefined Field-Temperature'].setValuesInStep(
    stepName='Step-ThermalExpansion', magnitudes=(120.0, ))


### Create job and submit
jobName = 'Job-TwoPlateContact'
mdb.Job(name=jobName, model='Model-1', description='', 
    type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0, queue=None, 
    memory=90, memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
    explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
    modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', scratch='', resultsFormat=ODB)
mdb.jobs['Job-TwoPlateContact'].submit(consistencyChecking=OFF)
mdb.jobs['Job-TwoPlateContact'].waitForCompletion()
print 'The analysis has been completed. On going postprocessing'


### Visualization
myOdb = visualization.openOdb(path = jobName + '.odb')
myViewport = session.Viewport(name='Viewport-Model-1', origin =(0.0, 0.0), width=150, height=120)
myViewport.setValues(displayedObject = myOdb)
myViewport.odbDisplay.display.setValues(plotState = CONTOURS_ON_DEF)
myViewport.odbDisplay.commonOptions.setValues(renderStyle = FILLED)
session.printToFile(fileName = 'Mises', format = PNG, canvasObjects = (myViewport, ))
print 'File Mises.png has been saved in the working directory. Please check!'


