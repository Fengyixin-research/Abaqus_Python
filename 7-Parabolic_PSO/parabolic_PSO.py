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
import random
import operator
import numpy
import math
from odbAccess import *
from deap import base
from deap import creator
from deap import tools
from deap import benchmarks

### Dimension of a single mesh element
lengthMesh = 1.0
### Number of mesh elements along each dimension of a single deformable element
nx, ny = 3, 1
### Dimensions of a single deformable element
xlength = nx * lengthMesh; ylength = ny * lengthMesh
### Number of deformable elements along each dimension of the global beam
nX, nY = 20, 3
# Total number of genes in the chromosome is nX*nY
sizePopulation, numIter = 40, 20
### Thermal expansion coefficients
expandActive, expandInactive = 0.0005, 0.000

### Register 'FitnessMax' property in toolbox
creator.create("FitnessMax", base.Fitness, weights=(-1.0,))
creator.create("Particle", list, fitness=creator.FitnessMax, speed=list, smin=None, smax=None, best=None)

# Position limit of particles in each dimension
lower, upper = 0, 5
# Dimension of particle
dimension = nX * nY
# Speed limit of particles
slimit = 0.1 * (upper-lower)
# Acceleration constants
Phi1, Phi2 = 2.0, 2.0
# Offset of the origin of the viewport
offset = 1.5
# iter & nmodel
iter, nmodel = None, None
### Stores the maximum Mises of the entire population in all iterations 
### fitness_record contains all the fitness values in the current population at each iteration
fitness_pop = []
fitness_record = []
Mises_max=[]
Mises_max_record=[]
best_individual = []

### Function to generate particles
def generate(size, pmin, pmax, smin, smax):
    # size: length of the particle, 
    part = creator.Particle(random.uniform(pmin, pmax) for _ in range(size)) 
    part.speed = [random.uniform(smin, smax) for _ in range(size)]
    # Speed lower and upper limit
    part.smin = smin
    part.smax = smax
    return part


### Function to update particle's speed and position
def updateParticle(part, best, phi1, phi2):
    # First compute the speed
    # then limit the speed values between smin and smax
    # random.uniform(a, b): generate random float number, including a, excluding b
    u1 = (random.uniform(0, phi1) for _ in range(len(part)))
    u2 = (random.uniform(0, phi2) for _ in range(len(part)))
    v_u1 = map(operator.mul, u1, map(operator.sub, part.best, part))
    v_u2 = map(operator.mul, u2, map(operator.sub, best, part))
    part.speed = list(map(operator.add, part.speed, map(operator.add, v_u1, v_u2)))
    for i, speed in enumerate(part.speed):
        if abs(speed) < part.smin:
            part.speed[i] = math.copysign(part.smin, speed)
        elif abs(speed) > part.smax:
            part.speed[i] = math.copysign(part.smax, speed)
    # Update the parts speed
    part[:] = list(map(operator.add, part, part.speed))
    #################################################
    # Limit the position of part positions
    for i, p in enumerate(part):
        if p < float(lower):
            part[i] = float(lower)
        elif p > float(upper)+0.9999:
            part[i] = float(upper)+0.9999
    #################################################
            

### Objective function that returns deviation of simulated shape from parabolic shape
def MinimizeDeviationFromParabolic(individual):
    # Declare variables as global
    global lengthMesh, nx, ny, nX, nY, xlength, ylength, iter, nmodel
    ### Round each element of individual to integer by eliminating decimal
    individual_int = [int(i) for i in individual]

    ### Create a new model named 'Model-Iter0-Num0'
    modelName = 'Model-Iter'+str(iter)+'-Num'+str(nmodel)
    myModel = mdb.Model(name= modelName)
    ### Create constrained sketch for all deformable elements
    for i in range(nX):
        for j in range(nY):
            s = mdb.models[modelName].ConstrainedSketch(name='__profile__', sheetSize=xlength*nX*1.2)
            s.setPrimaryObject(option=STANDALONE)
            s.rectangle(point1=(i*xlength, j*ylength), point2=((i+1)*xlength, (j+1)*ylength))
            # Create part named 'Part-Plate-i-j'
            partName = 'Part-Plate-'+str(i)+'-'+str(j)
            p = mdb.models[modelName].Part(name=partName, dimensionality=TWO_D_PLANAR, type=DEFORMABLE_BODY)
            p = mdb.models[modelName].parts[partName]
            p.BaseShell(sketch=s)
            s.unsetPrimaryObject()
            # Define vertices sets
            v = p.vertices
            location1 = (i*xlength, j*ylength, 0)  # bottom-left
            location2 = (i*xlength, (j+1)*ylength, 0)   # top-left
            location3 = ((i+1)*xlength, (j+1)*ylength, 0)   # top-right
            location4 = ((i+1)*xlength, j*ylength, 0)   # bottom-right
            verts1 = v.findAt((location1, ))
            verts2 = v.findAt((location2, ))
            verts3 = v.findAt((location3, ))
            verts4 = v.findAt((location4, ))
            p.Set(vertices=verts1, name='Set-1')
            p.Set(vertices=verts2, name='Set-2')
            p.Set(vertices=verts3, name='Set-3')
            p.Set(vertices=verts4, name='Set-4')
            del mdb.models[modelName].sketches['__profile__']
    
    ### Create multiple material properties
    # lower-upper: inactive-active
    for i in range(upper-lower+1):
        materialName = 'Material-'+str(i)
        mdb.models[modelName].Material(name=materialName)
        mdb.models[modelName].materials[materialName].Elastic(table=((210000.0, 0.3), ))
        expandCoefficient=(float(i-lower)/(upper-lower))*expandInactive+(1-float(i-lower)/(upper-lower))*expandActive
        mdb.models[modelName].materials[materialName].Expansion(table=((expandCoefficient, ), )) 
    ### Create multiple sections
    for i in range(upper-lower+1):
        sectionName = 'Section-'+str(i)
        materialName = 'Material-'+str(i)
        # Create solid section and assign material property
        mdb.models[modelName].HomogeneousSolidSection(name=sectionName, material=materialName, thickness=None)
    for i in range(nX):
        for j in range(nY):
            # Create partName named 'Part-Plate-i-j'
            partName = 'Part-Plate-'+str(i)+'-'+str(j)
            sectionName = 'Section-'+str(individual_int[j*nX+i])
            p = mdb.models[modelName].parts[partName]
            f = p.faces
            faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
            # Create a region
            region = regionToolset.Region(faces=faces)
            p.SectionAssignment(region=region, sectionName=sectionName, offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)

    ### Create assembly consisting of all plate instances
    a = mdb.models[modelName].rootAssembly
    a.DatumCsysByDefault(CARTESIAN)
    for i in range(nX):
        for j in range(nY):
            # Create partName named 'Part-Plate-i-j'
            partName = 'Part-Plate-'+str(i)+'-'+str(j)
            instanceName = partName+'-1'
            p = mdb.models[modelName].parts[partName]
            # Create indenpendent instances by setting 'dependent' to OFF
            a.Instance(name=instanceName, part=p, dependent=OFF)

    ### Generate mesh for the assembly
    # Create tuple of part instances
    partInstances = []
    for i in range(nX):
        for j in range(nY):
            instanceName = 'Part-Plate-'+str(i)+'-'+str(j)+'-1'
            partInstances.append(a.instances[instanceName])
    partInstances=tuple(partInstances)
    # Seed part instances
    a.seedPartInstance(regions=partInstances, size=lengthMesh, deviationFactor=0.1, minSizeFactor=0.1)
    # Generate mesh for partInstances based on previous seeding
    a.generateMesh(regions=partInstances)

    ### Create a new step named 'Step-ThermalExpansion' after initial step
    mdb.models[modelName].StaticStep(name='Step-ThermalExpansion', previous='Initial')

    ### Create interaction edges and define interaction relations 
    for i in range(nX):
        for j in range(nY):
            # Create instanceName named 'Part-Plate-i-j-1'
            instanceName = 'Part-Plate-'+str(i)+'-'+str(j)+'-1'
            s1 = a.instances[instanceName].edges
            edgeCenter1 = ((i+0.5)*xlength, (j+1.0)*ylength, 0)
            edgeCenter2 = ((i+0.0)*xlength, (j+0.5)*ylength, 0)
            edgeCenter3 = ((i+1.0)*xlength, (j+0.5)*ylength, 0)
            edgeCenter4 = ((i+0.5)*xlength, (j+0.0)*ylength, 0)
            edges1 = s1.findAt((edgeCenter1, ))
            edges2 = s1.findAt((edgeCenter2, ))
            edges3 = s1.findAt((edgeCenter3, ))
            edges4 = s1.findAt((edgeCenter4, ))
            edges1 = (edges1, ); edges2 = (edges2, ); edges3 = (edges3, ); edges4 = (edges4, )
            # Create surfaces for all 4 edges of current element
            # a is rootAssembly
            a.Surface(side1Edges=edges1, name='Surf-Plate-'+str(i)+'-'+str(j)+'-1')
            a.Surface(side1Edges=edges2, name='Surf-Plate-'+str(i)+'-'+str(j)+'-2')
            a.Surface(side1Edges=edges3, name='Surf-Plate-'+str(i)+'-'+str(j)+'-3')
            a.Surface(side1Edges=edges4, name='Surf-Plate-'+str(i)+'-'+str(j)+'-4')
    # Define interaction relations between lower and upper surfaces
    for i in range(nX):
        for j in range(nY-1):
            surfName1 = 'Surf-Plate-'+str(i)+'-'+str(j)+'-1'
            surfName2 = 'Surf-Plate-'+str(i)+'-'+str(j+1)+'-4'
            region1=a.surfaces[surfName1]
            region2=a.surfaces[surfName2]
            contactName = 'Constraint-'+str(i)+str(j)+str(i)+str(j+1)+'-14'
            mdb.models[modelName].Tie(name=contactName, master=region1, slave=region2, positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, thickness=ON)
    # Define interaction relations between left and right surfaces
    for i in range(nX-1):
        for j in range(nY):
            surfName1 = 'Surf-Plate-'+str(i)+'-'+str(j)+'-3'
            surfName2 = 'Surf-Plate-'+str(i+1)+'-'+str(j)+'-2'
            region1=a.surfaces[surfName1]
            region2=a.surfaces[surfName2]
            contactName = 'Constraint-'+str(i)+str(j)+str(i+1)+str(j)+'-32'
            mdb.models[modelName].Tie(name=contactName, master=region1, slave=region2, positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, thickness=ON)
    
    ### Create boundary conditions for the assembly
    a = mdb.models[modelName].rootAssembly
    for j in range(nY):
        instanceName = 'Part-Plate-0-'+str(j)+'-1'
        e1 = a.instances[instanceName].edges
        edgeCenter = (0.0, (j+0.5)*ylength, 0.0)
        edges1 = e1.findAt((edgeCenter, ))
        edges1 = (edges1, )
        setName = 'Set-Plate-0-'+str(j)+'-Fix'
        a.Set(edges=edges1, name=setName)
        # Set boundary conditions
        region = a.sets[setName]
        bcName = 'BC-Plate-0-'+str(j)+'-Fix'
        mdb.models[modelName].DisplacementBC(name=bcName, 
            createStepName='Initial', region=region, u1=SET, u2=UNSET, ur3=SET, 
            amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
    instanceName = 'Part-Plate-0-0-1'
    v1 = a.instances[instanceName].vertices
    vertexCenter = (0.0, 0.0, 0.0)
    v = v1.findAt((vertexCenter, ))
    v = (v, )
    setName = 'Set-Plate-Origin-Fix'
    a.Set(vertices=v, name=setName)
    region = a.sets[setName]
    bcName = 'BC-Origin-Fix'
    # Origin fix in x-direction, y-direction, xy-rotation
    mdb.models[modelName].DisplacementBC(name=bcName, 
            createStepName='Initial', region=region, u1=SET, u2=SET, ur3=SET, 
            amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)

    ### Define nodeSets that records the displacement information of the entire assembly
    # 1st: bottom left vertices of all deformable elements
    a = mdb.models[modelName].rootAssembly
    for i in range(nX):
        for j in range(nY):
            nodeName = 'Node-Plate-'+str(i)+'-'+str(j)
            instanceName = 'Part-Plate-'+str(i)+'-'+str(j)+'-1'
            v1 = a.instances[instanceName].vertices
            vertexCenter = (i*xlength, j*ylength, 0.0)
            v = v1.findAt((vertexCenter, ))
            v = (v, )
            a.Set(vertices=v, name=nodeName)
    # 2nd: top left vertices of elements in the top row
    for i in range(nX):
        nodeName = 'Node-Plate-'+str(i)+'-'+str(nY)
        instanceName = 'Part-Plate-'+str(i)+'-'+str(nY-1)+'-1'
        v1 = a.instances[instanceName].vertices
        vertexCenter = (i*xlength, nY*ylength, 0.0)
        v = v1.findAt((vertexCenter, ))
        v = (v, )
        a.Set(vertices=v, name=nodeName)
    # 3rd: bottom right vertices of elements in the rightmost column
    for j in range(nY):
        nodeName = 'Node-Plate-'+str(nX)+'-'+str(j)
        instanceName = 'Part-Plate-'+str(nX-1)+'-'+str(j)+'-1'
        v1 = a.instances[instanceName].vertices
        vertexCenter = (nX*xlength, j*ylength, 0.0)
        v = v1.findAt((vertexCenter, ))
        v = (v, )
        a.Set(vertices=v, name=nodeName)
    # 4th: top right vertex of the entire assembly
    nodeName = 'Node-Plate-'+str(nX)+'-'+str(nY)
    instanceName = 'Part-Plate-'+str(nX-1)+'-'+str(nY-1)+'-1'
    v1 = a.instances[instanceName].vertices
    vertexCenter = (nX*xlength, nY*ylength, 0.0)
    v = v1.findAt((vertexCenter, ))
    v = (v, )
    a.Set(vertices=v, name=nodeName)

    ### Define ambient temperature field
    a = mdb.models[modelName].rootAssembly
    allFaceSet = []
    for i in range(nX):
        for j in range(nY):
            instanceName = 'Part-Plate-'+str(i)+'-'+str(j)+'-1'
            f1 = a.instances[instanceName].faces
            faces = f1.getSequenceFromMask(mask=('[#1 ]', ), )
            allFaceSet.append(faces)
    allFaceSet =tuple(allFaceSet)
    setName = 'Set-allFace'
    a.Set(faces=allFaceSet, name=setName)
    region = a.sets[setName]
    # Set temperature for initial and final state
    initialT = 20.0
    finalT = 120.0
    mdb.models[modelName].Temperature(name='Predefined Field-Temperature', 
        createStepName='Initial', region=region, distributionType=UNIFORM, 
        crossSectionDistribution=CONSTANT_THROUGH_THICKNESS, magnitudes=(initialT, ))
    mdb.models[modelName].predefinedFields['Predefined Field-Temperature'].setValuesInStep(
        stepName='Step-ThermalExpansion', magnitudes=(finalT, ))

    ### Create job and submit. Set parallel computing
    jobName = 'Job-Iter'+str(iter)+'-Num'+str(nmodel)
    mdb.Job(name=jobName, model=modelName, description='', type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0, queue=None, 
        memory=90, memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
        modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', scratch='', resultsFormat=ODB, parallelizationMethodExplicit=DOMAIN, 
        numDomains=8, activateLoadBalancing=False, multiprocessingMode=THREADS, numCpus=8)     
    mdb.jobs[jobName].submit(consistencyChecking=OFF)
    mdb.jobs[jobName].waitForCompletion()
    print 'The analysis has been completed. On going postprocessing'

    ### Visualization
    myOdb = visualization.openOdb(path = jobName + '.odb')
    myViewport = session.Viewport(name='Viewport-'+modelName, origin =(0.0, 0.0), width=240, height=90)
    myViewport.setValues(displayedObject = myOdb)
    myViewport.view.setValues(viewOffsetY=offset*nY*ylength)  # Move center of the viewport upward
    myViewport.odbDisplay.display.setValues(plotState=(CONTOURS_ON_UNDEF, CONTOURS_ON_DEF, ))
    myViewport.odbDisplay.commonOptions.setValues(visibleEdges=NONE)
    myViewport.odbDisplay.superimposeOptions.setValues(visibleEdges=NONE)
    myViewport.odbDisplay.commonOptions.setValues(renderStyle = FILLED)
    myViewport.odbDisplay.setPrimaryVariable(variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U2'))
    myViewport.viewportAnnotationOptions.setValues(triad=OFF, title=OFF, state=OFF, compass=OFF)
    myViewport.viewportAnnotationOptions.setValues(legendFont='-*-verdana-medium-r-normal-*-*-80-*-*-p-*-*-*')
    myViewport.viewportAnnotationOptions.setValues(legendBox=OFF)
    fileName = 'U2-Iter'+str(iter)+'-Num'+str(nmodel)
    session.printToFile(fileName = fileName, format = PNG, canvasObjects = (myViewport, ))
    print 'File ' + fileName + ' has been saved in the working directory. Please check!'

    ### Output nodal displacements
    odb = openOdb(path = jobName + '.odb')
    assembly = odb.rootAssembly
    # First step + Last frame of the simulation
    step1 = odb.steps['Step-ThermalExpansion']
    lastFrame = step1.frames[-1]

    ##########################################
    # Define functions to calculate RMSE
    def SquaredError(xy):
        return (Shape(xy[0])-xy[1])**2

    # Defines the target shape
    def Shape(x):
        return (5*nY*ylength/((nX*xlength)**2))*x**2

    coordinates = []
    for i in range(nX+1):
        nodeName = 'NODE-PLATE-'+str(i)+'-'+str(0)
        nodes = assembly.nodeSets[nodeName]
        displacement = lastFrame.fieldOutputs['U']
        nodesDisplacement = displacement.getSubset(region = nodes)
        # coordinates stores the target coordinates of bottom nodes
        coordinates.append([i*xlength+nodesDisplacement.values[0].data[0], nodesDisplacement.values[0].data[1]])
    # Root Mean Square Error
    # Square root of sum of squares over number of individuals
    RMSE = (sum(map(SquaredError, coordinates))**0.5)/len(coordinates)

    ###################################################
    ### Read Mises
    maxValue = None
    stressOutputExists = False
    for step in odb.steps.values():
        print 'Searching steps: ', step.name
        for frame in step.frames:
            try:
                stress = frame.fieldOutputs['S']
                stressOutputExists=True
            except KeyError: 
                continue
            for stressValue in stress.values:
                if (not maxValue or stressValue.mises > maxValue.mises):
                    maxValue = stressValue
    Mises_max_record.append(maxValue.mises)
    ###################################################
    ### Must close current odb object in order to begin the next iteration
    odb.close()
    del mdb.models[modelName]
    del mdb.jobs[jobName]
    return RMSE,
    ##########################################

### Toolbox
toolbox = base.Toolbox()
# size=2: means the particle is in a 2-dimensional space
toolbox.register("particle", generate, size=dimension, pmin=float(lower), pmax=float(upper)+0.9999, smin=-slimit, smax=slimit)
toolbox.register("population", tools.initRepeat, list, toolbox.particle)
# acceleration constants
toolbox.register("update", updateParticle, phi1=Phi1, phi2=Phi2)
# toolbox.register("evaluate", benchmarks.h1)
toolbox.register("evaluate", MinimizeDeviationFromParabolic)

# variable best contains the best candidate ever found
# Size of the population is 5
pop = toolbox.population(n=sizePopulation)
stats = tools.Statistics(lambda ind: ind.fitness.values)
stats.register("avg", numpy.mean)
stats.register("std", numpy.std)
stats.register("min", numpy.min)
stats.register("max", numpy.max)

logbook = tools.Logbook()
logbook.header = ["gen", "evals"] + stats.fields

GEN = numIter
best = None

for g in range(GEN):
    global iter, nmodel
    iter = g
    fitness_record = []
    Mises_max_record=[]
    for i, part in enumerate(pop):
        nmodel = i 
        # Evaluate the fitness values of each object of the population
        part.fitness.values = toolbox.evaluate(part)

        ### Append fitness value of current particle
        fitness_record.append(part.fitness.values[0])

        # Finding the individual with minimum fitness value
        if not part.best or part.best.fitness < part.fitness:
            part.best = creator.Particle(part)
            part.best.fitness.values = part.fitness.values
        if not best or best.fitness < part.fitness:
            best = creator.Particle(part)
            best.fitness.values = part.fitness.values

    ### Stores the best individual in the current generation
    best_record = [int(i) for i in best]
    best_individual.append(best_record)
    fitness_pop.append(fitness_record)
    Mises_max.append(Mises_max_record)

    # Update the population with all parts and the best individual so far
    for part in pop:
        toolbox.update(part, best)
    
    # Gather all the fitnesses in one list and print the stats
    logbook.record(gen=g, evals=len(pop), **stats.compile(pop))
    print(logbook.stream)


import csv
import os
with open("fitness_pop.csv", "wb") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerows(fitness_pop)

with open("best_individual.csv", "wb") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['nX','nY','sizePopulation','numIter','expandActive','expandInactive','lower','upper'])
    writer.writerow([nX, nY, sizePopulation, numIter, expandActive, expandInactive, lower, upper])
    writer.writerows(best_individual)

with open("Mises_max.csv", "wb") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerows(Mises_max)

