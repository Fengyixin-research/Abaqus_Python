# -*- coding: mbcs -*-
# Do not delete the following import lines
from abaqus import *
from abaqusConstants import *
import __main__

def ViewportFormat():
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
    session.viewports['Viewport-Model-1'].odbDisplay.display.setValues(plotState=(
        CONTOURS_ON_DEF, ))
    session.viewports['Viewport-Model-1'].view.setValues(nearPlane=192.175, 
        farPlane=317.444, width=86.5755, height=36.117, viewOffsetX=0.603315, 
        viewOffsetY=5.29618)
    session.viewports['Viewport-Model-1'].odbDisplay.display.setValues(plotState=(
        CONTOURS_ON_UNDEF, CONTOURS_ON_DEF, ))
    session.viewports['Viewport-Model-1'].odbDisplay.commonOptions.setValues(
        visibleEdges=NONE)
    session.viewports['Viewport-Model-1'].odbDisplay.superimposeOptions.setValues(
        visibleEdges=NONE)
    session.viewports['Viewport-Model-1'].odbDisplay.setPrimaryVariable(
        variableLabel='U', outputPosition=NODAL, refinement=(INVARIANT, 
        'Magnitude'))
    session.viewports['Viewport-Model-1'].viewportAnnotationOptions.setValues(
        triad=OFF, title=OFF, state=OFF, compass=OFF)
    session.viewports['Viewport-Model-1'].viewportAnnotationOptions.setValues(
        legendFont='-*-verdana-medium-r-normal-*-*-80-*-*-p-*-*-*')
    session.viewports['Viewport-Model-1'].viewportAnnotationOptions.setValues(
        legendBox=OFF)


def ViewportOrigin():
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
    session.viewports['Viewport: 1'].view.setValues(nearPlane=192.175, 
        farPlane=317.444, width=86.5755, height=36.117, viewOffsetY=8.45917)


