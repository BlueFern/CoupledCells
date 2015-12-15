import re
import os
import sys
import glob

def natural_key(string_):
    """See http://www.codinghorror.com/blog/archives/001018.html"""
    return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', string_)]

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

cwd = os.getcwd()
print cwd

glob_pattern = os.path.join(cwd, 'solutionVTU/smc*.vtu')
print glob_pattern

files = sorted(glob.glob(glob_pattern), key=natural_key)
print files

# create a new 'XML Unstructured Grid Reader'
smc_data_t_ = XMLUnstructuredGridReader(FileName=files)
smc_data_t_.CellArrayStatus = ['SMC_Ca']

# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
renderView1.OrientationAxesVisibility = 0
renderView1.CenterAxesVisibility = 0

# uncomment following to set a specific view size
renderView1.ViewSize = [1482, 847]

# show data in view
smc_data_t_Display = Show(smc_data_t_, renderView1)
# trace defaults for the display properties.
smc_data_t_Display.ColorArrayName = [None, '']
smc_data_t_Display.ScalarOpacityUnitDistance = 230.0

# reset view to fit data
renderView1.ResetCamera()

# create a new 'Annotate Time Filter'
annotateTimeFilter1 = AnnotateTimeFilter(Input=smc_data_t_)
annotateTimeFilter1.Format = 't = %.2f'

# show data in view
annotateTimeFilter1Display = Show(annotateTimeFilter1, renderView1)
# trace defaults for the display properties.
annotateTimeFilter1Display.FontSize = 12

# Properties modified on smc_data_t_
smc_data_t_.CellArrayStatus = ['SMC_Ca']

# set scalar coloring
ColorBy(smc_data_t_Display, ('CELLS', 'SMC_Ca'))

# rescale color and/or opacity maps used to include current data range
smc_data_t_Display.RescaleTransferFunctionToDataRange(True)

# show color bar/color legend
smc_data_t_Display.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'SMCCa'
sMCCaLUT = GetColorTransferFunction('SMCCa')
sMCCaLUT.RGBPoints = [0.011199812077250126, 0.231373, 0.298039, 0.752941, 0.011200050109459181, 0.865003, 0.865003, 0.865003, 0.011200288141668235, 0.705882, 0.0156863, 0.14902]
sMCCaLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'SMCCa'
sMCCaPWF = GetOpacityTransferFunction('SMCCa')
sMCCaPWF.Points = [0.011199812077250126, 0.0, 0.5, 0.0, 0.011200288141668235, 1.0, 0.5, 0.0]
sMCCaPWF.ScalarRangeInitialized = 1

# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
sMCCaLUT.ApplyPreset('jet', True)

# Rescale transfer function
sMCCaLUT.RescaleTransferFunction(0.0, 1.0)

# Rescale transfer function
sMCCaPWF.RescaleTransferFunction(0.0, 1.0)

# set active source
SetActiveSource(annotateTimeFilter1)

# current camera placement for renderView1
#renderView1.CameraPosition = [8750.0, -2300.0, 26400.0]
#renderView1.CameraFocalPoint = [8750.0, -2300.0, -16050.0]
#renderView1.CameraParallelScale = 11000.0

renderView1.CameraPosition = [7279.867724733483, 56.878043635632764, 7108.638010203612]
renderView1.CameraFocalPoint = [7279.867724733483, 56.878043635632764, 0.0]
renderView1.CameraParallelScale = 10229.418883053055

# Properties modified on animationScene1
animationScene1.PlayMode = 'Sequence'

# Properties modified on animationScene1
animationScene1.StartTime = 0

# Properties modified on animationScene1
animationScene1.EndTime = 1000

# Properties modified on animationScene1
animationScene1.NumberOfFrames = 900

# save animation images/movie
WriteAnimation(os.path.join(cwd, 'Animation_SMC_Ca.avi'), Magnification=1, FrameRate=15.0, Compression=True)

RenderAllViews()




