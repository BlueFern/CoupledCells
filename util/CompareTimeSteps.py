# -*- coding: utf-8 -*-
"""
Created on Mon Aug  3 15:24:30 2015

@author: sed59
"""
from vtk import *
import sys
import os
import re
import glob

iren = vtk.vtkRenderWindowInteractor()
mapper_0 = vtkDataSetMapper()
mapper_1 = vtkDataSetMapper()
reader_0 = vtkXMLUnstructuredGridReader()
reader_1 = vtkXMLUnstructuredGridReader()
base_files = []
test_files = []
time_step = 0
max_time_step = 99
min_time_step = 0

def update_mappers(filename_0, filename_1):
    """Read in two files and set them as input to the corresponding mappers."""
    
    reader_0.SetFileName(base_files[time_step])
    reader_0.Update()
    output_0 = reader_0.GetOutput()
    output_0.GetCellData().SetActiveScalars(sys.argv[4])
    mapper_0.SetInput(output_0)
    
    reader_1.SetFileName(filename_0)
    reader_1.Update()
    output_1 = reader_1.GetOutput()
    output_1.GetCellData().SetActiveScalars(sys.argv[4])
    mapper_1.SetInput(output_1)

def sortNicely(l): 
    """ Sort the given list in the way that humans expect."""
     
    convert = lambda text: int(text) if text.isdigit() else text 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    l.sort( key=alphanum_key )


def keypressCallbackFunction(obj, event):
    global time_step
    key= iren.GetKeySym()
    if (key == "Right" and time_step < max_time_step):
        time_step += 1
        
    elif (key == "Left" and time_step > min_time_step):
        time_step -= 1
    
    if time_step < len(base_files) and time_step < len(test_files):
        update_mappers(base_files[time_step],test_files[time_step])
        iren.Render()   

def main():
    global base_files, test_files
    
    if len(sys.argv) != 5:  
        sys.exit("Expected arguments: <base data dir> <new time-steps dir> <ec|smc> <attribute name>")
    
    render_window = vtk.vtkRenderWindow()
    render_window.SetSize(1000, 500)
    
    iren.SetRenderWindow(render_window)
    iren.AddObserver('KeyPressEvent', keypressCallbackFunction)

    # Create list of files from both directories.
    base_files = glob.glob(sys.argv[1] + '/' + sys.argv[3] + "_Data_t_*.vtu")
    test_files = glob.glob(sys.argv[2] + '/' + sys.argv[3] + "_Data_t_*.vtu")

    # Sorts the files numerically given the time-step.
    sortNicely(base_files)
    sortNicely(test_files)
    
    renderer_0, renderer_1 = vtk.vtkRenderer(),vtk.vtkRenderer()
    render_window.AddRenderer(renderer_0)
    render_window.AddRenderer(renderer_1)
    
    renderer_0.SetViewport(0,0,0.5,1)
    renderer_1.SetViewport(0.5,0,1,1)
    
    update_mappers(base_files[time_step], test_files[time_step])
    
    # Create the mappers that correspond to the objects of the vtk files
    # into graphics elements.
    
    mapper_0.ScalarVisibilityOn()
    mapper_0.SetScalarModeToUseCellData
   
    mapper_1.ScalarVisibilityOn()
    mapper_1.SetScalarModeToUseCellData()

    # Create the Actors
    actor_0, actor_1 = vtkActor(), vtkActor()

    actor_0.SetMapper(mapper_0)   
    actor_1.SetMapper(mapper_1)
    
    renderer_0.AddActor(actor_0)
    renderer_1.AddActor(actor_1)
 
    render_window.Render()
    iren.Start()

if __name__ == "__main__":
    main()
