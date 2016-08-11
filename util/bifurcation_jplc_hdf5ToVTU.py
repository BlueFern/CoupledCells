import os
import sys
import h5py
import argparse
import vtk
import numpy

# ['VTK_BUILD_VERSION', 'VTK_MAJOR_VERSION', 'VTK_MINOR_VERSION', 'VTK_SOURCE_VERSION', 'VTK_VERSION', 'vtkFastNumericConversion', 'vtkVersion']

"""
Blend JPLC data from HDF5 dataset with the EC geometry from a VTK dataset.

The idea here is to load the parent and daughter ATP datasets and simply
replace the ATP values with what we have in the HDF5 file.
"""

INPUT_FILES = [
["vtk/ec_mesh_parent.vtp", 'solution/jplc_1.h5'],
["vtk/ec_mesh_left_daughter.vtp", 'solution/jplc_2.h5'],
["vtk/ec_mesh_right_daughter.vtp", 'solution/jplc_3.h5']
]

numEcsAxially = 4
numEcsCirc = 20

circScale = 0
axialScale = 0

circQuads = 0
axialQuads = 0

def HDF5toVTK():
    ecAppend = vtk.vtkAppendFilter()
    
    for files in INPUT_FILES:
        ecVTKFileReader = vtk.vtkXMLPolyDataReader()
        ecVTKFileReader.SetFileName(files[0])
        ecVTKFileReader.Update()
        ec_mesh = ecVTKFileReader.GetOutput()
        
        if ecVTKFileReader.GetErrorCode() != 0:
            print "Error code", ecVTKFileReader.GetErrorCode()
            print ecVTKFileReader
            sys.exit()
    
        print "Processing file", files[1]
    
        fid = h5py.h5f.open(files[1])
        dset = h5py.h5d.open(fid, '/jplc')
        shape = dset.shape
        rdata = numpy.zeros((shape[0], shape[1]), dtype=numpy.float64)
        dset.read(h5py.h5s.ALL, h5py.h5s.ALL, rdata)
        arr = rdata.ravel()
    
        jplcArray = vtk.vtkDoubleArray()
        jplcArray.SetNumberOfValues(numEcsAxially * numEcsCirc * circQuads * axialQuads)
        jplcArray.SetName('JPLC')
           
        requiredQuads = []      
        
        for axialquad in range(axialQuads):
            for circQuad in range(circQuads):
                if circQuad % circScale == 0 and axialquad % axialScale == 0:
                    requiredQuads.append(circQuad + (axialquad * circQuads))
           
    
        val = 0
        for i in requiredQuads: 
            quadOffset = i * numEcsAxially * numEcsCirc
        
            for extraAxial in range(0, axialScale):
                axialOffset = extraAxial * numEcsCirc * numEcsAxially * circQuads
    	    
                for j in range(0, numEcsAxially): #smc rows
        	         rowOffset = j * numEcsCirc
        	        
        	         for extraCirc in range(0, circScale):
        	             circOffset = extraCirc * numEcsAxially * numEcsCirc               
        	            
        	             for k in range(0, numEcsCirc): #smc cols
        	                
        	                 cellId = quadOffset + rowOffset + circOffset + axialOffset + k
        	                 jplcArray.SetValue(cellId,arr[val])
        	                 val += 1

    
        ec_mesh.GetCellData().SetScalars(jplcArray)
        if vtk.VTK_MAJOR_VERSION < 6:
            ecAppend.AddInput(ec_mesh)
        else:
            ecAppend.AddInputData(ec_mesh)
    
    ecAppend.Update()
    outputDataset = ecAppend.GetOutput()
    
    outputFileName = 'solution/jplc_input.vtu'
    print 'Writing file', os.path.abspath(outputFileName)
    
    jplcWriter = vtk.vtkXMLUnstructuredGridWriter()
    jplcWriter.SetFileName(outputFileName)
    if vtk.VTK_MAJOR_VERSION < 6:
        jplcWriter.SetInput(outputDataset)
    else:
        jplcWriter.SetInputData(outputDataset)
    jplcWriter.Update()

if __name__ == "__main__":
    argParser = argparse.ArgumentParser()  
    argParser.add_argument('circQuads', type=int, help='Number of circumferential quads for each branch.')
    argParser.add_argument('axialQuads', type=int, help='Number of axial quads for each branch.')
    
    argParser.add_argument('circScale', type=int, help='Number of circumferential quads joined.')
    argParser.add_argument('axialScale', type=int, help='Number of axial quads joined.')
    args = argParser.parse_args()
        
    circScale = args.circScale
    axialScale = args.axialScale
    circQuads = args.circQuads
    axialQuads = args.axialQuads

    HDF5toVTK()

