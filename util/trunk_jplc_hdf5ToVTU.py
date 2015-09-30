import sys
import h5py
import vtk
import numpy

"""
Blend JPLC data from HDF5 dataset with the EC geometry from a VTK dataset.

The idea here is to load the parent and daughter ATP datasets and simply
replace the ATP values with what we have in the HDF5 file.
"""

INPUT_FILES = [
["vtk/ec_mesh_parent.vtp", 'solution/jplc_0.h5']
]

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
    rdata = numpy.zeros((shape[0], shape[1]), dtype = numpy.float64)
    dset.read(h5py.h5s.ALL, h5py.h5s.ALL, rdata)

    arr = rdata.ravel()

    jplcArray = vtk.vtkDoubleArray()
    jplcArray.SetName('JPLC')

    for val in arr:
        jplcArray.InsertNextValue(val)

    ec_mesh.GetCellData().SetScalars(jplcArray)
    ecAppend.AddInput(ec_mesh)

ecAppend.Update()
outputDataset = ecAppend.GetOutput()

jplcWriter = vtk.vtkXMLUnstructuredGridWriter()
jplcWriter.SetFileName('solution/jplc_input.vtu')
jplcWriter.SetInput(outputDataset)
jplcWriter.Update()