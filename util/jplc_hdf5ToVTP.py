import sys
import h5py
import vtk
import numpy

"""
Blend JPLC data from HDF5 dataset with the EC geometry from a VTK dataset.

WARNING: The ordering of the cells is to be corrected.
"""

# The idea here is to load the parent and daughter ATP datasets and simply
# replace the ATP values with what we have in the HDF5 file.

# TODO: Sort the correct ordering of the cells.

INPUT_FILES = [
["vtk/ec_mesh_parent.vtp", 'solution/jplc_1.h5'],
["vtk/ec_mesh_left_daughter.vtp", 'solution/jplc_2.h5'],
["vtk/ec_mesh_right_daughter.vtp", 'solution/jplc_3.h5']
]

ecAppend = vtk.vtkAppendPolyData()

for files in INPUT_FILES:
    ecVTKFileReader = vtk.vtkXMLPolyDataReader()
    ecVTKFileReader.SetFileName(files[0])
    ecVTKFileReader.Update()
    ec_mesh = ecVTKFileReader.GetOutput()
    
    if ecVTKFileReader.GetErrorCode() != 0:
        print "Error code", ecVTKFileReader.GetErrorCode()
        print ecVTKFileReader
        sys.exit()

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

jplcWriter = vtk.vtkXMLPolyDataWriter()
jplcWriter.SetFileName('solution/jplc_input.vtp')
jplcWriter.SetInput(outputDataset)
jplcWriter.Update()

