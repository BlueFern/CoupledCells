import sys
import h5py
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
    if vtk.VTK_MAJOR_VERSION < 6:
        ecAppend.AddInput(ec_mesh)
    else:
        ecAppend.AddInputData(ec_mesh)

ecAppend.Update()
outputDataset = ecAppend.GetOutput()

jplcWriter = vtk.vtkXMLUnstructuredGridWriter()
jplcWriter.SetFileName('solution/jplc_input.vtu')
if vtk.VTK_MAJOR_VERSION < 6:
    jplcWriter.SetInput(outputDataset)
else:
    jplcWriter.SetInputData(outputDataset)
jplcWriter.Update()

