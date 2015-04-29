import sys
import h5py
import vtk

"""
Blend JPLC data from HDF5 dataset with the EC geometry from a VTK dataset.

WARNING: The ordering of the cells is to be corrected.
"""

# The idea here is to load the parent and daughter ATP datasets and simply
# replace the ATP values with what we have in the HDF5 file.

# TODO: Sort the correct ordering of the cells.

ecVTKFiles = [
"vtk/ec_mesh_parent.vtp",
"vtk/ec_mesh_left_daughter.vtp",
"vtk/ec_mesh_right_daughter.vtp",
]

ecAppend = vtk.vtkAppendPolyData()

for ecVTKFile in ecVTKFiles:
    ecVTKFileReader = vtk.vtkXMLPolyDataReader()
    ecVTKFileReader.SetFileName(ecVTKFile)
    ecVTKFileReader.Update()
    
    if ecVTKFileReader.GetErrorCode() != 0:
        print "Error code", ecVTKFileReader.GetErrorCode()
        print ecVTKFileReader
        sys.exit()

    print ecVTKFileReader.GetOutput().GetNumberOfCells()
    ecAppend.AddInput(ecVTKFileReader.GetOutput())

ecAppend.Update()
outputDataset = ecAppend.GetOutput()

jplcFile = h5py.File('solution/jplc.h5', 'r')
print jplcFile

jplcData = jplcFile['jplc']

rows, cols = jplcData.shape

jplcArray = vtk.vtkDoubleArray()
jplcArray.SetName('JPLC')

print rows, cols
for r in range(rows):
   for c in range(cols):
       jplcArray.InsertNextValue(jplcData[r][c])

outputDataset.GetCellData().SetScalars(jplcArray)

jplcWriter = vtk.vtkXMLPolyDataWriter()
jplcWriter.SetFileName('solution/Agonist_map.vtp')
jplcWriter.SetInput(outputDataset)
jplcWriter.Update()