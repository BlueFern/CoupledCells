import h5py
import numpy
import vtk

FILE_VTP = 'vtk/ec_mesh_parent.vtp'

reader = vtk.vtkXMLPolyDataReader()
reader.SetFileName(FILE_VTP)
reader.Update()

ec_mesh = reader.GetOutput()

FILE_HDF5 = 'solution/jplc_1.h5'

fid = h5py.h5f.open(FILE_HDF5)
dset = h5py.h5d.open(fid, '/jplc')
rdata = numpy.zeros((8, 80), dtype = numpy.float64)
dset.read(h5py.h5s.ALL, h5py.h5s.ALL, rdata)

arr = rdata.ravel()

jplcArray = vtk.vtkDoubleArray()
jplcArray.SetName('JPLC')

for val in arr:
	jplcArray.InsertNextValue(val)

ec_mesh.GetCellData().SetScalars(jplcArray)

writer = vtk.vtkXMLPolyDataWriter()
writer.SetFileName('agonist_parent.vtp')
writer.SetInput(ec_mesh)
writer.Update()



