import sys
import h5py
import numpy
import vtk

"""
This script is a very rough prototype, hence there's a lot of code repetition.

Read EC cell data from HDF files and combine them with geometry.

For each time step there are three HDF5 files, one file per branch.

The number of time steps to proces is to be specified as a command line argument.
"""

H5_FILE_BASE_NAME = 'solution/ec_data_t_'
VTP_FILE_BASE_NAME = 'solution/ec_data_t_'

INPUT_EC_MESH_FILES = [
'vtk/ec_mesh_parent.vtp',
]

attributes = ['EC_Ca', 'EC_Ca_coupling', 'EC_IP3', 'EC_IP3_coupling', 'EC_SR', 'EC_Vm', 'EC_Vm_coupling']

    
def read_array(arrays, h5_file_name, dataset_name):
    fid = h5py.h5f.open(h5_file_name)
    
    dset = h5py.h5d.open(fid, dataset_name)
    shape = dset.shape
    rdata = numpy.zeros((shape[0], shape[1]), dtype=numpy.float64)
    dset.read(h5py.h5s.ALL, h5py.h5s.ALL, rdata)

    arr = rdata.ravel()

    i = 0
    for val in arr:
        arrays[i % len(attributes)].InsertNextValue(val)
        i += 1

def run():
    INPUT_EC_MESHES = []

    # Read input EC meshes.
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(INPUT_EC_MESH_FILES[0])
    reader.Update()

    INPUT_EC_MESHES += [reader.GetOutput()]

    for time_step in range(int(sys.argv[1]), int(sys.argv[2])):
        
        array = []
        for att in attributes:
            array.append(vtk.vtkDoubleArray())
            array[-1].SetName(att)
                
        append_filter = vtk.vtkAppendFilter()

        # PARENT.
        mesh_parent = vtk.vtkPolyData()
        mesh_parent.DeepCopy(INPUT_EC_MESHES[0])

        h5_file_parent_base = H5_FILE_BASE_NAME + str(time_step) + '_b_1_' + 'x' + '.h5'
        print "Processing file", h5_file_parent_base


        for i in range(int(sys.argv[3])):
            h5_file_parent = h5_file_parent_base[:-4] + str(i) + h5_file_parent_base[-3:]
            
            read_array(array, h5_file_parent, "data")
                
        for i in range(len(attributes)):
            mesh_parent.GetCellData().AddArray(array[i])

        # Append parent.
        append_filter.AddInput(mesh_parent)
        append_filter.Update()

        # Write the result.
        vtp_file = VTP_FILE_BASE_NAME + str(time_step) + '.vtu'
        writer = vtk.vtkXMLUnstructuredGridWriter()
        writer.SetFileName(vtp_file)
        writer.SetInput(append_filter.GetOutput())
        writer.Update()

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print "Expected arguments: the number of time steps to process and the number of writers"
    else:
        run()

