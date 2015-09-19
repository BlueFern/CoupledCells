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
'vtk/ec_mesh_left_daughter.vtp',
'vtk/ec_mesh_right_daughter.vtp'
]

def read_array(h5_file_name, dataset_name):
    fid = h5py.h5f.open(h5_file_name)
    dset = h5py.h5d.open(fid, dataset_name)
    shape = dset.shape
    rdata = numpy.zeros((shape[0], shape[1]), dtype = numpy.float64)
    dset.read(h5py.h5s.ALL, h5py.h5s.ALL, rdata)

    arr = rdata.ravel()

    array = vtk.vtkDoubleArray()

    for val in arr:
        array.InsertNextValue(val)

    return array

def run():
    INPUT_EC_MESHES = []

    # Read input EC meshes.
    for in_file in INPUT_EC_MESH_FILES:
        reader = vtk.vtkXMLPolyDataReader()
        reader.SetFileName(in_file)
        reader.Update()

        INPUT_EC_MESHES += [reader.GetOutput()]

    for time_step in range(int(sys.argv[1])):
        append_filter = vtk.vtkAppendFilter()

        # PARENT.
        mesh_parent = vtk.vtkPolyData()
        mesh_parent.DeepCopy(INPUT_EC_MESHES[0])

        h5_file_parent = H5_FILE_BASE_NAME + str(time_step) + '_b_1.h5'

        ca_array_parent = read_array(h5_file_parent, '/EC_Ca')
        ca_array_parent.SetName('EC_Ca')
        mesh_parent.GetCellData().AddArray(ca_array_parent)

        ca_cpl_array_parent = read_array(h5_file_parent, '/EC_Ca_coupling')
        ca_cpl_array_parent.SetName('EC_Ca_coupling')
        mesh_parent.GetCellData().AddArray(ca_cpl_array_parent)

        ip3_array_parent = read_array(h5_file_parent, '/EC_IP3')
        ip3_array_parent.SetName('EC_IP3')
        mesh_parent.GetCellData().AddArray(ip3_array_parent)

        ip3_cpl_array_parent = read_array(h5_file_parent, '/EC_IP3_coupling')
        ip3_cpl_array_parent.SetName('EC_IP3_coupling')
        mesh_parent.GetCellData().AddArray(ip3_cpl_array_parent)

        vm_array_parent = read_array(h5_file_parent, '/EC_Vm')
        vm_array_parent.SetName('EC_Vm')
        mesh_parent.GetCellData().AddArray(vm_array_parent)

        vm_cpl_array_parent = read_array(h5_file_parent, '/EC_Vm_coupling')
        vm_cpl_array_parent.SetName('EC_Vm_coupling')
        mesh_parent.GetCellData().AddArray(vm_cpl_array_parent)

        sr_array_parent = read_array(h5_file_parent, '/EC_SR')
        sr_array_parent.SetName('EC_SR')
        mesh_parent.GetCellData().AddArray(sr_array_parent)

        # LEFT.
        mesh_left = vtk.vtkPolyData()
        mesh_left.DeepCopy(INPUT_EC_MESHES[1])

        h5_file_left = H5_FILE_BASE_NAME + str(time_step) + '_b_2.h5'

        ca_array_left = read_array(h5_file_left,'/EC_Ca')
        ca_array_left.SetName('EC_Ca')
        mesh_left.GetCellData().AddArray(ca_array_left)

        ca_cpl_array_left = read_array(h5_file_left, '/EC_Ca_coupling')
        ca_cpl_array_left.SetName('EC_Ca_coupling')
        mesh_left.GetCellData().AddArray(ca_cpl_array_left)

        ip3_array_left = read_array(h5_file_left, '/EC_IP3')
        ip3_array_left.SetName('EC_IP3')
        mesh_left.GetCellData().AddArray(ip3_array_left)

        ip3_cpl_array_left = read_array(h5_file_left, '/EC_IP3_coupling')
        ip3_cpl_array_left.SetName('EC_IP3_coupling')
        mesh_left.GetCellData().AddArray(ip3_cpl_array_left)

        vm_array_left = read_array(h5_file_left, '/EC_Vm')
        vm_array_left.SetName('EC_Vm')
        mesh_left.GetCellData().AddArray(vm_array_left)

        vm_cpl_array_left = read_array(h5_file_left, '/EC_Vm_coupling')
        vm_cpl_array_left.SetName('EC_Vm_coupling')
        mesh_left.GetCellData().AddArray(vm_cpl_array_left)

        sr_array_left = read_array(h5_file_left, '/EC_SR')
        sr_array_left.SetName('EC_SR')
        mesh_left.GetCellData().AddArray(sr_array_left)

        # RIGHT.
        mesh_right = vtk.vtkPolyData()
        mesh_right.DeepCopy(INPUT_EC_MESHES[2])

        h5_file_right = H5_FILE_BASE_NAME + str(time_step) + '_b_3.h5'

        ca_array_right = read_array(h5_file_right, '/EC_Ca')
        ca_array_right.SetName('EC_Ca')
        mesh_right.GetCellData().AddArray(ca_array_right)

        ca_cpl_array_right = read_array(h5_file_right, '/EC_Ca_coupling')
        ca_cpl_array_right.SetName('EC_Ca_coupling')
        mesh_right.GetCellData().AddArray(ca_cpl_array_right)

        ip3_array_right = read_array(h5_file_right, '/EC_IP3')
        ip3_array_right.SetName('EC_IP3')
        mesh_right.GetCellData().AddArray(ip3_array_right)

        ip3_cpl_array_right = read_array(h5_file_right, '/EC_IP3_coupling')
        ip3_cpl_array_right.SetName('EC_IP3_coupling')
        mesh_right.GetCellData().AddArray(ip3_cpl_array_right)

        vm_array_right = read_array(h5_file_right, '/EC_Vm')
        vm_array_right.SetName('EC_Vm')
        mesh_right.GetCellData().AddArray(vm_array_right)

        vm_cpl_array_right = read_array(h5_file_right, '/EC_Vm_coupling')
        vm_cpl_array_right.SetName('EC_Vm_coupling')
        mesh_right.GetCellData().AddArray(vm_cpl_array_right)

        sr_array_right = read_array(h5_file_right, '/EC_SR')
        sr_array_right.SetName('EC_SR')
        mesh_right.GetCellData().AddArray(sr_array_right)

        # Append parent, left, right.
        append_filter.AddInput(mesh_parent)
        append_filter.AddInput(mesh_left)
        append_filter.AddInput(mesh_right)
        append_filter.Update()

        # Write the result.
        vtp_file = VTP_FILE_BASE_NAME + str(time_step) + '.vtu'
        writer = vtk.vtkXMLUnstructuredGridWriter()
        writer.SetFileName(vtp_file)
        writer.SetInput(append_filter.GetOutput())
        writer.Update()

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print "Expected arguments: the number of time steps to process."
    else:
        run()

