# This script is a very rough prototype, hence there's a lot of code repetition.

# Read EC cell data from HDF files and combine them with geometry.

# For each time step there are three HDF5 files, one file per branch.

# The number of time steps to proces is to be specified as a command line argument.

import sys
import h5py
import numpy
import vtk

H5_FILE_BASE_NAME = 'ec_data_t_'
VTP_FILE_BASE_NAME = 'ec_data_t_'

INPUT_EC_MESH_FILES = ['../vtk/ec_mesh_parent.vtp','../vtk/ec_mesh_left_daughter.vtp', '../vtk/ec_mesh_right_daughter.vtp']

def read_array(h5_file_name, dataset_name):
    fid = h5py.h5f.open(h5_file_name)
    dset = h5py.h5d.open(fid, dataset_name)
    rdata = numpy.zeros((8, 80), dtype = numpy.float64)
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
        append_filter = vtk.vtkAppendPolyData()

        # PARENT.
        mesh_parent = vtk.vtkPolyData()
        mesh_parent.DeepCopy(INPUT_EC_MESHES[0])

        h5_file_parent = H5_FILE_BASE_NAME + str(time_step) + '_b_1.h5'

        ca_array_parent = read_array(h5_file_parent, '/ec_Ca')
        ca_array_parent.SetName('Ca')
        mesh_parent.GetCellData().AddArray(ca_array_parent)

        ca_cpl_array_parent = read_array(h5_file_parent, '/ec_cpl_Ca')
        ca_cpl_array_parent.SetName('cpl_Ca')
        mesh_parent.GetCellData().AddArray(ca_cpl_array_parent)

        ip3_array_parent = read_array(h5_file_parent, '/ec_IP3')
        ip3_array_parent.SetName('IP3')
        mesh_parent.GetCellData().AddArray(ip3_array_parent)

        ip3_cpl_array_parent = read_array(h5_file_parent, '/ec_cpl_IP3')
        ip3_cpl_array_parent.SetName('cpl_IP3')
        mesh_parent.GetCellData().AddArray(ip3_cpl_array_parent)

        vm_array_parent = read_array(h5_file_parent, '/ec_Vm')
        vm_array_parent.SetName('Vm')
        mesh_parent.GetCellData().AddArray(vm_array_parent)

        vm_cpl_array_parent = read_array(h5_file_parent, '/ec_cpl_Vm')
        vm_cpl_array_parent.SetName('cpl_Vm')
        mesh_parent.GetCellData().AddArray(vm_cpl_array_parent)


        sr_array_parent = read_array(h5_file_parent, '/ec_SR')
        sr_array_parent.SetName('Vm')
        mesh_parent.GetCellData().AddArray(sr_array_parent)

        # LEFT.
        mesh_left = vtk.vtkPolyData()
        mesh_left.DeepCopy(INPUT_EC_MESHES[1])

        h5_file_left = H5_FILE_BASE_NAME + str(time_step) + '_b_2.h5'

        ca_array_left = read_array(h5_file_left,'/ec_Ca')
        ca_array_left.SetName('Ca')
        mesh_left.GetCellData().AddArray(ca_array_left)

        ca_cpl_array_left = read_array(h5_file_left, '/ec_cpl_Ca')
        ca_cpl_array_left.SetName('cpl_Ca')
        mesh_left.GetCellData().AddArray(ca_cpl_array_left)

        ip3_array_left = read_array(h5_file_left, '/ec_IP3')
        ip3_array_left.SetName('IP3')
        mesh_left.GetCellData().AddArray(ip3_array_left)

        ip3_cpl_array_left = read_array(h5_file_left, '/ec_cpl_IP3')
        ip3_cpl_array_left.SetName('cpl_IP3')
        mesh_left.GetCellData().AddArray(ip3_cpl_array_left)

        vm_array_left = read_array(h5_file_left, '/ec_Vm')
        vm_array_left.SetName('Vm')
        mesh_left.GetCellData().AddArray(vm_array_left)

        vm_cpl_array_left = read_array(h5_file_left, '/ec_cpl_Vm')
        vm_cpl_array_left.SetName('cpl_Vm')
        mesh_left.GetCellData().AddArray(vm_cpl_array_left)


        sr_array_left = read_array(h5_file_left, '/ec_SR')
        sr_array_left.SetName('Vm')
        mesh_left.GetCellData().AddArray(sr_array_left)

        # RIGHT.
        mesh_right = vtk.vtkPolyData()
        mesh_right.DeepCopy(INPUT_EC_MESHES[2])

        h5_file_right = H5_FILE_BASE_NAME + str(time_step) + '_b_3.h5'

        ca_array_right = read_array(h5_file_right, '/ec_Ca')
        ca_array_right.SetName('Ca')
        mesh_right.GetCellData().AddArray(ca_array_right)

        ca_cpl_array_right = read_array(h5_file_right, '/ec_cpl_Ca')
        ca_cpl_array_right.SetName('cpl_Ca')
        mesh_right.GetCellData().AddArray(ca_cpl_array_right)

        ip3_array_right = read_array(h5_file_right, '/ec_IP3')
        ip3_array_right.SetName('IP3')
        mesh_right.GetCellData().AddArray(ip3_array_right)

        ip3_cpl_array_right = read_array(h5_file_right, '/ec_cpl_IP3')
        ip3_cpl_array_right.SetName('cpl_IP3')
        mesh_right.GetCellData().AddArray(ip3_cpl_array_right)

        vm_array_right = read_array(h5_file_right, '/ec_Vm')
        vm_array_right.SetName('Vm')
        mesh_right.GetCellData().AddArray(vm_array_right)

        vm_cpl_array_right = read_array(h5_file_right, '/ec_cpl_Vm')
        vm_cpl_array_right.SetName('cpl_Vm')
        mesh_right.GetCellData().AddArray(vm_cpl_array_right)


        sr_array_right = read_array(h5_file_right, '/ec_SR')
        sr_array_right.SetName('Vm')
        mesh_right.GetCellData().AddArray(sr_array_right)

	# Append parent, left, right.
        append_filter.AddInput(mesh_parent)
        append_filter.AddInput(mesh_left)
        append_filter.AddInput(mesh_right)
        append_filter.Update()

        # Write the result.
	vtp_file = VTP_FILE_BASE_NAME + str(time_step) + '.vtp'
        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName(vtp_file)
        writer.SetInput(append_filter.GetOutput())
        writer.Update()

if __name__ == "__main__":
    if sys.argc != 2:
        "The number of time steps to process is to be provided as the only command line argument."
        return

    run()

