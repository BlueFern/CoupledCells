# This script is a very rough prototype, hence there's a lot of code repetition.

# Read SMC cell data from HDF files and combine them with geometry.

# For each time step there are three HDF5 files, one file per branch.

# The number of time steps to proces is to be specified as a command line argument.

import sys
import h5py
import numpy
import vtk

H5_FILE_BASE_NAME = 'smc_data_t_'
VTP_FILE_BASE_NAME = 'smc_data_t_'

INPUT_SMC_MESH_FILES = ['../vtk/smc_mesh_parent.vtp','../vtk/smc_mesh_left_daughter.vtp', '../vtk/smc_mesh_right_daughter.vtp']

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
    INPUT_SMC_MESHES = []

    # Read input SMC meshes.
    for in_file in INPUT_SMC_MESH_FILES:
        reader = vtk.vtkXMLPolyDataReader()
        reader.SetFileName(in_file)
        reader.Update()

        INPUT_SMC_MESHES += [reader.GetOutput()]

    for time_step in range(int(sys.argv[1])):
        append_filter = vtk.vtkAppendFilter()

        # PARENT.
        mesh_parent = vtk.vtkPolyData()
        mesh_parent.DeepCopy(INPUT_SMC_MESHES[0])

        h5_file_parent = H5_FILE_BASE_NAME + str(time_step) + '_b_1.h5'

        ca_array_parent = read_array(h5_file_parent, '/SMC_Ca')
        ca_array_parent.SetName('SMC_Ca')
        mesh_parent.GetCellData().AddArray(ca_array_parent)

        ca_cpl_array_parent = read_array(h5_file_parent, '/SMC_Ca_coupling')
        ca_cpl_array_parent.SetName('SMC_Ca_coupling')
        mesh_parent.GetCellData().AddArray(ca_cpl_array_parent)

        ip3_array_parent = read_array(h5_file_parent, '/SMC_IP3')
        ip3_array_parent.SetName('SMC_IP3')
        mesh_parent.GetCellData().AddArray(ip3_array_parent)

        ip3_cpl_array_parent = read_array(h5_file_parent, '/SMC_IP3_coupling')
        ip3_cpl_array_parent.SetName('SMC_IP3_coupling')
        mesh_parent.GetCellData().AddArray(ip3_cpl_array_parent)

        vm_array_parent = read_array(h5_file_parent, '/SMC_Vm')
        vm_array_parent.SetName('SMC_Vm')
        mesh_parent.GetCellData().AddArray(vm_array_parent)

        vm_cpl_array_parent = read_array(h5_file_parent, '/SMC_Vm_coupling')
        vm_cpl_array_parent.SetName('SMC_Vm_coupling')
        mesh_parent.GetCellData().AddArray(vm_cpl_array_parent)

        sr_array_parent = read_array(h5_file_parent, '/SMC_SR')
        sr_array_parent.SetName('SMC_SR')
        mesh_parent.GetCellData().AddArray(sr_array_parent)

        w_array_parent = read_array(h5_file_parent, '/SMC_w')
        w_array_parent.SetName('SMC_w')
        mesh_parent.GetCellData().AddArray(w_array_parent)

        # LEFT.
        mesh_left = vtk.vtkPolyData()
        mesh_left.DeepCopy(INPUT_SMC_MESHES[1])

        h5_file_left = H5_FILE_BASE_NAME + str(time_step) + '_b_2.h5'

        ca_array_left = read_array(h5_file_left,'/SMC_Ca')
        ca_array_left.SetName('SMC_Ca')
        mesh_left.GetCellData().AddArray(ca_array_left)

        ca_cpl_array_left = read_array(h5_file_left, '/SMC_Ca_coupling')
        ca_cpl_array_left.SetName('SMC_Ca_coupling')
        mesh_left.GetCellData().AddArray(ca_cpl_array_left)

        ip3_array_left = read_array(h5_file_left, '/SMC_IP3')
        ip3_array_left.SetName('SMC_IP3')
        mesh_left.GetCellData().AddArray(ip3_array_left)

        ip3_cpl_array_left = read_array(h5_file_left, '/SMC_IP3_coupling')
        ip3_cpl_array_left.SetName('SMC_IP3_coupling')
        mesh_left.GetCellData().AddArray(ip3_cpl_array_left)

        vm_array_left = read_array(h5_file_left, '/SMC_Vm')
        vm_array_left.SetName('SMC_Vm')
        mesh_left.GetCellData().AddArray(vm_array_left)

        vm_cpl_array_left = read_array(h5_file_left, '/SMC_Vm_coupling')
        vm_cpl_array_left.SetName('SMC_Vm_coupling')
        mesh_left.GetCellData().AddArray(vm_cpl_array_left)

        sr_array_left = read_array(h5_file_left, '/SMC_SR')
        sr_array_left.SetName('SMC_SR')
        mesh_left.GetCellData().AddArray(sr_array_left)

        w_array_left = read_array(h5_file_left, '/SMC_w')
        w_array_left.SetName('SMC_w')
        mesh_left.GetCellData().AddArray(w_array_left)

        # RIGHT.
        mesh_right = vtk.vtkPolyData()
        mesh_right.DeepCopy(INPUT_SMC_MESHES[2])

        h5_file_right = H5_FILE_BASE_NAME + str(time_step) + '_b_3.h5'

        ca_array_right = read_array(h5_file_right, '/SMC_Ca')
        ca_array_right.SetName('SMC_Ca')
        mesh_right.GetCellData().AddArray(ca_array_right)

        ca_cpl_array_right = read_array(h5_file_right, '/SMC_Ca_coupling')
        ca_cpl_array_right.SetName('SMC_Ca_coupling')
        mesh_right.GetCellData().AddArray(ca_cpl_array_right)

        ip3_array_right = read_array(h5_file_right, '/SMC_IP3')
        ip3_array_right.SetName('SMC_IP3')
        mesh_right.GetCellData().AddArray(ip3_array_right)

        ip3_cpl_array_right = read_array(h5_file_right, '/SMC_IP3_coupling')
        ip3_cpl_array_right.SetName('SMC_IP3_coupling')
        mesh_right.GetCellData().AddArray(ip3_cpl_array_right)

        vm_array_right = read_array(h5_file_right, '/SMC_Vm')
        vm_array_right.SetName('SMC_Vm')
        mesh_right.GetCellData().AddArray(vm_array_right)

        vm_cpl_array_right = read_array(h5_file_right, '/SMC_Vm_coupling')
        vm_cpl_array_right.SetName('SMC_Vm_coupling')
        mesh_right.GetCellData().AddArray(vm_cpl_array_right)

        sr_array_right = read_array(h5_file_right, '/SMC_SR')
        sr_array_right.SetName('SMC_SR')
        mesh_right.GetCellData().AddArray(sr_array_right)

        w_array_right = read_array(h5_file_right, '/SMC_w')
        w_array_right.SetName('SMC_w')
        mesh_right.GetCellData().AddArray(w_array_right)

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

