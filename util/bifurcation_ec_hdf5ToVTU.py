import os
import h5py
import numpy
import argparse
import vtk

"""
This script is a very rough prototype, hence there's a lot of code repetition.

Read EC cell data from HDF files and combine them with geometry.

For each time step there are three HDF5 files, one file per branch.

The number of time steps to proces is to be specified as a command line argument.
"""

H5_FILE_BASE_NAME = 'solution/ec_data_t_'
VTU_FILE_BASE_NAME = 'solution/ec_data_t_'

INPUT_EC_MESH_FILES = [
    'vtk/ec_mesh_parent.vtp',
    'vtk/ec_mesh_left_daughter.vtp',
    'vtk/ec_mesh_right_daughter.vtp'
]

numEcsAxially = 4
numEcsCirc = 20

numSmcsAxially = 52
numSmcsCirc = 4

circScale = 0
axialScale = 0

circQuads = 0
axialQuads = 0


def read_array(h5_file_name, dataset_name):
    fid = h5py.h5f.open(h5_file_name)
    dset = h5py.h5d.open(fid, dataset_name)
    shape = dset.shape
    rdata = numpy.zeros((shape[0], shape[1]), dtype=numpy.float64)
    dset.read(h5py.h5s.ALL, h5py.h5s.ALL, rdata)

    arr = rdata.ravel()
    
    array = vtk.vtkDoubleArray()
    array.SetNumberOfValues(numEcsAxially * numEcsCirc * circQuads * axialQuads)
       
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
            
            for j in range(0, numEcsAxially): #Ec rows
                rowOffset = j * numEcsCirc
                
                for extraCirc in range(0, circScale):
                    circOffset = extraCirc * numEcsAxially * numEcsCirc               
                    
                    for k in range(0, numEcsCirc): #Ec cols
                        
                        cellId = quadOffset + rowOffset + circOffset + axialOffset + k
                        array.SetValue(cellId,arr[val])
                        val += 1

    return array


def HDF5toVTK(start, end):
    INPUT_EC_MESHES = []

    # Read input EC meshes.
    for in_file in INPUT_EC_MESH_FILES:
        reader = vtk.vtkXMLPolyDataReader()
        reader.SetFileName(in_file)
        reader.Update()

        INPUT_EC_MESHES += [reader.GetOutput()]

    for time_step in range(start, end + 1):
        append_filter = vtk.vtkAppendFilter()

        # PARENT.
        mesh_parent = vtk.vtkPolyData()
        mesh_parent.DeepCopy(INPUT_EC_MESHES[0])

        h5_file_parent = H5_FILE_BASE_NAME + str(time_step) + '_b_1.h5'
        print "Processing file", h5_file_parent

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
        print "Processing file", h5_file_left

        ca_array_left = read_array(h5_file_left, '/EC_Ca')
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
        print "Processing file", h5_file_right

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
        if vtk.VTK_MAJOR_VERSION < 6:
            append_filter.AddInput(mesh_parent)
            append_filter.AddInput(mesh_left)
            append_filter.AddInput(mesh_right)
        else:
            append_filter.AddInputData(mesh_parent)
            append_filter.AddInputData(mesh_left)
            append_filter.AddInputData(mesh_right)
        append_filter.Update()

        # Write the result.
        vtu_file = VTU_FILE_BASE_NAME + str(time_step) + '.vtu'
        print 'Writing file', os.path.abspath(vtu_file)

        writer = vtk.vtkXMLUnstructuredGridWriter()
        writer.SetFileName(vtu_file)
        if vtk.VTK_MAJOR_VERSION < 6:
            writer.SetInput(append_filter.GetOutput())
        else:
            writer.SetInputData(append_filter.GetOutput())
        writer.Update()


if __name__ == "__main__":   
    argParser = argparse.ArgumentParser(
        description='Fuse Coupled Cells simulation (bifurcation) data in HDF5 format with vessel geomtery data in VTK format and save the result as VTU data.')
    argParser.add_argument('start', type=int, help='Start time')
    argParser.add_argument('end', type=int, help='End time')
    
    argParser.add_argument('circQuads', type=int, help='Number of circumferential quads for each branch.')
    argParser.add_argument('axialQuads', type=int, help='Number of axial quads for each branch.')
    
    argParser.add_argument('circScale', type=int, help='Number of circumferential quads joined.')
    argParser.add_argument('axialScale', type=int, help='Number of axial quads joined.')
    args = argParser.parse_args()
        
    circScale = args.circScale
    axialScale = args.axialScale

    circQuads = args.circQuads
    axialQuads = args.axialQuads

    HDF5toVTK(args.start, args.end)
