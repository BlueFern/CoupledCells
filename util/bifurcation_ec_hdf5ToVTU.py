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


attributes = ['EC_Ca', 'EC_Ca_coupling', 'EC_IP3', 'EC_IP3_coupling', 'EC_SR', 'EC_Vm', 'EC_Vm_coupling']


def read_array(arrays, h5_file_name, dataset_name, branch):
    fid = h5py.h5f.open(h5_file_name)
    
    dset = h5py.h5d.open(fid, dataset_name)
    shape = dset.shape
    rdata = numpy.zeros((shape[0], shape[1]), dtype=numpy.float64)
    dset.read(h5py.h5s.ALL, h5py.h5s.ALL, rdata)

    arr = rdata.ravel()

    i = 0
    for val in arr:
        arrays[branch][i % len(attributes)].InsertNextValue(val)
        i += 1


def HDF5toVTK(start, end, writers):

    INPUT_EC_MESHES = []

    # Read input EC meshes.
    for in_file in INPUT_EC_MESH_FILES:
        reader = vtk.vtkXMLPolyDataReader()
        reader.SetFileName(in_file)
        reader.Update()

        INPUT_EC_MESHES += [reader.GetOutput()]

    for time_step in range(start, end + 1):
        
        arrays = [[],[],[]]
        for b in range(3):
            for att in attributes:
                arrays[b].append(vtk.vtkDoubleArray())
                arrays[b][-1].SetName(att)
            
        append_filter = vtk.vtkAppendFilter()

        # PARENT.
        branch = 0
        mesh_parent = vtk.vtkPolyData()
        mesh_parent.DeepCopy(INPUT_EC_MESHES[0])
        

        h5_file_parent_base = H5_FILE_BASE_NAME + str(time_step) + '_b_1_' + 'x' + '.h5'
        
        print "Processing file", h5_file_parent_base
        
        
        for i in range(writers):
            h5_file_parent = h5_file_parent_base[:-4] + str(i) + h5_file_parent_base[-3:]
            
            read_array(arrays, h5_file_parent,"data", branch)
                
        for i in range(len(attributes)):
            mesh_parent.GetCellData().AddArray(arrays[branch][i])
        
        
         # LEFT.
        branch = 1
        mesh_left = vtk.vtkPolyData()
        mesh_left.DeepCopy(INPUT_EC_MESHES[1])

        h5_file_left_base = H5_FILE_BASE_NAME + str(time_step) + '_b_2_' + 'x' + '.h5'
        print "Processing file", h5_file_left_base
        
        for i in range(writers):
            h5_file_left = h5_file_left_base[:-4] + str(i) + h5_file_left_base[-3:]
            
            read_array(arrays, h5_file_left,"data", branch)
            
        for i in range(len(attributes)):
            mesh_left.GetCellData().AddArray(arrays[branch][i])
            
         # RIGHT.
        branch = 2
        mesh_right = vtk.vtkPolyData()
        mesh_right.DeepCopy(INPUT_EC_MESHES[2])

        h5_file_right_base = H5_FILE_BASE_NAME + str(time_step) + '_b_3_' + 'x' + '.h5'
        print "Processing file", h5_file_right_base
        
        for i in range(writers):
            h5_file_right = h5_file_right_base[:-4] + str(i) + h5_file_right_base[-3:]
            
            read_array(arrays, h5_file_right, "data", branch)
                
        for i in range(len(attributes)):
            mesh_right.GetCellData().AddArray(arrays[branch][i])
        
    
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
    argParser.add_argument('writers', type=int, help='Number of writers per branch')

    args = argParser.parse_args()    

    HDF5toVTK(args.start, args.end, args.writers)
