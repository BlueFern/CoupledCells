import os
import h5py
import numpy
import argparse
import vtk
import glob

numEcsAxially = 4
numEcsCirc = 20

numSmcsAxially = 52
numSmcsCirc = 4

numCells = {"ec" : [numEcsAxially, numEcsCirc], "smc" : [numSmcsAxially, numSmcsCirc]}

smc_attributes = ['SMC_Ca', 'SMC_Ca_coupling', 'SMC_IP3', 'SMC_IP3_coupling', 'SMC_Vm', 'SMC_Vm_coupling', 'SMC_SR', 'SMC_w']
ec_attributes = ['EC_Ca', 'EC_Ca_coupling', 'EC_IP3', 'EC_IP3_coupling', 'EC_SR', 'EC_Vm', 'EC_Vm_coupling', 'EC_Gprot']
attributes = {"ec" : ec_attributes, "smc" : smc_attributes}

ec_mesh_files = [
    'vtk/ec_mesh_parent.vtp',
    'vtk/ec_mesh_left_daughter.vtp',
    'vtk/ec_mesh_right_daughter.vtp'
]
smc_mesh_files = [
    'vtk/smc_mesh_parent.vtp',
    'vtk/smc_mesh_left_daughter.vtp',
    'vtk/smc_mesh_right_daughter.vtp'
]
input_mesh_files = {"ec" : ec_mesh_files, "smc" : smc_mesh_files}
    
base_names = {"ec" : 'solution/ec_data_t_', "smc" : 'solution/smc_data_t_', "atp" : 'solution/atp', "wss" : 'solution/wss'}


def reorder_species(species_array, reordered_array, cellType):
    """ Reorders the values in the species array given the quad to task ratio.
        Inserts this value into the reordered array."""
           
    # Given the joining of quads, calculate where the new larger quads
    # begin with respect to the position of the original ones.
    requiredQuads = []      
    
    for axialquad in range(axialQuads):
        for circQuad in range(circQuads):
            if circQuad % circScale == 0 and axialquad % axialScale == 0:
                requiredQuads.append(circQuad + (axialquad * circQuads))
                
    val = 0
    for i in requiredQuads: 
        quadOffset = i * numCells[cellType][0] * numCells[cellType][1]
    
        for extraAxial in range(0, axialScale):
            axialOffset = extraAxial * numCells[cellType][1] * numCells[cellType][0] * circQuads
            
            for j in range(0, numCells[cellType][0]): #smc rows
                rowOffset = j * numCells[cellType][1]
                
                for extraCirc in range(0, circScale):
                    circOffset = extraCirc * numCells[cellType][0] * numCells[cellType][1]               
                    
                    for k in range(0, numCells[cellType][1]): #smc cols
                        
                        cellId = quadOffset + rowOffset + circOffset + axialOffset + k
                        reordered_array.SetValue(cellId,species_array[val])
                        val += 1

def append_datasets(writers, h5_file_base, dataset_name):
    """Appends all data from a number of H5 files""" 
    
    # Create a dictionary of species arrays
    species_arrays = {}
    for attribute in attributes[output]:
        species_arrays[attribute] = []

    for writer in range(writers):
        h5_file_name = h5_file_base[:-4] + str(writer) + h5_file_base[-3:]
        
        fid = h5py.h5f.open(h5_file_name)
    
        dset = h5py.h5d.open(fid, dataset_name)
        shape = dset.shape
        rdata = numpy.zeros(shape[0], dtype=numpy.float64)
        dset.read(h5py.h5s.ALL, h5py.h5s.ALL, rdata)
    
        arr = rdata.ravel()
    
        for i in range(len(arr)): 
            species_arrays[attributes[output][i % len(attributes[output])]].append(arr[i])
    return species_arrays
        

def HDF5toVTKLumen():

    cellType = "ec"     # Both ATP and WSS maps use EC mesh    
    input_meshes = []

    # Read input EC meshes.
    for in_file in input_mesh_files[cellType]:
        reader = vtk.vtkXMLPolyDataReader()
        reader.SetFileName(in_file)
        reader.Update()

        input_meshes += [reader.GetOutput()]
                
        # Only add parent mesh for a tube (non-bifurcation)
        if branches == 1:
            break   
    
    atp_timesteps = len(glob.glob(base_names[output] + '_b_1_0_*'))
    print "Number of atp timesteps detected " + str(atp_timesteps)
    
    for timestep in range(atp_timesteps):
    
        append_filter = vtk.vtkAppendFilter()
            
        for branch in range(branches):
            
            species_array = []
            
            mesh = vtk.vtkPolyData()
            mesh.DeepCopy(input_meshes[branch])
    
            # The base input h5 filename given the branch and from which writer it came on said branch.
            h5_file_base = base_names[output] + '_b_' + str(branch + 1) + '_' + 'x_' + str(timestep) + '.h5'
            print "Processing file", h5_file_base
            for writer in range(writers):
                h5_file_name = h5_file_base[:-6] + str(writer) + h5_file_base[-5:]
                
                fid = h5py.h5f.open(h5_file_name)
            
                dset = h5py.h5d.open(fid, "data")
                shape = dset.shape
                rdata = numpy.zeros(shape[0], dtype=numpy.float64)
                dset.read(h5py.h5s.ALL, h5py.h5s.ALL, rdata)
                
                species_array += list(rdata.ravel())[:]
            
            
            reordered_array = vtk.vtkDoubleArray()
            reordered_array.SetName(output)
            reordered_array.SetNumberOfValues(numCells[cellType][0] * numCells[cellType][1] * circQuads * axialQuads)
            reorder_species(species_array, reordered_array, cellType)
            mesh.GetCellData().AddArray(reordered_array)
    
            append_filter.AddInputData(mesh)
            
        append_filter.Update()
    
        # Write the result.
        vtu_file = base_names[output] + '_' + str(timestep) + '.vtu'
        print 'Writing file', os.path.abspath(vtu_file)
    
        writer = vtk.vtkXMLUnstructuredGridWriter()
        writer.SetFileName(vtu_file)
        writer.SetInputData(append_filter.GetOutput())
        writer.Update()


def HDF5toVTKCells():

    input_meshes = []

    # Read input meshes.
    for in_file in input_mesh_files[output]:
        reader = vtk.vtkXMLPolyDataReader()
        reader.SetFileName(in_file)
        reader.Update()

        input_meshes += [reader.GetOutput()]
        
        # Only add parent mesh for a tube (non-bifurcation)
        if branches == 1:
            break   

    for time_step in range(args.start, args.end + 1):
        
        append_filter = vtk.vtkAppendFilter()
        
        for branch in range(branches):
            mesh = vtk.vtkPolyData()
            mesh.DeepCopy(input_meshes[branch])
            

            # The base input h5 filename given the branch and from which writer it came on said branch.
            h5_file_base = base_names[output] + str(time_step) + '_b_' + str(branch + 1) + '_' + 'x' + '.h5'
            print "Processing file", h5_file_base
            
            # Group all datasets of a branch at a specific time point given
            # the number of writers the data was split into.
            species_array = append_datasets(writers, h5_file_base, "data")
            
            # Loop through all attirbutes and append them to a new array in the 
            # correct order given the quad to task ratio.            
            for attribute in attributes[output]:
                reordered_array = vtk.vtkDoubleArray()
                reordered_array.SetName(attribute)
                reordered_array.SetNumberOfValues(numCells[output][0] * numCells[output][1] * circQuads * axialQuads)
                
                reorder_species(species_array[attribute], reordered_array, output)
                mesh.GetCellData().AddArray(reordered_array)
               
            append_filter.AddInputData(mesh)
            
        append_filter.Update()

        # Write the result.
        vtu_file = base_names[output] + str(time_step) + '.vtu'
        print 'Writing file', os.path.abspath(vtu_file)

        writer = vtk.vtkXMLUnstructuredGridWriter()
        writer.SetFileName(vtu_file)
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
    argParser.add_argument('output', type=str, help='EC, SMC, ATP, or WSS')
        
    # TODO: Take (most) values from simulation .ini file when it is created. Until then, lots of params...

    args = argParser.parse_args()    
    circScale = args.circScale
    axialScale = args.axialScale

    circQuads = args.circQuads
    axialQuads = args.axialQuads
    
    output = args.output.lower()
    
    print "Simulation output to convert: " + output
    
    if not glob.glob("solution/" + output + "*.h5"):
        print("No solution files for this cell type!")
        exit()

    # Find the number of writers and the number of branches given the 
    # output files.        
    writers = len(glob.glob("solution/ec*t_0_b_1*"))
    branches = len(glob.glob("solution/ec*t_0_b_*_0*"))

    print "Number of branches detected: " + str(branches)
    print "Number of writers detected: " + str(writers)
    
    if output in ["smc", "ec"]:
        HDF5toVTKCells()
        
    elif output in ["atp", "wss"]:
        HDF5toVTKLumen()

    else:
        print "Wrong output type specified..."
        exit()
