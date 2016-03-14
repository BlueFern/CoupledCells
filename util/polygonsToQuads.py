from vtk import *
import re
import glob
import argparse

def sortNicely(l): 
    """ Sort the given list in the way that humans expect."""
     
    convert = lambda text: int(text) if text.isdigit() else text 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    l.sort( key=alphanum_key )
    
    
def polygonsToQuads(path, cell_type, start, end):

    # Get all files
    base_files = glob.glob(path + cell_type + "*")
    sortNicely(base_files)
    
    #print "Found files:", base_files
    
    start_index = -1
    end_index = -1
    for i in range(len(base_files)):
        if str(start) + ".vtu" in base_files[i]:
            start_index = i
        if str(end) + ".vtu" in base_files[i]:
            end_index = i

    if start_index == -1:
        exit("First time-step not found in files")
    if end_index == -1:
        exit("Last time-step not found in files")
    

    for f in base_files[start_index:end_index + 1]:
        print "Correcting:", f
        
        
        
        # Read in file and covert to polydata
        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(f)
        reader.Update()
        
        geometryFilter = vtk.vtkGeometryFilter()
        geometryFilter.SetInputData(reader.GetOutput())
        geometryFilter.Update()
        
        polyData = geometryFilter.GetOutput()
        
        
        newCells = vtk.vtkCellArray()
        oldCellData = polyData.GetCellData()
        numberOfArrays = oldCellData.GetNumberOfArrays()
        
        oldCellDataArrays = []
        newCellDataArrays = []
        
        # Get all previous cell data and initialise new cell data arrays with their
        # appropriate names.
        for i in range(numberOfArrays):
            oldCellDataArrays.append(oldCellData.GetArray(i))
            newCellDataArrays.append(vtkDoubleArray())
            newCellDataArrays[i].SetName(oldCellDataArrays[i].GetName())
            
        
        oldCell = vtk.vtkIdList()
        newCell = vtk.vtkIdList()
        newCell.SetNumberOfIds(4)
        
        # Iterate over all cells. Create two quads from all 6-sided polygons 
        # and give both new quads the same cell data as the original polygon.
        for i in range(polyData.GetNumberOfCells()):
            
            if polyData.GetCell(i).GetNumberOfPoints() == 6:
                
                polyData.GetCellPoints(i, oldCell)
                
                newCell.SetId(0, oldCell.GetId(0))
                newCell.SetId(1, oldCell.GetId(1))
                newCell.SetId(2, oldCell.GetId(4))
                newCell.SetId(3, oldCell.GetId(5))
                newCells.InsertNextCell(newCell)
                
                newCell.SetId(0, oldCell.GetId(1))
                newCell.SetId(1, oldCell.GetId(2))
                newCell.SetId(2, oldCell.GetId(3))
                newCell.SetId(3, oldCell.GetId(4))
                newCells.InsertNextCell(newCell)
                
                for j in range(numberOfArrays):
                    newCellDataArrays[j].InsertNextTuple(oldCellDataArrays[j].GetTuple(i))
                    newCellDataArrays[j].InsertNextTuple(oldCellDataArrays[j].GetTuple(i))
                    
            # Already a quad
            else:
                newCells.InsertNextCell(polyData.GetCell(i))
                
                for j in range(numberOfArrays):
                    newCellDataArrays[j].InsertNextTuple(oldCellDataArrays[j].GetTuple(i))
                
        # Create new polydata with new cells and cell data.
        newPolyData = vtk.vtkPolyData()
        newPolyData.SetPoints(polyData.GetPoints())
        newPolyData.SetPolys(newCells)
        
        for j in range(numberOfArrays):
            newPolyData.GetCellData().AddArray(newCellDataArrays[j])
            
        # Write out new file
            
        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName(f[:-4] + "_new.vtp")
        writer.SetInputData(newPolyData)
        writer.Update()
        

        

if __name__ == "__main__":
    
    argParser = argparse.ArgumentParser(
    description='Convert 6 sided polygons to 2 quads with identical cell data in vtu files.')
    argParser.add_argument('path', type=str, help='path')
    argParser.add_argument('cell_type', type=str, help='cell type')
    argParser.add_argument('start', type=int, help='first time-step')
    argParser.add_argument('end', type=int, help='last time-step')

    args = argParser.parse_args()    
    polygonsToQuads(args.path, args.cell_type, args.start, args.end)
