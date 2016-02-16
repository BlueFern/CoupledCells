"""
This script corrects the artifacts shown in meshes due to openGL problems
rendering non-convex, non-planar polygons. It replaces such polygons
(6-sided) with two quads, and duplicates the cell data. The output files
are in the same directory with the same filename with "new" at the start.

This is intended to be used on the vtu files created after a simulation has
been run and the h5 files converted to vtu.

The script is potentially quite expensive for large meshes over long 
simulation times.

"""

from vtk import *
import re
import sys
import glob

BASE_PATTERN = "_data_t_*.vtu"

def sortNicely(l): 
    """ Sort the given list in the way that humans expect."""
     
    convert = lambda text: int(text) if text.isdigit() else text 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    l.sort( key=alphanum_key )
    
    
def main():

    files = glob.glob(sys.argv[1] + '/' + sys.argv[2] + BASE_PATTERN)
    sortNicely(files);
    
    reader = vtkXMLUnstructuredGridReader()
    polyData = vtkPolyData()
    
    
    for i in range(int(sys.argv[3]), int(sys.argv[4])):
        
        
        reader.SetFileName(files[i])
        reader.Update()
        polyData = reader.GetOutput()
        
        newPolyData = vtkPolyData()
        cells = vtkCellArray()
        
        oldCellData = polyData.GetCellData()
        
        numberOfArrays = oldCellData.GetNumberOfArrays()
        
        oldCellDataArrays = []
        newCellDataArrays = []
        
        # Get all previous cell data and initialise new cell data arrays with their
	   # appropriate names.
        for j in range(numberOfArrays):
            oldCellDataArrays.append(oldCellData.GetArray(j))
            newCellDataArrays.append(vtkDoubleArray())
            newCellDataArrays[-1].SetName(oldCellDataArrays[-1].GetName())
            
        oldCell = vtkIdList()
        newCell = vtkIdList()
        newCell.SetNumberOfIds(4)
        
        # Iterate over all cells. Create two quads from all 6-sided polygons and give both
	   # new quads the same cell data as the original polygon.
        for j in range (polyData.GetNumberOfCells()):
            
            if polyData.GetCell(j).GetNumberOfPoints() == 6:
                
                polyData.GetCellPoints(j, oldCell)
                
                newCell.SetId(0, oldCell.GetId(0))
                newCell.SetId(1, oldCell.GetId(1))
                newCell.SetId(2, oldCell.GetId(4))
                newCell.SetId(3, oldCell.GetId(5))
                cells.InsertNextCell(newCell)
                
                newCell.SetId(0, oldCell.GetId(1))
                newCell.SetId(1, oldCell.GetId(2))
                newCell.SetId(2, oldCell.GetId(3))
                newCell.SetId(3, oldCell.GetId(4))
                cells.InsertNextCell(newCell)
                
                for k in range(numberOfArrays):
                    newCellDataArrays[k].InsertNextTuple(oldCellDataArrays[k].GetTuple(j))
                    newCellDataArrays[k].InsertNextTuple(oldCellDataArrays[k].GetTuple(j))
                    
            # Already a quad, copy across.
            else:
                cells.InsertNextCell(polyData.GetCell(j))
                for k in range(numberOfArrays):
                    newCellDataArrays[k].InsertNextTuple(oldCellDataArrays[k].GetTuple(j))
        
        
        newPolyData.SetPoints(polyData.GetPoints())
        newPolyData.SetPolys(cells)
        
        for j in range(numberOfArrays):
            newPolyData.GetCellData().AddArray(newCellDataArrays[j])
        
        writer = vtkXMLPolyDataWriter()
        writer.SetInputData(newPolyData)
        writer.SetFileName(sys.argv[1] + "/new_" + files[i][9:-4] + ".vtp")
        print sys.argv[1] + "/new_" + files[i][9:-4] + ".vtp"
        writer.Update()


if __name__ == "__main__":
    
    if len(sys.argv) != 5:  
        sys.exit("Expected arguments: <data dir> <ec|smc> <start> <end>")

    main()