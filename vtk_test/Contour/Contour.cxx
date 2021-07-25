#include <math.h>
#include <vtkDataSet.h>
#include <vtkDataSetMapper.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkProperty.h>
#include <vtkPointSet.h>
#include <vtkPointData.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLPUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPolyData.h>
#include <vtkFloatArray.h>
#include <vtkMinimalStandardRandomSequence.h>
#include <vtkDelaunay2D.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkWarpScalar.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkProbeFilter.h>

int main(int argc, char* argv[])
{
  // parse command line arguments
  if (argc < 2)
  {
    std::cerr << "Usage: " << argv[0] << " Filename(.vtu) e.g. tetra.vtu"
              << std::endl;
    return EXIT_FAILURE;
  }

  std::string filename = argv[1];

  // read data from a vtu file
  std::cout << "Read data: " << filename << std::endl;
  // vtkNew<vtkXMLUnstructuredGridReader> reader;  // for vtu file
  vtkNew<vtkXMLPUnstructuredGridReader> reader; // for pvtu file
  reader->SetFileName(filename.c_str());
  reader->Update();
  //Read data points, GetOutput(): vtkUnstructuredGrid
  vtkSmartPointer<vtkPoints> vtk_points = reader->GetOutput()->GetPoints(); // coordinates
  int n_points = vtk_points->GetNumberOfPoints();  // use size of the vtk_points object
  std::cout << "n_points: " << n_points << std::endl;
  //read data sets
  vtkSmartPointer<vtkPointData> vtk_point_data = reader->GetOutputAsDataSet()->GetPointData();
  vtkSmartPointer<vtkDataArray> temperatures = vtk_point_data->GetArray("T"); // read temperature
  std::cout << "T[0] = " << temperatures->GetTuple(0) << std::endl;

  // read x and y from file
  // a polydata 
  vtkNew<vtkPolyData> randomPolyData;
  randomPolyData->SetPoints(vtk_points);  //set as original vector
  randomPolyData->GetPointData()->SetScalars(temperatures);
  
  // Triangulate the grid points. If you do not have a mesh (points
  // only), the output will not be interpolated!
  std::cout << "Triangulate the grid points" << std::endl;
  vtkNew<vtkDelaunay2D> randomDelaunay;
  randomDelaunay->SetInputData(randomPolyData);
  randomDelaunay->Update();
  
  //////// Setup outputs ////////
  // Output random points
  // Map the output zvalues to the z-coordinates of the data
  std::cout << "Output data" << std::endl;
  vtkNew<vtkWarpScalar> randomWarpScalar;
  randomWarpScalar->SetInputConnection(randomDelaunay->GetOutputPort());
  randomWarpScalar->Update();
  vtkNew<vtkXMLPolyDataWriter> randomWriter;
  randomWriter->SetFileName("Subduction2D.vtp");
  randomWriter->SetInputConnection(randomWarpScalar->GetOutputPort());
  randomWriter->Write();


  return EXIT_SUCCESS;
}
