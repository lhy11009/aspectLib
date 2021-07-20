#include <math.h>
#include <vtkActor.h>
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
  //todo read data sets
  vtkSmartPointer<vtkPointData> vtk_point_data = reader->GetOutputAsDataSet()->GetPointData();
  vtkSmartPointer<vtkDataArray> temperatures = vtk_point_data->GetArray("T"); // read temperature
  std::cout << "T[0] = " << temperatures->GetTuple(0) << std::endl;

  // read x and y from file
  // get a random z array
  std::cout << "Construct new polydata" << std::endl;
  vtkNew<vtkFloatArray> zvalues;
  zvalues->SetName("ZValues");
  vtkNew<vtkMinimalStandardRandomSequence> randomSequence;
  randomSequence->SetSeed(8775070);
  double maxHeight=1.0e6;
  for (unsigned int i = 0; i < n_points; ++i)
  {
    double z;
    // random position and radius
    z = randomSequence->GetRangeValue(0, maxHeight);
    randomSequence->Next();
    zvalues->InsertNextValue(z);
  }
  // a polydata 
  vtkNew<vtkPolyData> randomPolyData;
  randomPolyData->SetPoints(vtk_points);  //set as original vector
  //randomPolyData->GetPointData()->SetScalars(zvalues);
  //randomPolyData->GetPointData()->AddArray(temperatures);
  randomPolyData->GetPointData()->SetScalars(temperatures);
  // todo set field
  
  // Triangulate the grid points. If you do not have a mesh (points
  // only), the output will not be interpolated!
  std::cout << "Triangulate the grid points" << std::endl;
  vtkNew<vtkDelaunay2D> randomDelaunay;
  randomDelaunay->SetInputData(randomPolyData);
  randomDelaunay->Update();
  
  // Create a grid of points to interpolate over
  std::cout << "Generate new grid" << std::endl;
  const int rSize = 100;
  const int phiSize = 200;
  const double rMax = 6370e3;
  const double rMin = 2980e3;
  const double rIntr = (rMax - rMin) / double(rSize);
  const double phiMax = 1.064;
  const double phiMin = 0.0;
  const double phiIntr = (phiMax - phiMin) / double(phiSize);
  vtkNew<vtkPoints> gridPoints;
  for (unsigned int i = 0; i < rSize; i++)
  {
    for (unsigned int j = 0; j < phiSize; j++)
    {
      double val_r = rMin + i*rIntr;
      double val_phi = phiMin + j*phiIntr;
      double val_x = val_r * cos(val_phi);
      double val_y = val_r * sin(val_phi);

      gridPoints->InsertNextPoint(val_x, val_y, 0);
    }
  }
  // Create a dataset from the grid points
  vtkNew<vtkPolyData> gridPolyData;
  gridPolyData->SetPoints(gridPoints);
  
  // Perform the interpolation
  std::cout << "Perform interpolation onto the new grid" << std::endl;
  vtkNew<vtkProbeFilter> probeFilter;
  probeFilter->SetSourceConnection(randomDelaunay->GetOutputPort());
  probeFilter->SetInputData(
      gridPolyData); //
                     // Interpolate 'Source' at these points
  probeFilter->Update();
  // Map the output zvalues to the z-coordinates of the data so that
  // we get a surface, rather than a flat grid with interpolated
  // scalars.
  vtkNew<vtkWarpScalar> gridWarpScalar;
  gridWarpScalar->SetInputConnection(probeFilter->GetOutputPort());
  gridWarpScalar->Update();

  //////// Setup outputs ////////
  // Output random points
  // Map the output zvalues to the z-coordinates of the data
  std::cout << "Output data" << std::endl;
  vtkNew<vtkWarpScalar> randomWarpScalar;
  randomWarpScalar->SetInputConnection(randomDelaunay->GetOutputPort());
  randomWarpScalar->Update();
  vtkNew<vtkXMLPolyDataWriter> randomWriter;
  randomWriter->SetFileName("randomSurface.vtp");
  randomWriter->SetInputConnection(randomWarpScalar->GetOutputPort());
  randomWriter->Write();
  // Mesh the output grid points
  vtkNew<vtkDelaunay2D> gridDelaunay;
  gridDelaunay->SetInputConnection(gridWarpScalar->GetOutputPort());
  vtkNew<vtkXMLPolyDataWriter> gridWriter;
  gridWriter->SetFileName("gridSurface.vtp");
  gridWriter->SetInputConnection(gridDelaunay->GetOutputPort());
  gridWriter->Write();

  // write the original data set for comparison
  /*
  vtkSmartPointer<vtkUnstructuredGrid> un_grid = reader->GetOutput();
  vtkNew<vtkXMLUnstructuredGridWriter> unstructureWriter;
  unstructureWriter->SetFileName("unstructure.vtp");
  unstructureWriter->SetInputData(un_grid);
  unstructureWriter->Write();
  */

  return EXIT_SUCCESS;
}
