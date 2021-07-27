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
#include <vtkBandedPolyDataContourFilter.h>
#include <vtkContourFilter.h>
#include <vtkSimplePointsWriter.h>
#include <vtkTriangle.h>


class AspectVtk
{
    public:
        void readfile(std::string filename);
        // get poly data from reader
        void input_poly_data();
        // Triangulate the grid points
        void triangulate_grid();
        // Extract contour
        void extract_contour(std::string filename);
        // interpolate to uniform grid
        void interpolate_uniform_grid(std::string filename);
        // integrate on cells
        void integrate_cells();

    private:
        vtkSmartPointer<vtkXMLPUnstructuredGridReader> reader;
        // polydata read in from file
        vtkSmartPointer<vtkPolyData> iPolyData; 
        // triangulated data
        vtkSmartPointer<vtkDelaunay2D> iDelaunay2D;
};


void AspectVtk::readfile(std::string filename)
{
    std::cout << "Read file: " << filename << std::endl;
    reader = vtkSmartPointer<vtkXMLPUnstructuredGridReader>::New();
    reader->SetFileName(filename.c_str());
    reader->Update();
}


void AspectVtk::input_poly_data()
{
    std::cout << "Input poly data: " << std::endl;
    // point data
    vtkSmartPointer<vtkPoints> vtk_points = reader->GetOutput()->GetPoints(); // coordinates
    // field data
    vtkSmartPointer<vtkPointData> vtk_point_data = reader->GetOutputAsDataSet()->GetPointData();
    vtkSmartPointer<vtkDataArray> temperatures = vtk_point_data->GetArray("T"); // read temperature
    // set value of the private member
    iPolyData = vtkSmartPointer<vtkPolyData>::New();
    iPolyData->SetPoints(vtk_points);  //set as original vector
    iPolyData->GetPointData()->SetScalars(temperatures);
    // print information
    auto nl = iPolyData->GetNumberOfLines();
    auto np = iPolyData->GetNumberOfPoints();  
    auto npi = iPolyData->GetNumberOfPieces(); 
    auto nv = iPolyData->GetNumberOfVerts();  //debug
    std::cout << "PolyData info:" << std::endl
        << "number of lines: " << nl << ", "
        << "number of points: " << np << ", "
        << "number of pieces: " << npi << ", "
        << "number of verts: " << nv << std::endl;
}
  

void AspectVtk::triangulate_grid()
{
    std::cout << "Triangulate the grid points" << std::endl;
    iDelaunay2D = vtkSmartPointer<vtkDelaunay2D>::New();
    iDelaunay2D->SetInputData(iPolyData);
    iDelaunay2D->Update();
}


void AspectVtk::integrate_cells()
{
    std::cout << "Integrate on cells" << std::endl;
    vtkSmartPointer <vtkPolyData> tpolydata = iDelaunay2D->GetOutput();// polydata after trangulation
    std::cout << "\t cell type: " << tpolydata->GetCell(0)->GetCellType()
        << ", number of cells: " << tpolydata->GetNumberOfCells() << std::endl;  // check cell type
    // const double rin = 2890e3;  // rin
    double total_area = 0.0;
    for (vtkIdType i = 0; i < tpolydata->GetNumberOfCells(); i++)
    {
        vtkCell* cell = tpolydata->GetCell(i);  // cell in this polydata
        double p0[3], p1[3], p2[3], r0, r1, r2;
        cell->GetPoints()->GetPoint(0, p0);
        cell->GetPoints()->GetPoint(1, p1);
        cell->GetPoints()->GetPoint(2, p2);
        // r0 = sqrt(p0[0] * p0[0] + p0[1] * p0[1]);
        // r1 = sqrt(p1[0] * p1[0] + p1[1] * p1[1]);
        // r2 = sqrt(p2[0] * p2[0] + p2[1] * p2[1]);
        double area = vtkTriangle::TriangleArea(p0, p1, p2);
        total_area += area;
    }
    std::cout << "\t total area: " << total_area << std::endl;
}
        

void AspectVtk::extract_contour(std::string filename)
{
    std::cout << "Filter contour" << std::endl;
    vtkNew<vtkContourFilter> contour_filter;  // simple
    contour_filter->SetInputData(iDelaunay2D->GetOutput());  // set to the trangulated dataset
    contour_filter->Update();
    contour_filter->GenerateValues(1, 1173.0, 1173.0);
    contour_filter->Update();
    // write output 
    vtkNew<vtkSimplePointsWriter> writer;
    writer->SetInputData(contour_filter->GetOutput());
    writer->SetFileName(filename.c_str());
    writer->Update();
    writer->Write();
}
        

void AspectVtk::interpolate_uniform_grid(std::string filename)
{
    // Create a grid of points to interpolate over
    std::cout << "interpolate uniform grid" << std::endl;
    std::cout << "\tGenerate new grid" << std::endl;
    const int rSize = 1000;
    const int phiSize = 2000;
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
    std::cout << "\t Perform interpolation onto the new grid" << std::endl;
    vtkNew<vtkProbeFilter> probeFilter;
    probeFilter->SetSourceConnection(iDelaunay2D->GetOutputPort());
    probeFilter->SetInputData(gridPolyData); // Interpolate 'Source' at these points
    probeFilter->Update();
    
    // Map the output zvalues to the z-coordinates of the data so that
    // we get a surface, rather than a flat grid with interpolated
    // scalars.
    vtkNew<vtkWarpScalar> gridWarpScalar;
    gridWarpScalar->SetInputConnection(probeFilter->GetOutputPort());
    gridWarpScalar->Update();
    
    //output
    vtkNew<vtkXMLPolyDataWriter> writer;
    writer->SetFileName(filename.c_str());
    writer->SetInputConnection(gridWarpScalar->GetOutputPort());
    writer->Write();
}


int main(int argc, char* argv[])
{
  // parse command line arguments
  if (argc < 2)
  {
    std::cerr << "Usage: " << argv[0] << " Filename(.vtu) e.g. tetra.vtu"
              << std::endl;
    return EXIT_FAILURE;
  }
  // read file
  std::string filename = argv[1];
  AspectVtk aspect_vtk;
  aspect_vtk.readfile(filename);
  aspect_vtk.input_poly_data();
  aspect_vtk.triangulate_grid();
  aspect_vtk.integrate_cells();
  // aspect_vtk.extract_contour("contour.txt");
  // aspect_vtk.interpolate_uniform_grid("uniform2D.vtp");  // intepolation
  return EXIT_SUCCESS;
}
