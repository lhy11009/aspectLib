#include <math.h>
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
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
#include <vtkCellData.h>

class FileReader
{
    public:
        void read_horiz_avg(const std::string &filename,
        std::vector<double> & depths,
        std::vector<double> &Ts,
        std::vector<double> &densities);

    private:

};

struct Prm
{
    //todo, parameters
    double grav_acc;  // value for gravity acceleration
    double Ro;
    double slab_depth;

};

class AspectVtk
{
    public:
        // initialize, todo
        AspectVtk(){
            prm.grav_acc = 10.0;  // value for gravity acceleration
            prm.Ro = 6371e3;
            prm.slab_depth = 50e3;
        };
        // readfile
        void readfile(std::string filename);
        // read horizontal average file
        void read_horiz_avg(const std::string &filename);
        // get poly data from reader
        void input_poly_data();
        // Triangulate the grid points
        void triangulate_grid();
        // Extract contour
        vtkSmartPointer<vtkPolyData> extract_contour(const std::string field_name, const double contour_value, const std::string filename);
        // Output
        void output(const vtkSmartPointer<vtkPolyData> opolydata, const std::string filename);
        // interpolate to uniform grid
        void interpolate_uniform_grid(std::string filename);
        // integrate on cells
        void integrate_cells();
        // derive density difference
        void density_diff();
        // derive mow from blocking temperature; 
        void mow_from_blocking(const double blocking_T, const double blocking_P);
        // get data from horiz
        double get_from_horiz(double depth, const std::string &field, const bool fix_out_of_bound=false);

        vtkSmartPointer<vtkXMLPUnstructuredGridReader> reader;
        // polydata read in from file
        vtkSmartPointer<vtkPolyData> iPolyData; 
        // triangulated data
        vtkSmartPointer<vtkDelaunay2D> iDelaunay2D;
        // reader
        FileReader file_reader;
        // horizontal average
        std::vector<double> horiz_depths;
        std::vector<double> horiz_Ts;
        std::vector<double> horiz_densities;
        // todo parameter
        Prm prm;
};

void read_options(const std::string &option_filename, std::vector<std::string> &options);