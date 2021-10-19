#include <math.h>
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include "AspectVTK.h"

void QuickSort(const std::vector<double> & A, std::vector<size_t> & I, size_t lo, size_t hi);

class SlabAnalysis : public AspectVtk
{
    public:
        void prepare_slab(const std::vector<std::string> names);
};

//todo
struct SlabOutputs
{
    double trench_theta; // trench position
    double slab_depth;  // slab depth
    void set_moprh(double trench_theta, double slab_depth);
};

//todo
void SlabOutputs::set_moprh(double trench_theta, double slab_depth){
    this->trench_theta = trench_theta;
    this->slab_depth = slab_depth;
}


void SlabAnalysis::prepare_slab(const std::vector<std::string> names)
{
    std::cout << "Praparing the slab (new composition)" << std::endl;
    vtkSmartPointer <vtkPolyData> tpolydata = iDelaunay2D->GetOutput();// polydata after trangulation
    std::vector<vtkFloatArray*> compositions;
    for (auto &p: names)
    {
        compositions.push_back(dynamic_cast<vtkFloatArray*>(tpolydata->GetPointData()->GetArray(p.c_str())));
    }
    vtkNew<vtkFloatArray> slab_fields;
    slab_fields->DeepCopy(*(compositions.begin())); // the first field is copied into the new vector
    slab_fields->SetName("slab");
    for (vtkIdType i = 0; i < tpolydata->GetNumberOfPoints(); i++)
    {
        for (auto iter=compositions.begin()+1; iter<compositions.end(); iter++)
        {
            // second to later fields are added
            const double new_comp = (*iter)->GetTuple1(i);
            slab_fields->SetTuple1(i, slab_fields->GetTuple1(i) + new_comp);
        }
    }
    tpolydata->GetPointData()->AddArray(slab_fields);
    iDelaunay2D->Update();
}

void analyze_slab(vtkSmartPointer<vtkPolyData> c_poly_data, SlabOutputs & slab_outputs)
{
  // take a contour output and analyze the slab morphology
  vtkSmartPointer<vtkPoints> c_points = c_poly_data->GetPoints();
  // intiate pointers
  auto rs = std::make_shared<std::vector<double>>();
  auto thetas = std::make_shared<std::vector<double>>();
  auto I = std::make_shared<std::vector<size_t>>();
  for (vtkIdType id = 0; id < c_points->GetNumberOfPoints(); id++)
  {
    double *p = c_points -> GetPoint(id);
    double x = p[0];
    double y = p[1];
    double r = sqrt(x*x + y*y);
    double theta = acos(x/r);
    rs->push_back(r);
    thetas->push_back(theta);
    I->push_back(id);
    // std::cout<< "id: " << id << ",r: " << r << ", theta: " << theta << std::endl;  // debug
  }
  QuickSort(*rs, *I, 0, c_points->GetNumberOfPoints()-1);  // reorder by r, index saved in I
  // found trench point
  double ro = 6371e3;
  double find_depth = 5e3;
  double n_aver = 3;
  auto r_deviation = std::make_shared<std::vector<double>>();
  auto I_deviation = std::make_shared<std::vector<size_t>>();
  for (vtkIdType id = 0; id < c_points->GetNumberOfPoints(); id++)
  {
    r_deviation->push_back(abs((*rs)[id] - ro + find_depth));
    I_deviation->push_back(id);
  }
  QuickSort(*r_deviation, *I_deviation, 0, c_points->GetNumberOfPoints()-1);
  std::cout << c_points->GetNumberOfPoints() << ", " << (*I_deviation)[0] << std::endl;
  std::cout << (*thetas)[(*I_deviation)[0]] << ", " << (*thetas)[(*I_deviation)[1]] << ", " << (*thetas)[(*I_deviation)[2]] << std::endl;
  double sum = 0.0;
  for (int i=0; i < n_aver; i++)
  {
    sum += (*thetas)[(*I_deviation)[0]];
  }
  double trench_theta = sum / n_aver;
  std::cout << "Trench theta: " << trench_theta << std::endl;
  // find slab depth
  double slab_depth = ro - (*rs)[(*I)[0]];
  std::cout << "Slab depth: " << slab_depth << std::endl;
  // todo return results
  slab_outputs.set_moprh(trench_theta, slab_depth);
}


void analyze_temperature(SlabAnalysis & slab_analysis, SlabOutputs & slab_outputs, const std::string filename)
{
  //todo
  vtkNew<vtkPoints> gridPoints;
  double Ro = 6371e3;
  double depth = 1000e3;  //debug
  unsigned num = 100;
  for (unsigned i = 0; i < num; i++){
    // a vertical line
    double val_theta = slab_outputs.trench_theta;
    double val_r = Ro - i * 1.0 / num * (depth);
    double val_x = val_r * cos(val_theta);
    double val_y = val_r * sin(val_theta);
    gridPoints->InsertNextPoint(val_x, val_y, 0);
  }
  vtkSmartPointer<vtkPolyData> iPolyData = slab_analysis.interpolate_grid(gridPoints);
  std::vector<std::string> fields{"T", "spcrust"};
  slab_analysis.write_ascii(iPolyData, fields, filename);
}


void QuickSort(const std::vector<double> & A, std::vector<size_t> & I, size_t lo, size_t hi)
{
    if (lo < hi)
    {
        double pivot = A[I[lo + (hi - lo) / 2]];
        size_t t;
        size_t i = lo - 1;
        size_t j = hi + 1;
        while (1)
        {
            while (A[I[++i]] < pivot);
            while (A[I[--j]] > pivot);
            if (i >= j)
                break;
            t = I[i];
            I[i] = I[j];
            I[j] = t;
        }
        QuickSort(A, I, lo, j);
        QuickSort(A, I, j + 1, hi);
    }
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
  // read options from file
  std::vector<std::string> options;
  std::string option_file(argv[1]);
  read_options(argv[1], options);
  if (options.size() < 2)
  {
    std::cerr << "option file has to have 2 variables"
              << std::endl;
    return EXIT_FAILURE;
  }
  // 1. PVTU_FILE_PATH
  std::string filename = options[0]; // get dir
  size_t found=option_file.find_last_of("/\\");
  std::string target_dir = option_file.substr(0,found);
  found=filename.find_last_of("."); // get step
  size_t found1 = filename.find_last_of("-");
  std::string pvtu_step = filename.substr(found1+1, found - found1-1);
  std::cout << "pvtu_step: " << pvtu_step << std::endl;
  // 2. AVG_FILE_PATH
  std::string avg_filename = options[1];
  SlabAnalysis slab_analysis;
  slab_analysis.readfile(filename);
  slab_analysis.read_horiz_avg(avg_filename);
  slab_analysis.input_poly_data();
  slab_analysis.triangulate_grid();
  //slab_analysis.density_diff();
  //slab_analysis.mow_from_blocking(973.0, 12.5e9);  // 725 + 273 from Quinteros
  slab_analysis.prepare_slab({"spcrust", "spharz"});
  slab_analysis.integrate_cells();
  // slab_analysis.extract_contour("T", 1173.0, target_dir + "/" + "contour.txt");
  vtkSmartPointer<vtkPolyData> c_poly_data = slab_analysis.extract_contour("slab", 0.99, target_dir + "/" + "contour_slab_" + pvtu_step + ".txt"); //apply contour
  SlabOutputs slab_outputs;
  analyze_slab(c_poly_data, slab_outputs); // analyze slab morphology
  //todo
  analyze_temperature(slab_analysis, slab_outputs, target_dir + "/" + "wedge_temperature_" + pvtu_step + ".txt"); // output tempertature
  //aspect_vtk.interpolate_uniform_grid("uniform2D.vtp");  // intepolation
  return EXIT_SUCCESS;
}