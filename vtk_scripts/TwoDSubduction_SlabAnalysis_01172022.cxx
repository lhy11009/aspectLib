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

struct SlabOutputs
{
    double trench_theta; // trench position
    double slab_depth;  // slab depth
    // todo
    double theta100; // angle where the slab is 100 km depth
    double dip100;
    void set_moprh(double trench_theta, double slab_depth);
};

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

std::shared_ptr<std::vector<size_t>> find_points_at_radius(std::shared_ptr<std::vector<double>> rs, const double r, const unsigned n)
{
  // todo
  // found trench point
  auto r_deviation = std::make_shared<std::vector<double>>();
  auto I_deviation = std::make_shared<std::vector<size_t>>();
  auto I_points = std::make_shared<std::vector<size_t>>();
  for (size_t id = 0; id < rs->size(); id++)
  {
    r_deviation->push_back(abs((*rs)[id] - r));
    I_deviation->push_back(id);
  }
  QuickSort(*r_deviation, *I_deviation, 0, rs->size()-1);
  for (int i=0; i < n; i++){
    I_points->push_back((*I_deviation)[i]);
  }
  return I_points;
}

void analyze_slab(vtkSmartPointer<vtkPolyData> c_poly_data, SlabOutputs & slab_outputs, const double find_trench_depth)
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
  double n_aver = 3;
  double find1 = find_trench_depth;
  // todo
  double diff, diff1, trench_theta;
  const double trench_theta_ref = 0.63;
  do {
    auto Ip = find_points_at_radius(rs, ro-find1, n_aver);
    diff =  abs((*thetas)[(*Ip)[n_aver-1]] - (*thetas)[(*Ip)[0]]);
    double sum = 0.0;
    for (int i=0; i < n_aver; i++)
    {
      sum += (*thetas)[(*Ip)[i]];
    }
    trench_theta = sum / n_aver;
    diff1 = abs(trench_theta - trench_theta_ref);
    std::cout << find1 <<", " << (*thetas)[(*Ip)[0]] << ", " << (*thetas)[(*Ip)[n_aver-1]] << std::endl; //debug
    find1 += 1e3;  // look for a deeper depth
  }while(diff >0.02 || diff1 > 0.2);
  auto Ip = find_points_at_radius(rs, ro-find_trench_depth, n_aver);
  // find_points_at_radius(double r)
  std::cout << "Trench theta: " << trench_theta << std::endl;
  // find slab depth
  double slab_depth = ro - (*rs)[(*I)[0]];
  std::cout << "Slab depth: " << slab_depth << std::endl;
  // find points at 100e3 depth
  double find2 = 100e3;
  n_aver = 10;
  Ip = find_points_at_radius(rs, ro-find2, n_aver);
  double theta100 = (*thetas)[(*Ip)[0]];
  for (unsigned i=0; i < n_aver; i++){
    // get the biggest value for a point on the surface
    if ((*thetas)[(*Ip)[i]] > theta100){
      theta100 = (*thetas)[(*Ip)[i]];
    }
  }
  slab_outputs.theta100 = theta100;
  slab_outputs.dip100 = atan(ro*(theta100 - trench_theta)/find2);
  std::cout << "100km theta: " << theta100 << std::endl;
  std::cout << "100km dip: " << slab_outputs.dip100 << std::endl;
  slab_outputs.set_moprh(trench_theta, slab_depth);
}


void analyze_wedge_temperature100(SlabAnalysis & slab_analysis, SlabOutputs & slab_outputs, const std::string filename)
{
  // wedge temperature above where the slab is 100 km deep
  vtkNew<vtkPoints> gridPoints;
  double Ro = 6371e3;
  double depth = 120e3;  // this is set > 100 so as to check for slab compositions
  unsigned num = 120;
  for (unsigned i = 0; i < num; i++){
    // a vertical line
    double val_theta = slab_outputs.theta100;
    double val_r = Ro - i * 1.0 / num * (depth);
    double val_x = val_r * cos(val_theta);
    double val_y = val_r * sin(val_theta);
    gridPoints->InsertNextPoint(val_x, val_y, 0);
  }
  vtkSmartPointer<vtkPolyData> iPolyData = slab_analysis.interpolate_grid(gridPoints);
  std::vector<std::string> fields{"T", "slab"};
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
  if (options.size() < 3)
  {
    std::cerr << "option file has to have 3 variables:" << std::endl
    << "0: path of the pvtu file." << std::endl
    << "1: path of the depth_average file" << std::endl
    << "2: find_trench_depth"
    << std::endl;
    return EXIT_FAILURE;
  }
  // 0. PVTU_FILE_PATH
  std::string filename = options[0]; // get dir
  size_t found=option_file.find_last_of("/\\");
  std::string target_dir = option_file.substr(0,found);
  found=filename.find_last_of("."); // get step
  size_t found1 = filename.find_last_of("-");
  std::string pvtu_step = filename.substr(found1+1, found - found1-1);
  std::cout << "pvtu_step: " << pvtu_step << std::endl;
  // 1. AVG_FILE_PATH
  std::string avg_filename = options[1];
  // 2. find_trench_depth
  double find_trench_depth = std::stod(options[2]); 
  SlabAnalysis slab_analysis;
  slab_analysis.readfile(filename);
  slab_analysis.read_horiz_avg(avg_filename);
  slab_analysis.input_poly_data();
  slab_analysis.triangulate_grid();
  //slab_analysis.density_diff();
  //slab_analysis.mow_from_blocking(973.0, 12.5e9);  // 725 + 273 from Quinteros
  slab_analysis.prepare_slab({"spcrust", "spharz"});
  slab_analysis.integrate_cells();
  vtkSmartPointer<vtkPolyData> c_poly_data = slab_analysis.extract_contour("slab", 0.99, target_dir + "/" + "contour_slab_" + pvtu_step + ".txt"); //apply contour
  SlabOutputs slab_outputs;
  analyze_slab(c_poly_data, slab_outputs, find_trench_depth); // analyze slab morphology
  analyze_wedge_temperature100(slab_analysis, slab_outputs, target_dir + "/" + "wedge_T100_" + pvtu_step + ".txt"); // output tempertature
  //aspect_vtk.interpolate_uniform_grid("uniform2D.vtp");  // intepolation
  return EXIT_SUCCESS;
}