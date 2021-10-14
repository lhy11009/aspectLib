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

void analyze_slab(vtkSmartPointer<vtkPolyData> c_poly_data)
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
  QuickSort(*rs, *I, 0, c_points->GetNumberOfPoints());  // reorder by r, index saved in I
  // todo: found trench point
  for (auto &p:*I)
    std::cout << p << ", " << (*rs)[p] << std::endl;  // debug
}

//todo: found trench point

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
  std::string step = filename.substr(found1+1, found - found1-1);
  std::cout << "step: " << step << std::endl;
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
  slab_analysis.output(slab_analysis.iDelaunay2D->GetOutput(), target_dir + "/" + "output.vtp");
  slab_analysis.integrate_cells();
  // slab_analysis.extract_contour("T", 1173.0, target_dir + "/" + "contour.txt");
  vtkSmartPointer<vtkPolyData> c_poly_data = slab_analysis.extract_contour("slab", 0.99, target_dir + "/" + "contour_slab_" + step + ".txt"); //apply contour
  analyze_slab(c_poly_data); // analyze slab morphology
  //aspect_vtk.interpolate_uniform_grid("uniform2D.vtp");  // intepolation
  return EXIT_SUCCESS;
}