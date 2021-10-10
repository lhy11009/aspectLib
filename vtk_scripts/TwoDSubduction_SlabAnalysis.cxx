#include <math.h>
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include "AspectVTK.h"

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
  std::string filename = options[0];
  size_t found=option_file.find_last_of("/\\");
  std::string target_dir = option_file.substr(0,found);
  std::cout << "target_dir: " << target_dir << std::endl;
  // 2. AVG_FILE_PATH
  std::string avg_filename = options[1];
  // todo
  SlabAnalysis slab_analysis;
  slab_analysis.readfile(filename);
  slab_analysis.read_horiz_avg(avg_filename);
  slab_analysis.input_poly_data();
  slab_analysis.triangulate_grid();
  //slab_analysis.density_diff();
  //slab_analysis.mow_from_blocking(973.0, 12.5e9);  // 725 + 273 from Quinteros
  slab_analysis.prepare_slab({"spcrust", "spharz"});
  slab_analysis.output(target_dir + "output.vtp");
  slab_analysis.integrate_cells();
  slab_analysis.extract_contour(target_dir + "contour.txt");
  //aspect_vtk.interpolate_uniform_grid("uniform2D.vtp");  // intepolation
  return EXIT_SUCCESS;
}