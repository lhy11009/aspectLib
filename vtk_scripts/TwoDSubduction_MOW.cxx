#include <math.h>
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include "AspectVTK.h"

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
  read_options(argv[1], options);
  if (options.size() < 2)
  {
    std::cerr << "option file has to have 2 variables"
              << std::endl;
    return EXIT_FAILURE;
  }
  // 1. PVTU_FILE_PATH
  std::string filename = options[0];
  // 2. VTK_OUTPUT_DIR
  std::string output_dir = options[1];
  // todo
  AspectVtk aspect_vtk;
  aspect_vtk.readfile(filename);
  aspect_vtk.input_poly_data();
  aspect_vtk.triangulate_grid();
  //aspect_vtk.density_diff();
  aspect_vtk.mow_from_blocking(998.0, 12.5e9);  // 725 + 273 from Quinteros
  //aspect_vtk.prepare_slab({"spcrust", "spharz"});
  aspect_vtk.output(output_dir + '/' + "output.vtp");
  //aspect_vtk.integrate_cells();
  //aspect_vtk.extract_contour("contour.txt");
  //aspect_vtk.interpolate_uniform_grid("uniform2D.vtp");  // intepolation
  return EXIT_SUCCESS;
}