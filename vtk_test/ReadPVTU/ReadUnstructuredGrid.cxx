#include <vtkActor.h>
#include <vtkDataSetMapper.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkProperty.h>
#include <vtkPolyData.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLPUnstructuredGridReader.h>
#include <vtkXMLPUnstructuredDataReader.h>

int main(int argc, char* argv[])
{
  // parse command line arguments
  if (argc < 2)
  {
    std::cerr << "Usage: " << argv[0] << " Filename(.pvtu) e.g. tetra.vtu"
              << std::endl;
    return EXIT_FAILURE;
  }

  std::string filename = argv[1];

  // read all the data from the file
  vtkNew<vtkXMLPUnstructuredGridReader> reader;
  reader->SetFileName(filename.c_str());
  reader->Update();

  // get grid
  // vtkUnstructuredGrid* ugrid = reader->GetOutput();
  vtkNew<vtkPolyData> polyData;
  polyData->SetPoints(reader->GetOutput()->GetPoints());

  // read field data
  vtkNew<vtkXMLPUnstructuredGridReader> reader1;
  //vtkNew<vtkXMLPUnstructuredDataReader> DataReader;
  // DataReader->SetFileName(filename.c_str());
  // DataReader->Update();

/*
  vtkNew<vtkNamedColors> colors;

  // Create a mapper and actor
  vtkNew<vtkDataSetMapper> mapper;
  mapper->SetInputConnection(reader->GetOutputPort());
  mapper->ScalarVisibilityOff();

  vtkNew<vtkActor> actor;
  actor->SetMapper(mapper);
  actor->GetProperty()->EdgeVisibilityOn();
  actor->GetProperty()->SetLineWidth(2.0);
  actor->GetProperty()->SetColor(colors->GetColor3d("MistyRose").GetData());

  vtkNew<vtkProperty> backFace;
  backFace->SetColor(colors->GetColor3d("Tomato").GetData());
  actor->SetBackfaceProperty(backFace);

  // Create a renderer, render window, and interactor
  vtkNew<vtkRenderer> renderer;
  vtkNew<vtkRenderWindow> renderWindow;
  renderWindow->AddRenderer(renderer);
  vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
  renderWindowInteractor->SetRenderWindow(renderWindow);

  // Add the actor to the scene
  renderer->AddActor(actor);
  renderer->SetBackground(colors->GetColor3d("Wheat").GetData());

  // Render and interact
  renderWindow->SetSize(640, 480);

  renderWindow->Render();
  renderWindowInteractor->Start();
*/

  return EXIT_SUCCESS;
}
