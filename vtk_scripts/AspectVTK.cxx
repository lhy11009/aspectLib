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
#include <vtkGeometryFilter.h>
#include "AspectVTK.h"


std::map<std::string, std::string> BuildMap(std::ifstream &map_file){
    std::map<std::string, std::string> trans_map;
    std::string key, value;
    while (map_file>>key && getline(map_file, value)){
        if (value.size() > 1)
            trans_map[key] = value[1];  // remove the leading blank
        else 
            throw std::runtime_error("no rule defined for " + key);
    }
    return trans_map;
}

const std::string transform(std::string &s, std::map<std::string, std::string> &m){
    auto map_it = m.find(s);
    if (map_it != m.cend())
        return map_it->second;  // return the substitution.
    else
        return s;  // return the original string if not found
}


void WordTransform(std::ifstream &map_file, std::ifstream &input)
{
    // build map
    auto trans_map = BuildMap(map_file);
    std::string text;
    while (getline(input, text)){
        std::istringstream stream(text);
        std::string word;
        bool first_word = true;
        while (stream >> word)
        {
            if (first_word)
                first_word = false;
            else
                std::cout << ' ';
            std::cout << transform(word, trans_map);
        }
        std::cout << std::endl;
    }
}


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
    vtk_point_data->Update();
    // set value of the private member
    iPolyData = vtkSmartPointer<vtkPolyData>::New();
    iPolyData->SetPoints(vtk_points);  //set as original vector
    std::vector<std::string> array_names;
    array_names.push_back("T");
    array_names.push_back("spcrust");
    array_names.push_back("spharz");
    array_names.push_back("density");
    array_names.push_back("p");
    for (auto p = array_names.begin(); p < array_names.end(); p++)
    {
        if (p == array_names.begin())
            iPolyData->GetPointData()->SetScalars(vtk_point_data->GetArray(p->c_str()));
        else
            iPolyData->GetPointData()->AddArray(vtk_point_data->GetArray(p->c_str()));
    }
    iPolyData->GetPointData()->Update();
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


void AspectVtk::input_poly_data1()
{
    std::cout << "Input poly data: " << std::endl;
    // set value of the private member
    auto geometryFilter = vtkSmartPointer<vtkGeometryFilter>::New();
    geometryFilter->SetInputConnection(reader->GetOutputPort());
    geometryFilter->Update();
    iPolyData = geometryFilter->GetOutput();
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
    std::cout << "Triangulate the grid points, ";
    iDelaunay2D = vtkSmartPointer<vtkDelaunay2D>::New();
    iDelaunay2D->SetInputData(iPolyData);
    iDelaunay2D->Update();
    std::cout << "Done" << std::endl; // debug
}


vtkSmartPointer<vtkPolyData> AspectVtk::extract_contour(const std::string field_name, const double contour_value, const std::string filename)
{
    std::cout << "Filter contour" << std::endl;
    vtkNew<vtkContourFilter> contour_filter;  // simple
    // prepare poly data for contour
    vtkNew<vtkPolyData> cPolyData; 
    cPolyData->DeepCopy(iDelaunay2D->GetOutput());
    vtkSmartPointer<vtkPointData> c_vtk_point_data = iDelaunay2D->GetOutput()->GetPointData();
    cPolyData->GetPointData()->SetScalars(c_vtk_point_data->GetArray(field_name.c_str()));
    cPolyData->GetPointData()->Update();
    // draw contour 
    contour_filter->SetInputData(cPolyData);  // set to the trangulated dataset
    contour_filter->Update();
    contour_filter->GenerateValues(1, contour_value, contour_value);
    contour_filter->Update();
    // write output 
    vtkNew<vtkSimplePointsWriter> writer;
    writer->SetInputData(contour_filter->GetOutput());
    writer->SetFileName(filename.c_str());
    writer->Update();
    writer->Write();
    std::cout << "Generate output (contour): " << filename.c_str() << std::endl;
    // return value
    return contour_filter->GetOutput();
}
        

void AspectVtk::interpolate_uniform_grid(std::string filename)
{
    // Create a grid of points to interpolate over
    std::cout << "interpolate uniform grid" << std::endl;
    std::cout << "\tGenerate new grid" << std::endl;
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

vtkSmartPointer<vtkPolyData> AspectVtk::interpolate_grid(vtkSmartPointer<vtkPoints> gridPoints)
{
    // todo
    // Create a grid of points to interpolate over
    std::cout << "interpolate onto new grid" << std::endl;
    
    // Create a dataset from the grid points
    vtkNew<vtkPolyData> gridPolyData;
    gridPolyData->SetPoints(gridPoints);

    // Perform the interpolation
    vtkNew<vtkProbeFilter> probeFilter;
    probeFilter->SetSourceConnection(iDelaunay2D->GetOutputPort());
    probeFilter->SetInputData(gridPolyData); // Interpolate 'Source' at these points
    probeFilter->Update();
    
    vtkSmartPointer<vtkPolyData> iPolyData = probeFilter->GetPolyDataOutput();
    return iPolyData;
}

void AspectVtk::write_ascii(vtkSmartPointer<vtkPolyData> polyData, std::vector<std::string>& fields, const std::string filename)
{
    // todo
    std::cout << "AspectVtk::write_ascii" << std::endl;
    // output header
    std::string header = "# 1: x (m)\n# 2: y (m)";
    unsigned i = 3;
    for (auto &p: fields)
    {
        header += std::string("\n# ") + std::to_string(i) + std::string(": ") + p;
        i++;
    }
    // output data
    std::string output= "";
    for (vtkIdType i = 0; i < polyData->GetNumberOfPoints(); i++)
    {
        if (i>0)
            output += std::string("\n");  // append new line for new point
        double *xs = polyData->GetPoint(i);
        output += std::to_string(xs[0]) + std::string("\t") + std::to_string(xs[1]);
        for (auto &p: fields){
            const double val = polyData->GetPointData()->GetArray(p.c_str())->GetTuple1(i); 
            output += std::string("\t") + std::to_string(val);
        }
    }
    // write file
    std::ofstream fout(filename, std::ofstream::out);
    if (fout){
        fout << header << std::endl << output;
    }
    std::cout << "\tWrite ascii data file: " << filename << std::endl;
}

void AspectVtk::read_horiz_avg(const std::string &filename)
{
    std::cout << "Read horiz_avg file: " << filename << std::endl;
    file_reader.read_horiz_avg(filename, horiz_depths, horiz_Ts, horiz_densities);
}

double AspectVtk::get_from_horiz(double depth, const std::string &field, const bool fix_out_of_bound)
{
    const unsigned horiz_size = horiz_depths.size();
    if (horiz_size < 1){
        std::cerr << "\tNo data recorded in horiz, call read_horiz_avg first" << std::endl;
    }
    double value = 0.0;  // return value
    if (depth < horiz_depths[0] or depth > horiz_depths[horiz_size - 1])
    {
        if (fix_out_of_bound){
            // fix depth range out of bound
            if (depth < horiz_depths[0]){
                if (field.compare("T") == 0)
                    return horiz_Ts[0];
                else if (field.compare("density") == 0)
                    return horiz_densities[0];
                else{
                    std::cerr << "\tField is not found" << std::endl;
                    exit (1);
                }
            }
            else{
                if (field.compare("T") == 0)
                    return horiz_Ts[horiz_size-1];
                else if (field.compare("density") == 0)
                    return horiz_densities[horiz_size-1];
                else{
                    std::cerr << "\tField is not found" << std::endl;
                    exit (1);
                }
            }
        }
        else
        {
            std::cerr << "\tDepth out of range: " << depth << std::endl;
            exit (1);
        }
    }
    const double interv = horiz_depths[1] - horiz_depths[0];
    const unsigned indice = floor((depth - horiz_depths[0]) / interv);
    const double depth0 = horiz_depths[indice];
    const double depth1 = horiz_depths[indice + 1];
    if (!(depth >= depth0 && depth <= depth1))
    {
        std::cerr << "\tdepth not in the range guessed from the interval from 1st and 2nd entry in the horiz_depths vector, somthing wrong" 
        << std::endl;
        exit (1);
    }
    if (field.compare("T") == 0)
    {
        value = horiz_Ts[indice] + (horiz_Ts[indice+1] - horiz_Ts[indice]) * (depth - depth0) / interv;
    }
    else if (field.compare("density") == 0)
    {
        value = horiz_densities[indice] + (horiz_densities[indice+1] - horiz_densities[indice]) * (depth - depth0) / interv;
    }
    else{
        std::cerr << "\tField is not found" << std::endl;
        exit (1);
    }
    return value;
}
        
void AspectVtk::density_diff()
{
    std::cout << "Pull out differential density" << std::endl;
    const double Ro = 6371e3;
    vtkSmartPointer <vtkPolyData> tpolydata = iDelaunay2D->GetOutput();// polydata after trangulation
    vtkSmartPointer<vtkPoints> vtk_points = reader->GetOutput()->GetPoints(); // coordinates
    auto densities = dynamic_cast<vtkFloatArray*>(tpolydata->GetPointData()->GetArray("density")); // array for density
    vtkNew<vtkFloatArray> densities_diff;
    densities_diff->DeepCopy(densities);
    densities_diff->SetName("density_diff");
    for (vtkIdType i = 0; i < tpolydata->GetNumberOfPoints(); i++)
    {
        double p0[3];
        vtk_points->GetPoint(i, p0);
        double r0 = sqrt(p0[0] * p0[0] + p0[1] * p0[1]);
        const double density = densities->GetTuple1(i);
        const double horiz_density0 = get_from_horiz(Ro - r0, "density", true);
        densities_diff->SetTuple1(i, density - horiz_density0);
    }
    tpolydata->GetPointData()->AddArray(densities_diff);
    iDelaunay2D->Update();
    return;
}


void AspectVtk::mow_from_blocking(const double blocking_T, const double blocking_P)
{
    std::cout << "Pull out mow field" << std::endl;
    const double Ro = 6371e3;
    vtkSmartPointer <vtkPolyData> tpolydata = iDelaunay2D->GetOutput();// polydata after trangulation
    vtkSmartPointer<vtkPoints> vtk_points = reader->GetOutput()->GetPoints(); // coordinates
    auto densities = dynamic_cast<vtkFloatArray*>(tpolydata->GetPointData()->GetArray("density")); // array for density
    auto Ps = dynamic_cast<vtkFloatArray*>(tpolydata->GetPointData()->GetArray("p")); // array for density
    auto Ts = dynamic_cast<vtkFloatArray*>(tpolydata->GetPointData()->GetArray("T")); // array for density
    vtkNew<vtkFloatArray> mow_phases;
    mow_phases->SetName("meta_stable_olivine");
    for (vtkIdType i = 0; i < tpolydata->GetNumberOfPoints(); i++)
    {
        const double Pi = Ps->GetTuple1(i);
        const double Ti = Ts->GetTuple1(i);
        if (Pi > blocking_P && Ti < blocking_T)
            mow_phases->InsertNextValue(1.0);
        else
            mow_phases->InsertNextValue(0.0);
    }
    tpolydata->GetPointData()->AddArray(mow_phases);
    iDelaunay2D->Update();
    return;
}

void AspectVtk::set_geometry(const std::string geometry){
    this->geometry = geometry;
}
    
void AspectVtk::output(const vtkSmartPointer<vtkPolyData> opolydata, const std::string filename){
    std::cout << "Output data" << std::endl;
    // write output 
    vtkNew<vtkXMLPolyDataWriter> writer;
    writer->SetInputData(opolydata);
    writer->SetFileName(filename.c_str());
    writer->Update();
    writer->Write();
    std::cout << "output file: " << filename << std::endl;
}

void FileReader::read_horiz_avg(const std::string &filename,
        std::vector<double> & depths,
        std::vector<double> &Ts,
        std::vector<double> &densities)
{
    std::ifstream file(filename);

    if (file)
    {
        std::stringstream buffer;
        buffer << file.rdbuf();
        file.close();
        std::string temp;
        std::getline(buffer, temp); // get next line, file header
        std::cout << "\tread header" << std::endl;
        std::vector<int> indices(3, -1); // column indices
        int n_columns, n_rows;
        unsigned i = 0;
        while (temp.compare(0, 1, "#") == 0)
        {
            if (temp.size() > 5)
                if (temp.compare(temp.size()-5,5,"depth") == 0){
                    std::cout << "\t\tdepth matched: column" << i << std::endl;
                    indices[0] = i;
                }
            if (temp.size() > 11)
                if (temp.compare(temp.size()-11,11,"temperature") == 0){
                    std::cout << "\t\ttemperature matched: column" << i << std::endl;
                    indices[1] = i;
                }
            if (temp.size() > 17)
                if (temp.compare(temp.size()-17,17,"adiabatic_density") == 0){
                    std::cout << "\t\tadiabatic_density matched: column" << i << std::endl;
                    indices[2] = i;
                }
            std::getline(buffer, temp); 
            i++;
        }
        for (auto &p: indices) // check for header indices
            if(p < 0){
                std::cerr << "\tread_horiz_avg: missing some fields" << std::endl;
                exit (1);
            }
        std::istringstream iss(temp);  // read shape of data
        iss >> n_rows >> n_columns;
        std::cout << "\tread the shape of data: " << n_rows << " " << n_columns << std::endl;
        std::vector<double> row_values(n_columns);
        depths.resize(n_rows, 0.0);
        Ts.resize(n_rows, 0.0);
        densities.resize(n_rows, 0.0);
        i = 0;
        while(!buffer.eof()) // read data
        {
            for (unsigned j=0; j<n_columns; j++)
                buffer >> row_values[j];
            depths[i] = row_values[indices[0]];
            Ts[i] = row_values[indices[1]];
            densities[i] = row_values[indices[2]];
            i++;
        }
    }
    else
    {
        std::cerr << "File: " << filename << " cannot be opened" << std::endl;
    }
}

void read_options(const std::string &option_filename, std::vector<std::string> &options)
{
    std::ifstream in(option_filename);
    unsigned i;
    if (in.is_open())
    {
        std::string temp;
        while (in>>temp)
        {
            if (in.fail())
            {
                in.clear();
            }
            options.push_back(temp);
            i++;
        }
        in.close();
    }
    else std::cerr << "Unable to open file: " << option_filename << std::endl;
}
