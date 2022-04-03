#include <math.h>
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include "AspectVTK.h"

void QuickSort(const std::vector<double> & A, std::vector<size_t> & I, size_t lo, size_t hi);


struct SlabOutputs
{
    double trench_theta; // trench position
    double slab_depth;  // slab depth
    double theta100; // angle where the slab is 100 km depth
    double dip100;
    double total_buoyancy; // total buoyancy forces on the slab
    double m_poly_r_min;
    double m_poly_r_max;
    double slab_profile_interval; //interval in r with the slab profiles
    double slab_shallow_cutoff; // a shallow cutoff for the slab profiles
    unsigned m_poly_total;
    /*this is radius (or y) vs theta (or x) following the composition
    contour of crust + harz*/
    std::pair <std::shared_ptr<std::vector<double>>, std::shared_ptr<std::vector<double>>> composition_contour;
    std::shared_ptr<std::vector<size_t>> composition_contour_sorted_index;  // save the indexes of sorted profile so that radius increases
    std::shared_ptr<std::vector<size_t>> internal_cell_index;  // save the indexes of internal cells
    std::shared_ptr<std::vector<size_t>> slab_envelop_left_indexes;  // save the indexes of the envelop cells
    std::shared_ptr<std::vector<size_t>> slab_envelop_right_indexes;  // save the indexes of the envelop cells
    std::shared_ptr<std::vector<double>> in_slab_depths; // depths within the range of [slab_depth_shallow, slab_depth]
    std::shared_ptr<std::vector<double>> in_slab_buoyancies; // depths within the range of [slab_depth_shallow, slab_depth]
    /*this is radius (or y) vs theta (or x) following the slab surface*/
    std::pair <std::vector<double>, std::vector<double>> surface_profile;
    vtkSmartPointer<vtkPolyData> m_poly_data;  // reference field of mantle 
    SlabOutputs(){
      m_poly_r_min = 6371e3 - 2890e3 + 1e3;
      m_poly_r_max = 6371e3 - 1e3;
      m_poly_total = 3000;
      slab_profile_interval = 5e3;
      slab_shallow_cutoff = 50e3;
      this->composition_contour_sorted_index = std::make_shared<std::vector<size_t>>();
      this->internal_cell_index = std::make_shared<std::vector<size_t>>();
      this->in_slab_depths = std::make_shared<std::vector<double>>();
      this->slab_envelop_left_indexes = std::make_shared<std::vector<size_t>>();
      this->slab_envelop_right_indexes = std::make_shared<std::vector<size_t>>();
    }
    void export_profile(std::string filename, const std::string geometry);
    void export_forces_profile(std::string filename, const std::string geometry);
    double get_mantle_field(const double r, const std::string field);
};

void SlabOutputs::export_profile(std::string filename, const std::string geometry){
    std::string header = "# 1: x (m)\n# 2: y (m)";
    std::string output;
    bool is_first = true;
    for (unsigned i = 0; i<this->surface_profile.first.size(); i++){
        double theta = this->surface_profile.first[i];
        double r = this->surface_profile.second[i];
        double x, y;
        if (geometry == "chunk"){
          x = r * cos(theta); // convert to x, y
          y = r * sin(theta);
        }
        else if (geometry == "box"){
          x = theta;
          y = r;
        }
        if (is_first)
        {
          is_first = false;
        }
        else
        {
          output += std::string("\n");
        }
        output += std::to_string(x) + std::string("\t") + std::to_string(y);
    }
    // write file
    std::ofstream fout(filename, std::ofstream::out);
    if (fout){
        fout << header << std::endl << output << std::endl;
    }
    std::cout << "SlabOutputs::export_profile: Write ascii data file: " << filename << std::endl;
}


void SlabOutputs::export_forces_profile(std::string filename, const std::string geometry){
    std::string header = "# 1: depth (m)\n# 2: buoyancy (N/m)";
    std::string output;
    bool is_first = true;
    for (unsigned i = 0; i<this->in_slab_depths->size()-1; i++){
        double depth = ((*this->in_slab_depths)[i] + (*this->in_slab_depths)[i+1]) / 2.0;
        double buoyancy = (*this->in_slab_buoyancies)[i];
        if (is_first)
        {
          is_first = false;
        }
        else
        {
          output += std::string("\n");
        }
        output += std::to_string(depth) + std::string("\t") + std::to_string(buoyancy);
    }
    // write file
    std::ofstream fout(filename, std::ofstream::out);
    if (fout){
        fout << header << std::endl << output << std::endl;
    }
    std::cout << "SlabOutputs::export_profile: Write ascii data file: " << filename << std::endl;
}


double SlabOutputs::get_mantle_field(const double r, const std::string field)
{
    const double foo = (r - this->m_poly_r_min) / 
                       (this->m_poly_r_max - this->m_poly_r_min)*
                        this->m_poly_total;
    const unsigned i0 = floor(foo);
    const double residual = foo - i0;
    const double r0 = this->m_poly_r_min +
         i0 * (this->m_poly_r_max - this->m_poly_r_min) / this->m_poly_total ;
    const double r1 = this->m_poly_r_min +
         (i0 + 1) * (this->m_poly_r_max - this->m_poly_r_min) / this->m_poly_total ;
    const double val0 = this->m_poly_data->GetPointData()->GetArray(field.c_str())->GetTuple1(i0);
    const double val1 = this->m_poly_data->GetPointData()->GetArray(field.c_str())->GetTuple1(i0+1);
    double v = val0 * (1-residual) + val1 * residual;
    return v;
}



class SlabAnalysis : public AspectVtk
{
    public:
        void prepare_slab(const std::vector<std::string> names);
        void construct_slab_internal(SlabOutputs & slab_outputs, const double slab_threshold, const double T_threshold);
        void export_slab_points(std::string filename, const std::shared_ptr <std::vector<size_t>> Iic);
        void slab_buoyancy(SlabOutputs &slab_outputs);
        void interpolate_mantle_profile(SlabOutputs &slab_outputs);
        void export_mantle_reference_profile(const SlabOutputs &slab_outputs, std::string filename);
        void slab_envelop(SlabOutputs &slab_outputs);
};


void SlabAnalysis::export_mantle_reference_profile(const SlabOutputs &slab_outputs, std::string filename)
{
    std::vector<std::string> fields = {"density"};
    write_ascii(slab_outputs.m_poly_data, fields, filename);
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

void SlabAnalysis::construct_slab_internal(SlabOutputs & slab_outputs, const double slab_threshold, const double T_threshold)
{
    const double ro = 6371e3;
    const double slab_thickness_cutoff = 0.0; // used for cutoff while looking for internal points
    vtkSmartPointer <vtkPolyData> tpolydata = iDelaunay2D->GetOutput();// polydata after trangulation
    std::cout << "Construct slab internal" << std::endl;
    std::cout << "\tcell type: " << tpolydata->GetCell(0)->GetCellType()
        << ", number of cells: " << tpolydata->GetNumberOfCells() << std::endl;  // check cell type
    // const double rin = 2890e3;  // rin
    auto slab = dynamic_cast<vtkFloatArray*>(tpolydata->GetPointData()->GetArray("slab")); // array for slab
    auto Ts = dynamic_cast<vtkFloatArray*>(tpolydata->GetPointData()->GetArray("T")); // array for temperature
    for (vtkIdType i = 0; i < tpolydata->GetNumberOfCells(); i++)
    {
        vtkCell* cell = tpolydata->GetCell(i);  // cell in this polydata
        vtkIdList* idList = cell->GetPointIds();  // list of point ids in this cell
        // point and area
        double p0[3], p1[3], p2[3];
        cell->GetPoints()->GetPoint(0, p0);
        cell->GetPoints()->GetPoint(1, p1);
        cell->GetPoints()->GetPoint(2, p2);
        const double x = (p0[0] + p1[0] + p2[0]) / 3.0;  // x of cell center
        const double y = (p0[1] + p1[1] + p2[1]) / 3.0;  // y of cell center
        const double r = sqrt(x * x + y * y);
        const double theta = acos(x/r);
        if ( (ro - r) < (slab_outputs.slab_shallow_cutoff) ||
             (ro - r) > (slab_outputs.slab_depth + slab_thickness_cutoff))
        {
          continue;  // apply a radius cutoff for points that are too shallow or too deep
        }
        /* subducting slab composition */
        const double vslab = (slab->GetTuple1(idList->GetId(0)) +
                      slab->GetTuple1(idList->GetId(1)) +
                      slab->GetTuple1(idList->GetId(2))) / 3.0 ;
        if (vslab > slab_threshold){
          /*add points that has slab composition*/
          slab_outputs.internal_cell_index->push_back(i); // add this point to the slab internal
          continue;
        }
        /* temperature*/
        const double T = (Ts->GetTuple1(idList->GetId(0)) +
                      Ts->GetTuple1(idList->GetId(1)) +
                      Ts->GetTuple1(idList->GetId(2))) / 3.0 ;
        unsigned intercept = 0;
        double phi = 0.0;
        double phi_last =0.0;
        double theta_j, r_j, x_j, y_j;
        for (unsigned j; j< slab_outputs.surface_profile.first.size(); j++)
        {
          theta_j = slab_outputs.surface_profile.first[j];
          r_j = slab_outputs.surface_profile.second[j];
          x_j = r_j * cos(theta_j);
          y_j = r_j * sin(theta_j);
          if (r_j < r){
            break; // only look at points about cell
          }
          phi = atanh((x_j - x)/(y_j - y));  // angle from vertical up
          if (j > 0 && phi * phi_last < 0.0) {
            intercept++;
          }
          phi_last = phi;
        }
        if (intercept % 2 == 1 and T < T_threshold){
          slab_outputs.internal_cell_index->push_back(i);// add this point to the slab internal
        }
    }
    std::cout << "\tNumber of the slab internal cells: " << slab_outputs.internal_cell_index->size() << std::endl;
}


void SlabAnalysis::export_slab_points(std::string filename, const std::shared_ptr<std::vector<size_t>> Iic){
    std::cout <<  "SlabAnalysis::export_slab_internal" << std::endl;
    vtkSmartPointer <vtkPolyData> tpolydata = iDelaunay2D->GetOutput();// polydata after trangulation
    std::string header = "# 1: x (m)\n# 2: y (m)";
    std::string output = "";
    bool is_first = true;
    for (auto &i: *Iic){
        vtkCell* cell = tpolydata->GetCell(i);  // cell in this polydata
        vtkIdList* idList = cell->GetPointIds();  // list of point ids in this cell
        // point and area
        double p0[3], p1[3], p2[3];
        cell->GetPoints()->GetPoint(0, p0);
        cell->GetPoints()->GetPoint(1, p1);
        cell->GetPoints()->GetPoint(2, p2);
        const double x = (p0[0] + p1[0] + p2[0]) / 3.0;  // x of cell center
        const double y = (p0[1] + p1[1] + p2[1]) / 3.0;  // y of cell center
        if (is_first)
        {
          is_first = false;
        }
        else
        {
          output += std::string("\n");
        }
        output += std::to_string(x) + std::string("\t") + std::to_string(y);
    }
    // write file
    std::ofstream fout(filename, std::ofstream::out);
    if (fout){
        fout << header << std::endl << output << std::endl;
    }
    std::cout << "\tWrite ascii data file: " << filename << std::endl;
}


//todo
void SlabAnalysis::slab_buoyancy(SlabOutputs &slab_outputs)
{
    const double ro = 6371e3;
    std::cout << "SlabAnalysis:slab_buoyancy" << std::endl;
    const double grav_acc = 10.0;  // value for gravity acceleration
    std::shared_ptr<std::vector<size_t>> Iic = slab_outputs.internal_cell_index;  // get inputs from slab_outputs
    const double slab_depth = slab_outputs.slab_depth;
    vtkSmartPointer <vtkPolyData> tpolydata = iDelaunay2D->GetOutput();// polydata after trangulation
    std::cout << "\tcell type: " << tpolydata->GetCell(0)->GetCellType()
        << ", number of cells: " << tpolydata->GetNumberOfCells() << std::endl;  // check cell type
    // const double rin = 2890e3;  // rin
    auto slab = dynamic_cast<vtkFloatArray*>(tpolydata->GetPointData()->GetArray("slab")); // array for spcrust
    auto density_data = dynamic_cast<vtkFloatArray*>(tpolydata->GetPointData()->GetArray("density")); // array for density
    double total_area = 0.0;
    double total_gravity = 0.0;
    double total_buoyancy = 0.0;
    int out_lyers = 0;
    /* initiate the in_slab_buoyancy with a bunch of 0.0*/
    slab_outputs.in_slab_buoyancies = std::make_shared<std::vector<double>>(
                                    std::vector<double> (slab_outputs.in_slab_depths->size()-1,0.0)
                                    );
    for (auto &i: *Iic)
    {
        /* loop for slab internal cells*/
        vtkCell* cell = tpolydata->GetCell(i);  // cell in this polydata
        vtkIdList* idList = cell->GetPointIds();
        // point and area
        double p0[3], p1[3], p2[3], pc[3];
        cell->GetPoints()->GetPoint(0, p0);  // get nodal points
        cell->GetPoints()->GetPoint(1, p1);
        cell->GetPoints()->GetPoint(2, p2);
        double r0 = sqrt(p0[0] * p0[0] + p0[1] * p0[1]);
        double r1 = sqrt(p1[0] * p1[0] + p1[1] * p1[1]);
        double r2 = sqrt(p2[0] * p2[0] + p2[1] * p2[1]);
        const double x = (p0[0] + p1[0] + p2[0])/3.0;
        const double y = (p0[1] + p1[1] + p2[1])/3.0;
        const double r = sqrt(x * x + y * y);
        double area = vtkTriangle::TriangleArea(p0, p1, p2);
        // read from horiz_avg
        // const double horiz_density0 = get_from_horiz(Ro - r0, "density", true);
        // const double horiz_density1 = get_from_horiz(Ro - r1, "density", true);
        // const double horiz_density2 = get_from_horiz(Ro - r2, "density", true);
        // subducting slab
        double slab0 = slab->GetTuple1(idList->GetId(0));
        double slab1 = slab->GetTuple1(idList->GetId(1));
        double slab2 = slab->GetTuple1(idList->GetId(2));
        // density
        double density0 = density_data->GetTuple1(idList->GetId(0));
        double density1 = density_data->GetTuple1(idList->GetId(1));
        double density2 = density_data->GetTuple1(idList->GetId(2));
        // average value, future: volume from trangulation
        double slab_avg = (slab0 + slab1 + slab2) / 3.0;
        double density_avg = (density0 + density1 + density2) / 3.0;
        //double horiz_density_avg = (horiz_density0 + horiz_density1 + horiz_density2) / 3.0;
        double horiz_density_avg = slab_outputs.get_mantle_field(r, "density");
        // areas on composition
        double gravity = density_avg * grav_acc * area;
        double buoyancy = -(density_avg - horiz_density_avg) * grav_acc * area;
        // add up
        total_area += area;
        total_gravity += gravity;
        total_buoyancy += buoyancy;
        unsigned i_s = floor((ro - r - slab_outputs.slab_shallow_cutoff) / (slab_outputs.slab_profile_interval));
        if (i_s < slab_outputs.in_slab_depths->size()-1)
        {
          (*slab_outputs.in_slab_buoyancies)[i_s] += buoyancy;
        }
        else
        {
          out_lyers += 1;
        }
    }
    slab_outputs.total_buoyancy = total_buoyancy;
    std::cout<< "\tGet " << slab_outputs.internal_cell_index->size() <<
    " cells, with " << out_lyers <<" outlyers" << std::endl;
    std::cout << "total slab area (m^2): " << total_area << std::endl <<
    "total gravity (N/m): " << total_gravity << std::endl <<
    "total buoyancy (N/m): " << total_buoyancy <<std::endl;
}


void SlabAnalysis::interpolate_mantle_profile(SlabOutputs &slab_outputs)
{
    std::cout << "SlabAnalysis::interpolate_mantle_profile" << std::endl;
    vtkNew<vtkPoints> gridPoints;
    const double ro = slab_outputs.m_poly_r_max;
    const double ri = slab_outputs.m_poly_r_min;
    unsigned total = slab_outputs.m_poly_total;
    const double theta0 = 0.0017; // 1 deg
    std::vector<double> rs;
    for (unsigned i = 0; i < total; i++)
    {
      rs.push_back(ri + (ro-ri)/(total-1) * i);
    }
    for (unsigned int i = 0; i < total; i++)
    {
      double r = rs[i];
      double x = r * cos(theta0);
      double y = r * sin(theta0);
      gridPoints->InsertNextPoint(x, y, 0);
    }
    slab_outputs.m_poly_data = interpolate_grid(gridPoints);
}


void SlabAnalysis::slab_envelop(SlabOutputs &slab_outputs)
{
    /* devide the whole envelop into this segment*/
    const double ro = 6371e3;
    std::cout << "SlabAnalysis::slab_envelop" << std::endl;
    const unsigned n_segment = ceil((slab_outputs.slab_depth-slab_outputs.slab_shallow_cutoff) / slab_outputs.slab_profile_interval);
    std::vector<std::vector<unsigned>> ids (n_segment, std::vector<unsigned>());
    std::vector<std::vector<double>> thetas (n_segment, std::vector<double>());
    vtkSmartPointer <vtkPolyData> tpolydata = iDelaunay2D->GetOutput();// polydata after trangulation
    unsigned outlyer = 0;
    unsigned count = 0;
    std::cout << "\treorder the internal points" << std::endl;
    for (auto &i: *(slab_outputs.internal_cell_index))
    {
        /* First loop for internal points and reorganize them in terms of radius
        Note here i is the index of the cells.*/
        vtkCell* cell = tpolydata->GetCell(i);  // cell in this polydata
        vtkIdList* idList = cell->GetPointIds();  // list of point ids in this cell
        // point and area
        double p0[3], p1[3], p2[3];
        cell->GetPoints()->GetPoint(0, p0);
        cell->GetPoints()->GetPoint(1, p1);
        cell->GetPoints()->GetPoint(2, p2);
        const double x = (p0[0] + p1[0] + p2[0]) / 3.0;  // x of cell center
        const double y = (p0[1] + p1[1] + p2[1]) / 3.0;  // y of cell center
        const double r = sqrt(x * x + y * y);
        const double theta = acos(x/r);
        /* note the usage of floor and ceil here, in this way the max value of the index is n_segment - 1*/
        unsigned i_r = floor((ro - r - slab_outputs.slab_shallow_cutoff) / slab_outputs.slab_profile_interval);
        if (i_r <= n_segment-1)
        {
            thetas[i_r].push_back(theta); // record the theta of the cell in the segment
            ids[i_r].push_back(i); // record the index of the cell in the segment
        }
        else{
            outlyer ++;
        }
        count ++;
    }
    std::cout << "\tCount: " << count << ", Outlyer: " << outlyer << std::endl;
    std::cout << "\tderive the envelop" << std::endl;
    for (unsigned i_r = 0; i_r < n_segment; i_r++)
    {
        /* Then extract boundary points within each radius groups.
        Note here i_r is the index of segments*/
        auto &cell_ids = ids[i_r]; // ids of cell in this segment
        if (cell_ids.size()==0){
            continue;  //if no points is present, look into the next segment
        }
        auto I = std::make_shared<std::vector<size_t>>();
        for (unsigned id = 0; id < cell_ids.size(); id++)
        {
            I->push_back(id);
        }
        auto &theta = thetas[i_r];
        QuickSort(theta, *I, 0, I->size()-1);  // reorder by r, index saved in I
        size_t I_min = I->front();
        size_t I_max = I->back();
        slab_outputs.slab_envelop_left_indexes->push_back(cell_ids[I_min]);  // record the index of cell on the left
        slab_outputs.slab_envelop_right_indexes->push_back(cell_ids[I_max]); // and right of this segment.
    }

}


std::shared_ptr<std::vector<size_t>> find_points_at_radius(std::shared_ptr<std::vector<double>> rs, const double r,
                                                           const unsigned n)
{
  // found trench point, this would around the indexes of the point according to theire nearness to the radius provided
  auto r_deviation = std::make_shared<std::vector<double>>();
  auto I_deviation = std::make_shared<std::vector<size_t>>();
  for (size_t id = 0; id < rs->size(); id++)
  {
    r_deviation->push_back(abs((*rs)[id] - r));
    I_deviation->push_back(id);
  }
  QuickSort(*r_deviation, *I_deviation, 0, rs->size()-1);
  return I_deviation;
}

double find_slab_interal_at_depth(std::shared_ptr<std::vector<double>> rs, std::shared_ptr<std::vector<double>> thetas,
                                        const double theta_ref, const double radius_find_point,
                                        const int n_average)
{
  // This function would find a point at the slab internal by value of theta (or x) with a given profile and a give depth.
  // Note rs and thetas represent a contour the envelop the slab composition (i.e crust + harzburgite)
  // And we are looking at the points around a certain depth and then average their theta (or x) value.
  // By presenting an initial guess, we would save us from selecting points from the far end of the domain.
  auto Ip = find_points_at_radius(rs, radius_find_point, n_average);
  double sum = 0.0;
  int found = 0;
  for (int i=0; i < Ip->size(); i++){
    double temp = ((*thetas)[(*Ip)[i]]- theta_ref) / theta_ref;  // nearness relative to an initial guess
    if (abs(temp) < 0.2) // prevent selecting points from the far end of the domain
    {
      sum += (*thetas)[(*Ip)[i]];
      found++;
    }
    if (found == n_average)
      break;
  }
  if (found < n_average)
  {
    std::cerr << "find_slab_interal_at_depth: not enough points found for slab internal point at radius " << radius_find_point  << std::endl; // error message
    exit(1);
  }
  double theta = sum / n_average;
  return theta;
}


double find_slab_surface_at_depth(std::shared_ptr<std::vector<double>> rs, std::shared_ptr<std::vector<double>> thetas,
                                        const double theta_ref, const double radius_find_point,
                                        const int n_average)
{
  // This function would find a point at the slab surface by value of theta (or x) with a given profile and a give depth.
  // Note rs and thetas represent a contour the envelop the slab composition (i.e crust + harzburgite)
  // And we are looking at the points around a certain depth and find out the one with maximum theta (or x)
  // By presenting an initial guess, we would save us from selecting points from the far end of the domain.
  auto Ip = find_points_at_radius(rs, radius_find_point, n_average);
  double theta_surface = -1.0;
  bool found = false;
  for (int i=0; i < n_average; i++){
    // get the biggest value for a point on the surface
    double theta = (*thetas)[(*Ip)[i]];
    double temp = (theta- theta_ref) / theta_ref;  // nearness relative to an initial guess
    if (theta > theta_surface && abs(temp) < 0.2){
      theta_surface = (*thetas)[(*Ip)[i]];
      found = true;
    }
  }
  if (!found)
  {
    std::cerr << "find_slab_surface_at_depth_by_guess: no point found at radius " << radius_find_point << std::endl; // error message
    exit(1);
  }
  return theta_surface;
}


void reorder_contour_points(vtkSmartPointer<vtkPolyData> c_poly_data, SlabOutputs & slab_outputs, const std::string geometry="chunk")
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
    if (geometry == "chunk"){
      rs->push_back(r);
      thetas->push_back(theta);
    }
    else if (geometry == "box"){
      rs->push_back(y);
      thetas->push_back(x);
    }
    else
      throw std::runtime_error("Geometry must be either \"box\" or \"chunk\"");
    I->push_back(id);
  }
  QuickSort(*rs, *I, 0, c_points->GetNumberOfPoints()-1);  // reorder by r, index saved in I
  slab_outputs.composition_contour = make_pair(thetas, rs);  // save this to composition_contour
  slab_outputs.composition_contour_sorted_index = I; 
}


void analyze_slab(SlabOutputs & slab_outputs, const double theta_ref_trench, const std::string geometry="chunk")
{
  // parameters we use here
  double ro;
  if (geometry == "chunk"){
    ro = 6371e3;
  }
  else if (geometry == "box"){
    ro = 2890e3;
  }
  double depth_find_trench_point = 5e3; // parameters to find trench point
  const int n_average_for_trench_point = 3;
  double depth_find_100km = 100e3; // parameters to find slab surface at 100km depth
  const int n_average_for_100_km = 10;
  // found trench point
  auto thetas = slab_outputs.composition_contour.first;  // get profile of compositional contour from slab_outputs
  auto rs = slab_outputs.composition_contour.second;
  auto I = slab_outputs.composition_contour_sorted_index;
  const double trench_theta = find_slab_interal_at_depth(rs, thetas, theta_ref_trench, ro-depth_find_trench_point, n_average_for_trench_point);
  slab_outputs.trench_theta = trench_theta;
  std::cout << "Trench theta: " << trench_theta << std::endl;
  // find slab depth
  double slab_depth;
  slab_depth = ro - (*rs)[(*I)[0]];
  slab_outputs.slab_depth = slab_depth;
  double depth = slab_outputs.slab_shallow_cutoff;  // construct the array of in_slab_depths
  while(depth <= slab_depth){
    slab_outputs.in_slab_depths->push_back(depth);
    depth += slab_outputs.slab_profile_interval;
  }
  std::cout << "Slab depth: " << slab_depth << std::endl;
  // find points at 100e3 depth
  double theta100;
  theta100 = find_slab_surface_at_depth(rs, thetas, theta_ref_trench, ro-depth_find_100km, n_average_for_100_km);
  slab_outputs.theta100 = theta100;
  if (geometry == "chunk")
    slab_outputs.dip100 = atan(ro*(theta100 - trench_theta)/depth_find_100km);
  else if (geometry == "box")
    slab_outputs.dip100 = atan((theta100 - trench_theta)/depth_find_100km);
  std::cout << "100km theta: " << theta100 << std::endl;
  std::cout << "100km dip: " << slab_outputs.dip100 << std::endl;
}


void analyze_slab_surface_profile(SlabOutputs & slab_outputs, const std::string geometry="chunk")
{
  // parameters we use here
  double ro;
  if (geometry == "chunk"){
    ro = 6371e3;
  }
  else if (geometry == "box"){
    ro = 2890e3;
  }
  const int n_average = 10;
  // slab profile
  double profile_theta = slab_outputs.trench_theta; // use trench position as the initial guess
  const double profile_interval= 1e3;  // interval in profile points
  double profile_depth_start = 100e3;  // start point of profile
  auto thetas = slab_outputs.composition_contour.first;  // get profile of compositional contour from slab_outputs
  auto rs = slab_outputs.composition_contour.second;
  for (double profile_depth = profile_depth_start; profile_depth <= slab_outputs.slab_depth; profile_depth += profile_interval)
  {
    slab_outputs.surface_profile.second.push_back(ro-profile_depth);
    profile_theta = find_slab_surface_at_depth(rs, thetas, profile_theta, ro-profile_depth, n_average);
    slab_outputs.surface_profile.first.push_back(profile_theta);
  }
}




void analyze_wedge_temperature100(SlabAnalysis & slab_analysis, SlabOutputs & slab_outputs, const std::string filename, const std::string geometry="chunk")
{
  // wedge temperature above where the slab is 100 km deep
  vtkNew<vtkPoints> gridPoints;
  double ro;
  if (geometry == "chunk"){
    ro = 6371e3;
  }
  else if (geometry == "box"){
    ro = 2890e3;
  }
  double depth = 120e3;  // this is set > 100 so as to check for slab compositions
  unsigned num = 120;
  double val_theta, val_r, val_x, val_y;
  for (unsigned i = 0; i < num; i++){
    // a vertical line
    if (geometry == "chunk"){
      val_theta = slab_outputs.theta100;
      val_r = ro - i * 1.0 / num * (depth);
      val_x = val_r * cos(val_theta);
      val_y = val_r * sin(val_theta);
    }
    else if (geometry == "box"){
      val_x = slab_outputs.theta100;
      val_y = ro - i * 1.0 / num * (depth);
    }
    gridPoints->InsertNextPoint(val_x, val_y, 0);
  }
  vtkSmartPointer<vtkPolyData> iPolyData = slab_analysis.interpolate_grid(gridPoints);
  std::vector<std::string> fields{"T", "slab"};
  slab_analysis.write_ascii(iPolyData, fields, filename);
}


void QuickSort(const std::vector<double> & A, std::vector<size_t> & I, size_t lo, size_t hi)
{
    /* QuickSort by reordering the indexes, by an increasing order*/
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
  /* Parse options
  1. PVTU_FILE
  */
  std::string filename = options[0]; // get dir
  size_t found=option_file.find_last_of("/\\");
  std::string target_dir = option_file.substr(0,found);  // upper lever route, in this way, one need to provide the full path here.
  found=filename.find_last_of("."); // get step
  size_t found1 = filename.find_last_of("-");
  std::string pvtu_step = filename.substr(found1+1, found - found1-1);
  std::cout << "pvtu_step: " << pvtu_step << std::endl;
  /*2. VTK_HORIZ_FILE */
  std::string avg_filename = options[1];
  /* 3. THETA_REF_TRENCH
   if there is entry (at least 4 in the file), assign it to theta_ref_trench
   otherwise, assign a default value*/
  double theta_ref_trench = 0.63;
  if (options.size() > 2)
    theta_ref_trench = std::stod(options[2]);
  /* 4. geometry
    if there is an entry for geometry (at least 4 in total), assign it to geometry
    otherwise, assign 'chunk'*/
  std::string geometry;
  if (options.size() > 3)
    geometry = options[3];
  else
    geometry = "chunk";
  assert((geometry == "box") || (geometry == "chunk"));
  /* 4. geometry
    if there is an entry for geometry (at least 4 in total), assign it to geometry
    otherwise, assign 'chunk'*/
  std::string operation;
  if (options.size() > 4)
    operation = options[4];
  else
    operation = "morphology";
  assert((operation == "morphology") || (operation == "default"));
  // Initialiation & read option & data
  SlabAnalysis slab_analysis;
  slab_analysis.set_geometry(geometry);
  slab_analysis.readfile(filename);
  slab_analysis.read_horiz_avg(avg_filename);
  slab_analysis.input_poly_data();
  slab_analysis.triangulate_grid();  # implementation
  //slab_analysis.density_diff();
  //slab_analysis.mow_from_blocking(973.0, 12.5e9);  // 725 + 273 from Quinteros
  slab_analysis.prepare_slab({"spcrust", "spharz"}); // here we create a new field called "slab" from adding "spcrust" and "spharz".
  /* Here, the extract_contour function extract the contour of "slab" based on the second value give and output it to the path with the
  third value given.
  */
  vtkSmartPointer<vtkPolyData> c_poly_data = slab_analysis.extract_contour("slab", 0.99, target_dir + "/" + "contour_slab_" + pvtu_step + ".txt"); //apply contour
  SlabOutputs slab_outputs;
  reorder_contour_points(c_poly_data, slab_outputs, geometry); // sort the points on contour by radius (small to big)
  if (operation == "default"){
    /* by default, handle both morphology and buoyancy*/
    analyze_slab(slab_outputs, theta_ref_trench, geometry); // analyze slab morphology
    analyze_wedge_temperature100(slab_analysis, slab_outputs, target_dir + "/" + "wedge_T100_" + pvtu_step + ".txt", geometry); // output tempertature
    analyze_slab_surface_profile(slab_outputs, geometry);
    slab_outputs.export_profile(target_dir + "/" + "slab_surface_" + pvtu_step + ".txt", geometry);  // export profile of slab surface
    slab_analysis.construct_slab_internal(slab_outputs, 0.8, 1473.0);  // find indexes of cell for internal cells
    slab_analysis.export_slab_points(target_dir + "/" + "slab_internal_" + pvtu_step + ".txt",\
                                       slab_outputs.internal_cell_index); // export slab internal point
    slab_analysis.interpolate_mantle_profile(slab_outputs); // get a mantle reference profile
    slab_analysis.export_mantle_reference_profile(slab_outputs, target_dir + "/" + "mantle_reference_" + pvtu_step + ".txt");  // export the mantle profile
    slab_analysis.slab_buoyancy(slab_outputs); // analyze slab buoyancies
    slab_outputs.export_forces_profile(target_dir + "/" + "slab_forces_" + pvtu_step + ".txt", geometry);
    slab_analysis.slab_envelop(slab_outputs);
    slab_analysis.export_slab_points(target_dir + "/" + "slab_envelop_left_" + pvtu_step + ".txt",\
                                       slab_outputs.slab_envelop_left_indexes); // export slab envelop points
    slab_analysis.export_slab_points(target_dir + "/" + "slab_envelop_right_" + pvtu_step + ".txt",\
                                       slab_outputs.slab_envelop_right_indexes); // export slab envelop points
  }
  else if (operation == "morphology"){
    /* by "morphology", only look at the morphology of the slab*/
    analyze_slab(slab_outputs, theta_ref_trench, geometry); // analyze slab morphology
    analyze_wedge_temperature100(slab_analysis, slab_outputs, target_dir + "/" + "wedge_T100_" + pvtu_step + ".txt", geometry); // output tempertature
  }
  //aspect_vtk.interpolate_uniform_grid("uniform2D.vtp");  // intepolation
  return EXIT_SUCCESS;
}