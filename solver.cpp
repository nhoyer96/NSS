// Solver for the TOV equation
// Solves for a static and radially symmetric Neutron Star
//
// Written by
// Nils Hoyer           -- <n.hoyer@stud.uni-heidelberg.de>
// Maximilian Sasserath -- <sasserath@stud.uni-heidelberg.de>

// Libraries
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector> // NEW
#include <sstream> // NEW
#include <iomanip> // NEW

// Class to read in data table..
class CSVRow{
  public:
    std::string const & operator [](std::size_t index) const{
      return m_data[index];
    }
    std::size_t size() const{
      return m_data.size();
    }
    void ReadNextRow(std::istream & str){
      std::string line, cell;
      std::getline(str, line);
      std::stringstream linestream(line);

      m_data.clear();
      while(std::getline(linestream, cell, ',')){
        m_data.push_back(cell);
      }
      if(!linestream && cell.empty()) m_data.push_back("");
    }
  private:
    std::vector<std::string> m_data;
}row;

std::istream & operator >>(std::istream & str, CSVRow & data){
  data.ReadNextRow(str);
  return str;
}

// Function declarations
inline double GetDensity(double & pressure, double & kappa, double & gamma);
inline double GetPressure(double & density, double & kappa, double & gamma);
inline double GetPressure(const double & density, double & kappa, double & gamma);
double TOV(double radius, double & mass, double pressure, const double & c, const double & G, const double & pi, double & kappa, double & gamma);

double Euler(double & radius, double & mass, double & pressure, double & dr, const double & c, const double & G, const double & pi, double & kappa, double & gamma);
double Rk2(double & radius, double & mass, double & pressure, double & dr, const double & c, const double & G, const double & pi, double & kappa, double & gamma);
double Rk3(double & radius, double & mass, double & pressure, double & dr, const double & c, const double & G, const double & pi, double & kappa, double & gamma);
double Rk4(double & radius, double & mass, double & pressure, double & dr, const double & c, const double & G, const double & pi, double & kappa, double & gamma);

void Adapt(std::string & file_name, std::string & method, double & initial_guess, double & maximal_radius, double & total_mass, double & frac, double & dr, double & stop, double & dr_div, unsigned long long max_steps, double & mass_resolution, double & lower_boundary, double & upper_boundary, double & kappa, double & gamma, const double & c, const double & G, const double & pi);
std::vector < std::vector < double > > Single(std::string & file_name, std::string & method, double & initial_guess, double & maximal_radius, double & total_mass, double & frac, double & dr, double & stop, double & dr_div, unsigned long long max_steps, double & mass_resolution, double & lower_boundary, double & upper_boundary, double & kappa, double & gamma, const double & c, const double & G, const double & pi);
void Multiple(std::string & p_to_vary, double & radius_min, double & radius_max, double & mass_min, double & mass_max, unsigned long & nns, std::string & file_name, std::string & method, double & initial_guess, double & maximal_radius, double & total_mass, double & frac, double & dr, double & stop, double & dr_div, unsigned long long max_steps, double & mass_resolution, double & lower_boundary, double & upper_boundary, double & kappa, double & gamma, const double & c, const double & G, const double & pi);

int main(){
  const double c = 2.99792458e10; // speed of light
  const double G = 6.67259e-8; // gravitational constant
  const double pi = 3.1416926;

  // read in parameter from config file
  std::ifstream cf("cf");

  std::string parameter, value;
  std::vector < std::string > v_p, v_v;

  if(cf.is_open()){
    while(cf >> parameter >> value){
      v_p.push_back(parameter);
      v_v.push_back(value);
    }
  }
  cf.close();

  // Default parameters
  double maximal_radius = 20.0; // in km
  double total_mass = 2.2; // in solar masses
  double mass_resolution = 1e-3;
  double kappa = 1.98183e-6;
  double gamma = 2.75;
  double radius_min = 12.0;
  double radius_max = 22.0;
  double mass_min = 1.1;
  double mass_max = 2.1;
  unsigned long nns = 10;
  std::string method = "Rk4";
  std::string p_to_vary = "M";
  std::string option = "Single";

  // Hidden parameter
  std::string file_name = "Single.csv";
  std::string file_nameX = "Multiple.csv";
  double frac = 1e-4;
  double dr = 400.0;
  double stop = 1.0;
  double dr_div = 350.0;
  unsigned long long max_steps = 1e5;
  double lower_boundary = 5.0e13;
  double upper_boundary = 5.0e16;
  double initial_guess = 3.0e14;

  // update parameters set by user
  if(v_p.size() == v_v.size()){
    for(std::size_t i=0; i<v_p.size(); i++){
      if(v_p.at(i).compare("maximal_radius") == 0) maximal_radius = std::stod(v_v.at(i));
      if(v_p.at(i).compare("total_mass") == 0) total_mass = std::stod(v_v.at(i));
      if(v_p.at(i).compare("mass_resolution") == 0) mass_resolution = std::stod(v_v.at(i));
      if(v_p.at(i).compare("kappa") == 0) kappa = std::stod(v_v.at(i));
      if(v_p.at(i).compare("gamma") == 0) gamma = std::stod(v_v.at(i));
      if(v_p.at(i).compare("radius_max") == 0) radius_max = std::stod(v_v.at(i));
      if(v_p.at(i).compare("radius_min") == 0) radius_min = std::stod(v_v.at(i));
      if(v_p.at(i).compare("mass_min") == 0) mass_min = std::stod(v_v.at(i));
      if(v_p.at(i).compare("mass_max") == 0) mass_max = std::stod(v_v.at(i));
      if(v_p.at(i).compare("nns") == 0) nns = (unsigned long) std::stod(v_v.at(i));
      if(v_p.at(i).compare("mode") == 0) method = v_v.at(i);
      if(v_p.at(i).compare("to_vary") == 0) p_to_vary = v_v.at(i);
      if(v_p.at(i).compare("option") == 0) option = v_v.at(i);

      if(v_p.at(i).compare("file_nameS") == 0) file_name = v_v.at(i);
      if(v_p.at(i).compare("file_nameM") == 0) file_nameX = v_v.at(i);
      if(v_p.at(i).compare("frac") == 0) frac = std::stod(v_v.at(i));
      if(v_p.at(i).compare("dr") == 0) dr = std::stod(v_v.at(i));
      if(v_p.at(i).compare("stop") == 0) stop = std::stod(v_v.at(i));
      if(v_p.at(i).compare("dr_div") == 0) dr_div = std::stod(v_v.at(i));
      if(v_p.at(i).compare("max_steps") == 0) max_steps = (unsigned long long) std::stod(v_v.at(i));
      if(v_p.at(i).compare("lower_boundary") == 0) lower_boundary = std::stod(v_v.at(i));
      if(v_p.at(i).compare("upper_boundary") == 0) upper_boundary = std::stod(v_v.at(i));
      if(v_p.at(i).compare("initial_guess") == 0) initial_guess = std::stod(v_v.at(i));
    }
  }

  if(option.compare("Single") == 0) Adapt(file_name, method, initial_guess, maximal_radius, total_mass, frac, dr, stop, dr_div, max_steps, mass_resolution, lower_boundary, upper_boundary, kappa, gamma, c, G, pi);
  if(option.compare("Multiple") == 0) Multiple(p_to_vary, radius_min, radius_max, mass_min, mass_max, nns, file_nameX, method, initial_guess, maximal_radius, total_mass, frac, dr, stop, dr_div, max_steps, mass_resolution, lower_boundary, upper_boundary, kappa, gamma, c, G, pi);



  return 0;
}


// Function definitions
inline double GetDensity(double & pressure, double & kappa, double & gamma){
  return pow(pressure / kappa, 1.0 / gamma);
}
inline double GetPressure(double & density, double & kappa, double & gamma){
  return kappa * pow(density, gamma);
}
inline double GetPressure(const double & density, double & kappa, double & gamma){ // needed for initial estimate only
  return kappa * pow(density, gamma);
}
double TOV(double radius, double & mass, double pressure, const double & c, const double & G, const double & pi, double & kappa, double & gamma){
  double t1 = -G / pow(radius, 2.0);
  double t2 = GetDensity(pressure, kappa, gamma) + pressure / pow(c, 2.0);
  double t3 = mass + 4.0 * pi * pow(radius, 3.0) * pressure / pow(c, 2.0);
  double t4 = 1.0 - 2.0 * G * mass / pow(c, 2.0) / radius;
  return t1 * t2 * t3 / t4;
}

double Euler(double & radius, double & mass, double & pressure, double & dr, const double & c, const double & G, const double & pi, double & kappa, double & gamma){
  // The explicit Euler method
  //
  // 1 | 0
  // =====
  //   | 1
  return dr * TOV(radius+0.0, mass, pressure+0.0, c, G, pi, kappa, gamma);
}
double Rk2(double & radius, double & mass, double & pressure, double & dr, const double & c, const double & G, const double & pi, double & kappa, double & gamma){
  // The Runge-Kutta method of second order
  //
  // 0   |
  // 0.5 | 0.5
  // =============
  //     | 0   | 1

  double k1 = dr * TOV(radius+0.0,          mass, pressure+0.0,          c, G, pi, kappa, gamma);
  double k2 = dr * TOV(radius + dr/2.0, mass, pressure + k1/2.0, c, G, pi, kappa, gamma);
  return k2;
}
double Rk3(double & radius, double & mass, double & pressure, double & dr, const double & c, const double & G, const double & pi, double & kappa, double & gamma){
  // The Runge-Kutta method of third order
  //
  // 0   |
  // 0.5 | 0.5 |
  // 1   | -1  |  2  |
  // =====================
  //     | 1/6 | 4/6 | 1/6

  double k1 = dr * TOV(radius+0.0,          mass, pressure+0.0,               c, G, pi, kappa, gamma);
  double k2 = dr * TOV(radius + dr/2.0, mass, pressure + k1/2.0,      c, G, pi, kappa, gamma);
  double k3 = dr * TOV(radius + dr,     mass, pressure - k1 + 2.0*k2, c, G, pi, kappa, gamma);
  return (k1 + 4.0*k2 + k3)/6.0;
}
double Rk4(double & radius, double & mass, double & pressure, double & dr, const double & c, const double & G, const double & pi, double & kappa, double & gamma){
  // The classical Runge-Kutta method of fourth order
  //
  // 0   |
  // 0.5 | 0.5 |
  // 0.5 | 0   | 0.5 |
  // 1   | 0   | 0   | 1   |
  // ===========================
  //     | 1/6 | 2/6 | 2/6 | 1/6

  double k1 = dr * TOV(radius + 0.0,    mass, pressure + 0.0,    c, G, pi, kappa, gamma);
  double k2 = dr * TOV(radius + dr/2.0, mass, pressure + k1/2.0, c, G, pi, kappa, gamma);
  double k3 = dr * TOV(radius + dr/2.0, mass, pressure + k2/2.0, c, G, pi, kappa, gamma);
  double k4 = dr * TOV(radius + dr/2.0, mass, pressure + k3,     c, G, pi, kappa, gamma);
  return (k1 + 2.0*k2 + 2.0*k3 + k4)/6.0;
}


void Adapt(std::string & file_name, std::string & method, double & initial_guess, double & maximal_radius, double & total_mass, double & frac, double & dr, double & stop, double & dr_div, unsigned long long max_steps, double & mass_resolution, double & lower_boundary, double & upper_boundary, double & kappa, double & gamma, const double & c, const double & G, const double & pi){
  double mass_ratio = 100.0; // something other than pm below
  double bl = lower_boundary; // will be adapted
  double bc = initial_guess; // will be adapted
  double bu = upper_boundary; // will be adapted

  double pm = 0.0; // previous mass

  std::vector < std::vector < double > > results; // from Single-function

  // another break criteria
  double break_res = 0.1;

 // if((residual1 < break_res) || (residual2 < break_res)) result.back().at(4) = -2000.0;

  while(true){
    results = Single(file_name, method, bc, maximal_radius, total_mass, frac, dr, stop, dr_div, max_steps, mass_resolution, bl, bu, kappa, gamma, c, G, pi);

    mass_ratio = results.back().at(2) / (total_mass * 1.9989e33);

    results.clear();

    // break loop if mr is finde
    if(std::abs(mass_ratio - 1.0) <= mass_resolution) break; // here we will save data. This should be done in the calculation above already to save computation time..

    // check if something went wrong..
    if(mass_ratio <= 0.01 || (mass_ratio == pm) || std::isnan(mass_ratio)) break;

    // if residual1 or residual2 is below break res, stop iteration.
    double residual1 = std::abs(GetDensity(results.front().at(3), kappa, gamma) / upper_boundary - 1.0);
    double residual2 = std::abs(GetDensity(results.front().at(3), kappa, gamma) / lower_boundary - 1.0);

    if(residual1 < break_res || residual2 < break_res) break;

    // set new previous mass
    pm = mass_ratio;

    // Set new boundaries
    if(mass_ratio - 1.0 < 0.0){
      bl = bc;
      bc = sqrt(bl*bu);
    }
    if(mass_ratio - 1.0 > 0.0){
      bu = bc;
      bc = sqrt(bl*bu);
    }
  }

  return ;
}


std::vector < std::vector < double > > Single(std::string & file_name, std::string & method, double & initial_guess, double & maximal_radius, double & total_mass, double & frac, double & dr, double & stop, double & dr_div, unsigned long long max_steps, double & mass_resolution, double & lower_boundary, double & upper_boundary, double & kappa, double & gamma, const double & c, const double & G, const double & pi){
  // Parameters for single calculation
  double central_pressure = GetPressure(initial_guess, kappa, gamma);
  double mass_ = 0.0;
  double radius_ = 0.0;
  double pre1 = central_pressure;
  double pre2 = central_pressure;
  double pre3 = 0.0;
  double nsteps = 0.0;
  double dr_ = dr;

  // Change the number of digits to print out
  std::setprecision(10);

  std::vector < std::vector < double > > result;

  while(true){
    std::vector < double > temp;
    temp.push_back(nsteps);
    temp.push_back(radius_);
    temp.push_back(mass_);
    temp.push_back(pre1);
    temp.push_back(GetDensity(pre1, kappa, gamma));
    temp.push_back(initial_guess);
    temp.push_back(dr_);

    result.push_back(temp);

    if(pre1 <= frac*pre2) {
      dr_ /= dr_div;
      pre2 = pre1;
    }

    radius_ += dr_;
    mass_ += 4.0 * pi * pow(radius_, 2.0) * dr_ * GetDensity(pre1, kappa, gamma);
    nsteps += 1.0;

    if(method.compare("Euler") == 0) pre3 = Euler(radius_, mass_, pre1, dr_, c, G, pi, kappa, gamma);
    if(method.compare("Rk2") == 0) pre3 = Rk2(radius_, mass_, pre1, dr_, c, G, pi, kappa, gamma);
    if(method.compare("Rk3") == 0) pre3 = Rk3(radius_, mass_, pre1, dr_, c, G, pi, kappa, gamma);
    if(method.compare("Rk4") == 0) pre3 = Rk4(radius_, mass_, pre1, dr_, c, G, pi, kappa, gamma);

    if((pre1 < stop) || (radius_ / 1e5 >= maximal_radius) || (static_cast < unsigned long long > (nsteps) > max_steps)) break;
    else pre1 += pre3;
  }

  if(std::abs(result.back().at(2) / (total_mass * 1.9989e33)) - 1.0 <= mass_resolution){
    std::ofstream ofile(file_name);
    ofile << "nsteps,radius,mass,pressure,density,initial density,dr," << std::endl;
    for(std::size_t i=0; i<result.size(); i++){
      for(std::size_t j=0; j<result.at(i).size(); j++){
        ofile << result.at(i).at(j) << ",";
      }
      ofile << std::endl;
    }
  }

  return result;
}


void Multiple(std::string & p_to_vary, double & radius_min, double & radius_max, double & mass_min, double & mass_max, unsigned long & nns, std::string & file_name, std::string & method, double & initial_guess, double & maximal_radius, double & total_mass, double & frac, double & dr, double & stop, double & dr_div, unsigned long long max_steps, double & mass_resolution, double & lower_boundary, double & upper_boundary, double & kappa, double & gamma, const double & c, const double & G, const double & pi){
  // Output file
  std::ofstream ofile(file_name);
  ofile << "nsteps,radius,mass,pressure,density,central density,dr,," << std::endl;

  // parameter
  std::vector < std::vector < std::string > > results, temp2;
  std::string file_name2 = "temorary.csv";
  double nsteps2 = 0;
  double radius2 = 0;
  double mass2 = 0;
  double pressure2 = 0;
  double density2 = 0;
  double dr2 = 0;

  if(p_to_vary.compare("M") == 0){
    double increment = (mass_max - mass_min) / static_cast < double > (nns);

    for(std::size_t i=0; i<nns+1; i++){
      std::vector < std::string > temp;
      double new_mass = mass_max - i * increment;

      Adapt(file_name2, method, initial_guess, maximal_radius, new_mass, frac, dr, stop, dr_div, max_steps, mass_resolution, lower_boundary, upper_boundary, kappa, gamma, c, G, pi);

      // open temporary file, move to second to last line and read in last line
      std::ifstream ifile(file_name2);
      while(ifile >> row){
        for(std::size_t i=0; i<row.size(); i++) temp.push_back(row[i]);
        temp2.push_back(temp);
        temp.clear();
      }
      for(std::size_t j=0; j<temp2.back().size(); j++) temp.push_back(temp2.back().at(j));
      results.push_back(temp);

      temp.clear();
      ifile.close();
    }
  }

  if(p_to_vary.compare("R") == 0){
    double increment = (radius_max - radius_min) / static_cast < double > (nns);

    for(std::size_t i=0; i<nns+1; i++){
      std::vector < std::string > temp;
      double new_radius = radius_max - i * increment;

      Adapt(file_name2, method, initial_guess, new_radius, total_mass, frac, dr, stop, dr_div, max_steps, mass_resolution, lower_boundary, upper_boundary, kappa, gamma, c, G, pi);

      // open temporary file, move to second to last line and read in last line
      std::ifstream ifile(file_name2);
      while(ifile >> row){
        for(std::size_t i=0; i<row.size(); i++) temp.push_back(row[i]);
        temp2.push_back(temp);
        temp.clear();
      }
      for(std::size_t j=0; j<temp2.back().size(); j++) temp.push_back(temp2.back().at(j));
      results.push_back(temp);

      temp.clear();
      ifile.close();
    }
  }


  // store results in output file
  for(std::size_t i=0; i<results.size(); i++){
    for(std::size_t j=0; j<results.at(i).size(); j++){
      ofile << results.at(i).at(j) << ",";
    }
    ofile << std::endl;
  }

  ofile.close();

  return ;
}
