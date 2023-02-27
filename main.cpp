#include "grid.hpp"
#include "highfive/H5DataSet.hpp"
#include "highfive/H5Exception.hpp"
#include <array>
#include <boost/multi_array.hpp>
#include <fstream>
#include <highfive/H5Easy.hpp>
#include <highfive/H5File.hpp>
#include <iomanip>
#include <iostream>
#include <omp.h>
#include <stdexcept>
#include <string>
#include <unistd.h>
#include <vector>

void mainloop(double J, double beta, const int N, size_t L,
              std::string const FILENAME, const size_t multihitparam) {
  std::cout << "beta: " << beta << std::endl;

  
  Grid3D grid(L, L, L, J, beta, multihitparam, N, FILENAME);
  // std::cout << grid << std::endl;

  std::string Groupname = "size_" + std::to_string(L) + "/";
  std::string Subgroup = "beta_" + std::to_string(beta) + "/";
  std::string Snapshots = "snapshots/";
  std::vector<double> mag;
  std::vector<double> energy;
  mag.reserve(N);
  energy.reserve(N);
  for (int i = 0; i < N; i++) {
    grid.hb_sweep();
    /*TODO: Move Vector of magnetization and energy to grid class -> Write to
     * file in destructor*/
    mag.push_back(grid.calculate_magnetization());
    energy.push_back(grid.calculate_energy());

    if (i % 2000 == 0) {
      boost::multi_array<double, 4> snapshot = grid.snapshot();
#pragma omp critical
      {
        HighFive::File file(FILENAME, HighFive::File::ReadWrite);
        HighFive::DataSet data_set = file.createDataSet<double>(
            Groupname + Subgroup + Snapshots + std::to_string(i),
            HighFive::DataSpace::From(snapshot));
        data_set.write(snapshot);
        //   H5Easy::dump(file, Groupname + Subgroup + Snapshots +
        //   std::to_string(i), snapshot,
        //                H5Easy::DumpMode::Overwrite);
      }
    }
  }

#pragma omp critical
  {
    try {

      H5Easy::File file(FILENAME, H5Easy::File::ReadWrite);

      H5Easy::dump(file, Groupname + Subgroup + "Magnetization", mag,
                   H5Easy::DumpMode::Overwrite);
      H5Easy::dump(file, Groupname + Subgroup + "Energy", energy,
                   H5Easy::DumpMode::Overwrite);

    } catch (HighFive::Exception const &err) {
      std::cout << err.what() << std::endl;
      H5Easy::File file(FILENAME, H5Easy::File::ReadWrite);

      H5Easy::dump(file, Groupname + Subgroup + "Magnetization_new", mag);
      H5Easy::dump(file, Groupname + Subgroup + "Energy_new", energy);
    }
  }
}

int main() {
  double J = 1.0;
  double beta = 0.5;
  const int N = 1e4;
  size_t L = 11;
  const std::string FILENAME = "heisenberg_model_test_new.h5";
  std::cout << std::setprecision(10);
  const size_t multihitparam = 7;
  size_t steps = 12;

  {
    H5Easy::File file(FILENAME, H5Easy::File::ReadWrite | H5Easy::File::Create);
  }
  //   std::cout << grid.snapshot() << std::endl;


#pragma omp parallel for
  for (int i = 0; i < steps; i++) {

    if(steps > 12 && steps % 12 != 0){
    throw std::runtime_error("For Batches larger than N, use a number, that is divisible by 12.");
  }


    double beta = 0.60 + 0.03*i;

    Grid3D grid(L, L, L, J, beta, multihitparam, N, FILENAME);
    grid.mainloop();
    // sleep for 1 sec

   
  }

  // // std::cout << grid << std::endl;
  // std::cout << grid.calculate_magnetization() << std::endl;
  // // open fstream
  // std::ofstream mag_file;
  // std::ofstream energy_file;
  // mag_file.open("mag.dat");
  // energy_file.open("energy.dat");
  // for (int i = 0; i < N; i++) {
  //     grid.hb_sweep();
  //     mag_file << grid.calculate_magnetization() << "\n";
  //     energy_file << grid.calculate_energy() << "\n";

  // }

  // mag_file.close();
  // energy_file.close();

  //   mainloop(J, beta, N, L, FILENAME);

  return 0;
}