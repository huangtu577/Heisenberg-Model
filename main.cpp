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
// #include <omp.h>
#include <stdexcept>
#include <string>
#include <unistd.h>
#include <vector>



int main() {
  double J = 1.0;
  const int N = 1e5;
  size_t L = 11;
  const std::string FILENAME = "heisenberg_model_test.h5";
  std::cout << std::setprecision(10);
  const size_t multihitparam = 7;
  size_t steps = 12;
  const size_t CORES = 12;

  double start_value = 0.55;
  double increment = 0.02;

  {
    H5Easy::File file(FILENAME, H5Easy::File::ReadWrite | H5Easy::File::Create);
  }
  // Grid3D grid(L, L, L, J, beta, multihitparam, N, FILENAME);
  // grid.mainloop(true);
  // return 0;

  //   std::cout << grid.snapshot() << std::endl;

  // create a list of L values containing 25, 30, 35, 40, 45, 50, 55, 60
  // std::vector<size_t> L_list = {15, 25, 30, 35, 40};
  std::vector<size_t> L_list = {15};

  for (size_t L : L_list) {
    std::cout << "Starting L = " << L << std::endl;

#pragma omp parallel for
    for (int i = 0; i < steps; i++) {

      double beta = start_value + increment * i;

      Grid3D grid(L, L, L, J, beta, multihitparam, N, FILENAME);
      grid.mainloop(true);
    }
    std::cout << "Finished L = " << L << std::endl;
  }
  return 0;
}