#include "grid.hpp"
#include "highfive/H5DataSet.hpp"
#include "highfive/H5Exception.hpp"
#include <array>
#include <boost/multi_array.hpp>
#include <csignal>
#include <exception>
#include <fstream>
#include <highfive/H5Easy.hpp>
#include <highfive/H5File.hpp>
#include <iomanip>
#include <iostream>
#include <map>
// #include <omp.h>
#include <stdexcept>
#include <string>
#include <unistd.h>
#include <utility>
#include <vector>

void signalHandler(int signum) {
  std::cout << "Interrupt signal (" << signum << ") received.\n";
  globals::Stop = true;
}

int main() {
  double J = 1.0;
  const int N = 1e6;
  const std::string FILENAME = "heisenberg_model_wolff_updated.h5";
  std::cout << std::setprecision(10);
  const size_t multihitparam = 7;
  const size_t steps = 8;
  const size_t CORES = 8;
  const size_t INSTANCE_NUM{0};
  const size_t INSTANCES{3};

  // std::vector<Grid3D*> grids;

  /* TODO: Need to make a smaller max value and finer steps around the critical
   * temperatur*/
  double increment = 0.0038 / 3;
  double start_value = 0.67;

  if (INSTANCE_NUM >= INSTANCES) {
    throw std::runtime_error("Instance number is too high");
  }

  std::signal(SIGTERM, signalHandler);

  std::cout << "Lowest beta: " << start_value
            << "Highest beta: " << start_value + increment * steps << std::endl;

  if (steps % CORES != 0) {
    throw std::runtime_error(
        "For Batches larger than N, use a number, that is divisible by " +
        std::to_string(CORES));
  }

  {
    H5Easy::File file(FILENAME, H5Easy::File::ReadWrite | H5Easy::File::Create);
  }

  std::vector<size_t> L_list = {16, 21, 24, 27, 29};

  std::vector<std::pair<size_t, double>> L_list_beta;
  for (size_t L : L_list) {
    for (int i = 0; i < steps; i++) {
      double beta =
          start_value + increment * i * INSTANCES + increment * INSTANCE_NUM;
      L_list_beta.push_back(std::make_pair(L, beta));
    }
  }

#pragma omp parallel for schedule(static, 1)
  for (auto &[L, beta] : L_list_beta) {
    time_t now = time(0);
    char *dt = ctime(&now);
    std::cout << "Start: Local date and time is: " << dt
              << "\nStarting L = " << L << "Beta= " << beta << std::endl;

    Grid3D grid(L, L, L, J, beta, multihitparam, N, FILENAME);
    grid.mainloop(true);

    now = time(0);
    dt = ctime(&now);
    std::cout << "Finished L = " << L << "Beta= " << beta
              << "\nEnd: Local date and time is: " << dt << std::endl;
  }

  return 0;
}