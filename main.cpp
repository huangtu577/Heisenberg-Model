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
  const int N = 5e6;
  const std::string FILENAME = "heisenberg_model.h5";
  std::cout << std::setprecision(10);
  const size_t multihitparam = 7;
  size_t steps = 16;
  const size_t CORES = 8;

  double start_value = 0.55;
  double increment = 0.02;

  std::cout << "Lowest beta: " << start_value << "Highest beta: " << start_value + increment * steps << std::endl;
  

  if(steps % CORES != 0){
    throw std::runtime_error("For Batches larger than N, use a number, that is divisible by " + std::to_string(CORES));
  }

  {
    H5Easy::File file(FILENAME, H5Easy::File::ReadWrite | H5Easy::File::Create);
  }
  //   std::cout << grid.snapshot() << std::endl;

  // create a list of L values containing 25, 30, 35, 40, 45, 50, 55, 60
  std::vector<size_t> L_list = {20, 30, 35, 40};

  

  for(size_t L: L_list){
    //Print out the time
    time_t now = time(0);
    char* dt = ctime(&now);
    std::cout << "Start: Local date and time is: " << dt << std::endl;
    
    std::cout << "Starting L = " << L << std::endl;
    
    
    
    #pragma omp parallel for
      for (int i = 0; i < steps; i++) {
        

        double beta = start_value + increment *i;

        Grid3D grid(L, L, L, J, beta, multihitparam, N, FILENAME);
        grid.mainloop();

      
      }
      std::cout << "Finished L = " << L << std::endl;
      now = time(0);
      dt = ctime(&now);
      std::cout << "End: Local date and time is: " << dt << std::endl;
  }
  return 0;
}