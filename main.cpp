#include "grid.hpp"
#include "highfive/H5DataSet.hpp"
#include "highfive/H5Exception.hpp"
#include <array>
#include <boost/multi_array.hpp>
#include <exception>
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
#include <csignal>



void signalHandler(int signum) {
  std::cout << "Interrupt signal (" << signum << ") received.\n";
  globals::Stop = true;

  sleep(20e3);
  
 
  
  exit(signum);
  
}


int main() {
  double J = 1.0;
  const int N = 1e6;
  const std::string FILENAME = "heisenberg_model_wolff_test.h5";
  std::cout << std::setprecision(10);
  const size_t multihitparam = 7;
  const size_t steps = 24;
  const size_t CORES = 12;

  // std::vector<Grid3D*> grids;

  double start_value = 0.55;
  double increment = 0.02;

  std::signal(SIGTERM, signalHandler);

  std::cout << "Lowest beta: " << start_value << "Highest beta: " << start_value + increment * steps << std::endl;
  

  if(steps % CORES != 0){
    throw std::runtime_error("For Batches larger than N, use a number, that is divisible by " + std::to_string(CORES));
  }

  {
    H5Easy::File file(FILENAME, H5Easy::File::ReadWrite | H5Easy::File::Create);
  }
  //   std::cout << grid.snapshot() << std::endl;

  // create a list of L values containing 25, 30, 35, 40, 45, 50, 55, 60
  std::vector<size_t> L_list = {20};

  

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


        grid.mainloop(true);
        

      
      }
      std::cout << "Finished L = " << L << std::endl;
      now = time(0);
      dt = ctime(&now);
      std::cout << "End: Local date and time is: " << dt << std::endl;
  }
  return 0;
}