#pragma ONCE
#include "highfive/H5Easy.hpp"
#include "highfive/H5File.hpp"
#include "spin.hpp"
#include <algorithm>
#include <array>
#include <boost/multi_array.hpp>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <numbers>
#include <random>
#include <string>
#include <vector>
#include <stack>

class Grid {
protected:
  size_t dim;
  std::random_device rd;
  std::mt19937 gen{rd()};
  double J;
  double beta;
  size_t multihitparam{1};
  size_t N;
  const std::string filename;

  double acceptance_rate{0};

  std::string groupname, subgroup, snapshots;

  Grid(size_t dim, double J, double beta, size_t multihitparam, size_t N,
       std::string filename)
      : dim{dim}, J{J}, beta(beta), multihitparam(multihitparam), N{N},
        filename(filename) {
    assert(multihitparam > 0 && multihitparam % 2 == 1);

    subgroup = "beta_" + std::to_string(beta) + "/";
    snapshots = "snapshots/";
  };

public:
  std::uniform_real_distribution<double> dis_phi{0.0, 2 * std::numbers::pi};
  std::uniform_real_distribution<> dis_theta{0.0, std::numbers::pi};
  std::uniform_real_distribution<double> dis_probability{0.0, 1.0};

  // virtual ~Grid() = 0;
};

class Grid3D : Grid {
public:
  typedef boost::multi_array<Spin3D, 3> array_type;

  array_type grid;

  std::uniform_int_distribution<> dis_x;
  std::uniform_int_distribution<> dis_y;
  std::uniform_int_distribution<> dis_z;

  std::vector<double> magnetization, energy;

  double Z;

  Grid3D(const size_t x_val, const size_t y_val, const size_t z_val,
         const double J, const double beta, const size_t multihitparam,
         const size_t N, const std::string filename, bool random = true)
      : x(x_val), y(y_val), z(z_val), grid(boost::extents[z_val][y_val][x_val]),
        dis_x(0, x_val - 1), dis_y(0, y_val - 1),
        dis_z(0, z_val - 1), Grid{3, J, beta, multihitparam, N, filename} {

    groupname = "size_" + std::to_string(x) + "/";

    // const size_t x{x_val}, y{y_val}, z{z_val};
    // boost::array<array_type::index, 3> shape = {
    //     {static_cast<long>(x), static_cast<long>(y), static_cast<long>(z)}};
    // array_type grid(shape);
    // std::cout << grid.data() << std::endl;

    // Fill the array with random Values
    for (array_type::index i = 0; i < z; ++i) {
      for (array_type::index j = 0; j < y; ++j) {
        for (array_type::index k = 0; k < x; ++k) {
          if (random){
            this->grid[i][j][k].update_coords(dis_theta(gen), dis_phi(gen));

          }
          else{
            this->grid[i][j][k].update_coords(0, 0);
          }
          // std::cout << this->grid[i][j][k] << std::endl;
        }
      }
    }

    this->magnetization.reserve(this->N);
    this->energy.reserve(this->N);

    // std::cout << *grid.data() << std::endl;

    // grid = grid_tmp;
    this->Z = this->calculate_partition_function();
    std::cout << "Grid constructor: L: " << this->x << ", Beta: "<< beta << std::endl;
  }

  ~Grid3D() {

#pragma omp critical
    {
      try {

        H5Easy::File file(this->filename, H5Easy::File::ReadWrite);

        H5Easy::dump(file, this->groupname + this->subgroup + "Magnetization",
                     this->magnetization, H5Easy::DumpMode::Overwrite);
        H5Easy::dump(file, this->groupname + this->subgroup + "Energy",
                     this->energy, H5Easy::DumpMode::Overwrite);

      } catch (HighFive::Exception const &err) {
        std::cout << err.what() << std::endl;
        H5Easy::File file(this->filename, H5Easy::File::ReadWrite);

        H5Easy::dump(file,
                     this->groupname + this->subgroup + "Magnetization_new",
                     this->magnetization);
        H5Easy::dump(file, this->groupname + this->subgroup + "Energy_new",
                     this->energy);
      }
    }

    std::cout << "Grid destructor" << std::endl;
    std::cout << "Acceptance rate: " << this->acceptance_rate/(N*this->x*this->y*this->z) << std::endl;
    std::cout << "Accepted Changes: " << this->acceptance_rate << std::endl;
  }

  /* Solution to this Integral comes from Wolfram Alpha*/
  double calculate_partition_function() {
    return   std::pow(2*std::numbers::pi * std::sinh(this->beta * this->J) /
                               (this->beta * this->J),
                           1);
  }

  void hb_sweep() {

    for (int i = 0; i < this->z; i++) {
      for (int j = 0; j < this->y; j++) {
        for (int k = 0; k < this->x; k++) {

          Spin3D s(dis_theta(gen), dis_phi(gen));

          double E{0.0};
          E -= this->J * (s * this->grid[(i + 1 + this->z) % this->z][j][k] +
                          s * this->grid[(i - 1 + this->z) % this->z][j][k] +
                          s * this->grid[i][(j + 1 + this->y) % this->y][k] +
                          s * this->grid[i][(j - 1 + this->y) % this->y][k] +
                          s * this->grid[i][j][(k + 1 + this->x) % this->x] +
                          s * this->grid[i][j][(k - 1 + this->x) % this->x]);
          double p{std::exp(-beta * J * E) / this->Z};
          // std::cout << E << std::endl;
          assert(0 <= p && p <= 1);
          double r = dis_probability(gen);

          for (int n = 0; n < this->multihitparam; ++n) {
            if (r < p) {
              this->acceptance_rate++;
              this->grid[i][j][k] = s;
              break;
            }
          }
        }
      }
    }

    this->magnetization.push_back(this->calculate_magnetization());
    this->energy.push_back(this->calculate_energy());


  }

  void wolf_sweep() {
    const Spin3D r(dis_theta(gen), dis_phi(gen));
    std::vector<std::array<size_t, 3>> visited{};
    visited.reserve(this->x*this->y*this->z);

    std::stack<std::array<size_t, 3>> stack;

    /* Random start sites*/
    double x(dis_x(gen)), y(dis_y(gen)), z(dis_z(gen));

    Spin3D s = this->grid[static_cast<int>(z)][static_cast<int>(y)]
                         [static_cast<int>(x)];
    spin_flip_wolf(r, this->grid[static_cast<int>(z)][static_cast<int>(y)]
                         [static_cast<int>(x)]);
    stack.push({static_cast<size_t>(z), static_cast<size_t>(y),
                static_cast<size_t>(x)});

    while (!stack.empty()) {
      std::array<size_t, 3> current{stack.top()};
      if(std::find(visited.begin(), visited.end(), current) != visited.end()){
        stack.pop();
        continue;
      }
      visited.push_back(current);
      stack.pop();
      Spin3D current_s = this->grid[current[0]][current[1]][current[2]];

      /* Calculate NN Interaction*/
      Spin3D neighbor;

      std::array<std::array<size_t, 3>, 6> nn = calc_nn(current);

      for (int i = 0; i < 6; i++) {
        neighbor = this->grid[nn[i][0]][nn[i][1]][nn[i][2]];
        double q{dis_probability(gen)};
        double p{1 - std::exp(2 * this->beta * this->J * (current_s * r)*(neighbor*r))};

        // std::cout << p << std::endl;

        
          
          if (q < p) {
            this->acceptance_rate++;
            spin_flip_wolf(r, this->grid[nn[i][0]][nn[i][1]][nn[i][2]]);
            stack.push(nn[i]);
          }
        }
      }
      this->magnetization.push_back(this->calculate_magnetization());
      this->energy.push_back(this->calculate_energy());
    }
  

    std::array<std::array<size_t, 3>, 6> calc_nn(
        std::array<size_t, 3> current) {
      std::array<std::array<size_t, 3>, 6> nn;
      nn[0] = {(current[0] + 1 + this->z) % this->z, current[1], current[2]};
      nn[1] = {(current[0] - 1 + this->z) % this->z, current[1], current[2]};
      nn[2] = {current[0], (current[1] + 1 + this->y) % this->y, current[2]};
      nn[3] = {current[0], (current[1] - 1 + this->y) % this->y, current[2]};
      nn[4] = {current[0], current[1], (current[2] + 1 + this->x) % this->x};
      nn[5] = {current[0], current[1], (current[2] - 1 + this->x) % this->x};
      return nn;
    }

    /* Performs the Spin Flip of s*/
    void spin_flip_wolf(const Spin3D &r, Spin3D &s) {
      double a = 2 * (s * r);
      // s.update_coords(s.spin[0] - a*r.spin[0], s.spin[1] - a*r.spin[1]);
      // return;

      s.spin[0] -= a * r.spin[0];
      s.spin[1] -= a * r.spin[1];
    }

    void mainloop(bool wolff=false) {

      if(!wolff){

      for (int i = 0; i < this->N; i++) {
        this->hb_sweep();

        if (i % 100000 == 0) {
          boost::multi_array<double, 4> snapshot = this->snapshot();
        #pragma omp critical
        {
          HighFive::File file(this->filename, HighFive::File::ReadWrite);
          HighFive::DataSet data_set = file.createDataSet<double>(
              this->groupname + this->subgroup + this->snapshots +
                  std::to_string(i),
              HighFive::DataSpace::From(snapshot));
          data_set.write(snapshot);
        }
      }
    }
      }
      else{
        for (int i = 0; i < this->N; i++) {
        this->wolf_sweep();

        if (i % 100000 == 0) {
          boost::multi_array<double, 4> snapshot = this->snapshot();
        #pragma omp critical
        {
          HighFive::File file(this->filename, HighFive::File::ReadWrite);
          HighFive::DataSet data_set = file.createDataSet<double>(
              this->groupname + this->subgroup + this->snapshots +
                  std::to_string(i),
              HighFive::DataSpace::From(snapshot));
          data_set.write(snapshot);
        }
      }
    }
      }

  }

  double calculate_magnetization() {
    double m{0.0};
    double x{0.0}, y{0.0}, z{0.0};
    Spin3D s;
    for (array_type::index i = 0; i < this->z; ++i) {
      for (array_type::index j = 0; j < this->y; ++j) {
        for (array_type::index k = 0; k < this->x; ++k) {
          s = this->grid[i][j][k];
          x += std::sin(s.spin[0]) * std::cos(s.spin[1]);
          y += std::sin(s.spin[0]) * std::sin(s.spin[1]);
          z += std::cos(s.spin[0]);
        }
      }
    }
    /* Calculate the direction of Magnetization*/
    m = std::abs(1 / static_cast<double>(this->x * this->y * this->z) *
                 std::sqrt(std::pow(x, 2) + std::pow(y, 2) + std::pow(z, 2)));
    return m;
  }

  double calculate_energy() {
    double h{0};
    for (array_type::index i = 0; i < this->z; ++i) {
      for (array_type::index j = 0; j < this->y; ++j) {
        for (array_type::index k = 0; k < this->x; ++k) {
          h -= this->J/2 * (this->grid[i][j][k] *
                              this->grid[(i + 1 + this->z) % this->z][j][k] +
                          this->grid[i][j][k] *
                              this->grid[(i - 1 + this->z) % this->z][j][k] +
                          this->grid[i][j][k] *
                              this->grid[i][(j + 1 + this->y) % this->y][k] +
                          this->grid[i][j][k] *
                              this->grid[i][(j - 1 + this->y) % this->y][k] +
                          this->grid[i][j][k] *
                              this->grid[i][j][(k + 1 + this->x) % this->x] +
                          this->grid[i][j][k] *
                              this->grid[i][j][(k - 1 + this->x) % this->x]);
        }
      }
    }
    return 1 / static_cast<double>(this->x * this->y * this->z) * h;
  }

  boost::multi_array<double, 4> snapshot() {
    boost::multi_array<double, 4> snapshot(
        boost::extents[this->z][this->y][this->x][2]);
    for (array_type::index i = 0; i < this->z; ++i) {
      for (array_type::index j = 0; j < this->y; ++j) {
        for (array_type::index k = 0; k < this->x; ++k) {
          snapshot[i][j][k][0] = this->grid[i][j][k].spin[0];
          snapshot[i][j][k][1] = this->grid[i][j][k].spin[1];
        }
      }
    }
    return snapshot;
  }

protected:
  const size_t x, y, z;

  friend std::ostream &operator<<(std::ostream &out, const Grid3D &grid);
  };

std::ostream &operator<<(std::ostream &out, const Grid3D &grid) {
  for (Grid3D::array_type::index i = 0; i < grid.z; ++i) {
    out << "[";
    for (Grid3D::array_type::index j = 0; j < grid.y; ++j) {
      out << "[";
      for (Grid3D::array_type::index k = 0; k < grid.x; ++k) {
        out << grid.grid[i][j][k];
        // out << "x = " << k << ", y = " << j << ", z = " << i << " ";
      }
      out << "]\n";
    }
    out << "]\n";
  }
  // out << "Grid3D: " << grid.grid << std::endl;
  return out;
}
