#pragma ONCE
#include <array>
#include <cstddef>
#include <ostream>
#include <cmath>

class Spin{
    public:
        Spin(size_t dim) : dim{dim} {};

    protected:
        size_t dim;
};

class Spin3D : Spin{
    public:
    
    Spin3D(): Spin{3}{
        spin = {0, 0};
    };

    Spin3D(double theta, double phi): Spin{3}{
        spin = {theta, phi};
    };

    void update_coords(double theta, double phi){
        spin[0] = theta;
        spin[1] = phi;
    }
    friend std::ostream& operator << (std::ostream &out, const Spin3D &spin);

    /* Dot Product between two vectors, yields cos(omega) with omega, the angle in between*/
    double const operator*(const Spin3D &rhs){
        return std::sin(spin[0]) * std::sin(rhs.spin[0]) * std::cos(spin[1] - rhs.spin[1]) + std::cos(spin[0]) * std::cos(rhs.spin[0]);
    }

    Spin3D  operator+= (const Spin3D &rhs){
        
        this->spin[0] += rhs.spin[0];
        this->spin[1] += rhs.spin[1];
        return *this;
    }

    

        std::array<double, 2> spin;
    protected:
};

std::ostream& operator << (std::ostream &out, const Spin3D &spin){
    
        out << "Spin3D: " << spin.spin[0] << ", " << spin.spin[1] << " ";
        return out;
    }