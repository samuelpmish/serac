#include "common/mesh_utils.hpp"

#include <fstream>

int main() {
  {
    std::ofstream outfile("disk.mesh");
    serac::DiskMesh(1000)->Print(outfile);
    outfile.close();
  }

  {
    std::ofstream outfile("ball.mesh");
    serac::BallMesh(6000)->Print(outfile);
    outfile.close();
  }

  {
    std::ofstream outfile("square.mesh");
    serac::CuboidMesh(20, 20)->Print(outfile);
    outfile.close();
  }

  {
    std::ofstream outfile("cube.mesh");
    serac::CuboidMesh(20, 20, 20)->Print(outfile);
    outfile.close();
  }
}