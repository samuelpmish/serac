// Copyright (c) 2019-2020, Lawrence Livermore National Security, LLC and
// other Serac Project Developers. See the top-level LICENSE file for
// details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "mesh_utils.hpp"

#include <fstream>

#include "common/logger.hpp"
#include "common/terminator.hpp"
#include "fmt/fmt.hpp"

namespace serac {

std::shared_ptr<mfem::ParMesh> buildParallelMesh(const std::string& mesh_file, const int refine_serial,
                                                 const int refine_parallel, const MPI_Comm comm)
{
  // Get the MPI rank for logging purposes
  int rank = 0;
  MPI_Comm_rank(comm, &rank);

  // Open the mesh
  std::string msg = fmt::format("Opening mesh file: {0}", mesh_file);
  SLIC_INFO_ROOT(rank, msg);
  std::ifstream imesh(mesh_file);

  if (!imesh) {
    serac::logger::flush();
    std::string err_msg = fmt::format("Can not open mesh file: {0}", mesh_file);
    SLIC_ERROR_ROOT(rank, err_msg);
    serac::exitGracefully();
  }

  auto mesh = std::make_unique<mfem::Mesh>(imesh, 1, 1, true);
  imesh.close();

  // mesh refinement if specified in input
  for (int lev = 0; lev < refine_serial; lev++) {
    mesh->UniformRefinement();
  }

  // create the parallel mesh
  auto par_mesh = std::make_shared<mfem::ParMesh>(comm, *mesh);
  for (int lev = 0; lev < refine_parallel; lev++) {
    par_mesh->UniformRefinement();
  }

  return par_mesh;
}

void squish(mfem::Mesh& mesh)
{
  int num_vertices = mesh.GetNV();
  int dim          = mesh.SpaceDimension();

  mfem::Vector vertices;
  mesh.GetVertices(vertices);
  mfem::Vector vertex(dim);
  for (int i = 0; i < num_vertices; i++) {
    for (int d = 0; d < dim; d++) {
      vertex(d) = vertices[d * num_vertices + i];
    }

    double l1 = vertex.Norml1();
    double l2 = vertex.Norml2();
    vertex *= (l2 < 1.0e-6) ? 0.0 : (l1 / l2);

    for (int d = 0; d < dim; d++) {
      vertices[d * num_vertices + i] = vertex(d);
    }
  }
  mesh.SetVertices(vertices);
}

std::unique_ptr<mfem::Mesh> DiskMesh(int approx_number_of_elements)
{
  static constexpr int dim                   = 2;
  static constexpr int num_elems             = 4;
  static constexpr int num_vertices          = 5;
  static constexpr int num_boundary_elements = 4;

  static constexpr double vertices[num_vertices][dim] = {{0, 0}, {1, 0}, {0, 1}, {-1, 0}, {0, -1}};
  static constexpr int    triangles[num_elems][3]     = {{1, 2, 0}, {2, 3, 0}, {3, 4, 0}, {4, 1, 0}};
  static constexpr int    segments[num_elems][2]      = {{1, 2}, {2, 3}, {3, 4}, {4, 1}};

  auto mesh = std::make_unique<mfem::Mesh>(dim, num_vertices, num_elems, num_boundary_elements);

  for (auto vertex : vertices) {
    mesh->AddVertex(vertex);
  }
  for (auto triangle : triangles) {
    mesh->AddTriangle(triangle);
  }
  for (auto segment : segments) {
    mesh->AddBdrSegment(segment);
  }
  mesh->FinalizeTriMesh();

  while (mesh->GetNE() < (0.5 * approx_number_of_elements)) {
    mesh->UniformRefinement();
  }

  squish(*mesh);

  return mesh;
}

std::unique_ptr<mfem::Mesh> BallMesh(int approx_number_of_elements)
{
  static constexpr int dim                   = 3;
  static constexpr int num_elems             = 8;
  static constexpr int num_vertices          = 7;
  static constexpr int num_boundary_elements = 8;

  static constexpr double vertices[num_vertices][dim] = {{0, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, 0, -1},
                                                         {0, 0, 1}, {0, -1, 0}, {1, 0, 0}};
  static constexpr int    triangles[num_elems][3]     = {{4, 5, 6}, {4, 6, 2}, {4, 2, 1}, {4, 1, 5},
                                                  {5, 1, 3}, {5, 3, 6}, {3, 1, 2}, {6, 3, 2}};
  static constexpr int    tetrahedra[num_elems][4]    = {{0, 4, 5, 6}, {0, 4, 6, 2}, {0, 4, 2, 1}, {0, 4, 1, 5},
                                                   {0, 5, 1, 3}, {0, 5, 3, 6}, {0, 3, 1, 2}, {0, 6, 3, 2}};

  auto mesh = std::make_unique<mfem::Mesh>(dim, num_vertices, num_elems, num_boundary_elements);

  for (auto vertex : vertices) {
    mesh->AddVertex(vertex);
  }
  for (auto tetrahedron : tetrahedra) {
    mesh->AddTet(tetrahedron);
  }
  for (auto triangle : triangles) {
    mesh->AddBdrTriangle(triangle);
  }
  mesh->FinalizeTetMesh();

  while (mesh->GetNE() < (0.25 * approx_number_of_elements)) {
    mesh->UniformRefinement();
  }

  squish(*mesh);

  return mesh;
}

std::unique_ptr<mfem::Mesh> RectangleMesh(int elements_in_x, int elements_in_y)
{
  return std::make_unique<mfem::Mesh>(elements_in_x, elements_in_y, mfem::Element::QUADRILATERAL, true);
}

std::unique_ptr<mfem::Mesh> CuboidMesh(int elements_in_x, int elements_in_y, int elements_in_z)
{
  return std::make_unique<mfem::Mesh>(elements_in_x, elements_in_y, elements_in_z, mfem::Element::HEXAHEDRON, true);
}

}  // namespace serac
