

#include <gtest/gtest.h>

#include <iostream>
#include <utility>
#include <vector>

#include "common/common.hpp"
#include "mfem.hpp"

namespace serac {

TEST(entity_set, all_entities_1d)
{
  constexpr int DIVISIONS = 10;
  // Double as each element has two vertices
  constexpr int EXPECTED_VERT_DOFS = DIVISIONS * 2;

  mfem::Mesh               mesh(DIVISIONS);
  mfem::H1_FECollection    fec(1, 1);
  mfem::FiniteElementSpace space(&mesh, &fec);

  EntitySet all_vertices(space, EntitySet::Vertex, [](const mfem::Element&) { return true; });
  // There are no edges or faces in 1D

  EXPECT_EQ(EXPECTED_VERT_DOFS, all_vertices.dofs().Size());
}

TEST(entity_set, all_entities_2d)
{
  constexpr int DIVISIONS = 10;
  // Double for two triangles per quad, triple as each element (triangle) has three vertices
  constexpr int EXPECTED_VERT_DOFS = 2 * DIVISIONS * DIVISIONS * 3;
  constexpr int EXPECTED_EDGE_DOFS = 2 * DIVISIONS * DIVISIONS * 3 * 2;

  mfem::Mesh               mesh(DIVISIONS, DIVISIONS, mfem::Element::TRIANGLE);
  mfem::H1_FECollection    fec(1, 2);
  mfem::FiniteElementSpace space(&mesh, &fec);

  EntitySet all_vertices(space, EntitySet::Vertex, [](const mfem::Element&) { return true; });
  EntitySet all_edges(space, EntitySet::Edge, [](const mfem::Element&) { return true; });
  // There are no faces in 2D

  EXPECT_EQ(EXPECTED_VERT_DOFS, all_vertices.dofs().Size());
  EXPECT_EQ(EXPECTED_EDGE_DOFS, all_edges.dofs().Size());
}

TEST(entity_set, all_entities_3d)
{
  constexpr int DIVISIONS = 10;
  // Each element (hex/cube) has eight vertices, twelve edges, six faces (but shared)
  constexpr int EXPECTED_VERT_DOFS = DIVISIONS * DIVISIONS * DIVISIONS * 8;
  constexpr int EXPECTED_EDGE_DOFS = DIVISIONS * DIVISIONS * DIVISIONS * 12 * 2;
  constexpr int EXPECTED_FACE_DOFS = DIVISIONS * DIVISIONS * DIVISIONS * 6 * 4;

  mfem::Mesh               mesh(DIVISIONS, DIVISIONS, DIVISIONS, mfem::Element::HEXAHEDRON);
  mfem::H1_FECollection    fec(1, 2);
  mfem::FiniteElementSpace space(&mesh, &fec);

  EntitySet all_vertices(space, EntitySet::Vertex, [](const mfem::Element&) { return true; });
  EntitySet all_edges(space, EntitySet::Edge, [](const mfem::Element&) { return true; });
  EntitySet all_faces(space, EntitySet::Face, [](const mfem::Element&) { return true; });

  EXPECT_EQ(EXPECTED_VERT_DOFS, all_vertices.dofs().Size());
  EXPECT_EQ(EXPECTED_EDGE_DOFS, all_edges.dofs().Size());
  EXPECT_EQ(EXPECTED_FACE_DOFS, all_faces.dofs().Size());
}

}  // namespace serac
