

#include <iostream>
#include <utility>
#include <vector>

#include "mfem.hpp"

using namespace mfem;

class EntitySet {
public:
  virtual Array<int> getDofs() const = 0;
};

class EdgeSet : EntitySet {
public:
  template <typename Pred>
  EdgeSet(const FiniteElementSpace& space, Pred pred)
  {
    const Mesh& mesh = *space.GetMesh();
    Array<int>  edges;
    Array<int>  edge_indices;
    Array<int>  orientation;

    for (int i = 0; i < mesh.GetNE(); i++) {
      auto ele = mesh.GetElement(i);
      if (pred(*ele)) {
        mesh.GetElementEdges(i, edges, orientation);
        edge_indices.Append(edges);
      }
    }

    Array<int> edge_dofs;
    for (const auto edge_index : edge_indices) {
      space.GetEdgeDofs(edge_index, edge_dofs);
      dofs.Append(edge_dofs);
    }
  }

  Array<int> getDofs() const override { return dofs; }

private:
  Array<int> dofs;
};

int main()
{
  constexpr int     n_divisions = 10;
  std::vector<Mesh> meshes;
  meshes.emplace_back(n_divisions);
  meshes.emplace_back(n_divisions, n_divisions, Element::Type::TRIANGLE);
  meshes.emplace_back(n_divisions, n_divisions, n_divisions, Element::Type::HEXAHEDRON);
  auto fec = H1_FECollection(2);

  for (auto&& mesh : meshes) {
    std::cout << "----MESH-----" << std::endl;
    std::cout << "num verts " << mesh.GetNV() << std::endl;
    std::cout << "num eles " << mesh.GetNE() << std::endl;
    std::cout << "num bound eles " << mesh.GetNBE() << std::endl;
    std::cout << "num edges " << mesh.GetNEdges() << std::endl;
    std::cout << "num faces " << mesh.GetNFaces() << std::endl;
    std::cout << "num face elements " << mesh.GetNumFaces() << std::endl;
    auto space   = FiniteElementSpace(&mesh, &fec);
    auto edgeset = EdgeSet(space, [](const Element& ele) { return true; });
  }
}
