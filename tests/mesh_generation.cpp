// Copyright (c) 2019-2020, Lawrence Livermore National Security, LLC and
// other Serac Project Developers. See the top-level LICENSE file for
// details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include <gtest/gtest.h>

#include "common/mesh_utils.hpp"

TEST(meshgen, successful_creation)
{
  // the disk and ball meshes don't exactly hit the number
  // of elements specified, they refine to get as close as possible
  ASSERT_EQ(serac::DiskMesh(1000)->GetNE(), 1024);
  ASSERT_EQ(serac::BallMesh(6000)->GetNE(), 4096);
  ASSERT_EQ(serac::RectangleMesh(20, 20)->GetNE(), 400);
  ASSERT_EQ(serac::CuboidMesh(20, 20, 20)->GetNE(), 8000);
}