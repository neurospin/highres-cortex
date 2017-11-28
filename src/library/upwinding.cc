/*
Copyright Télécom ParisTech (2015).

Contributor: Yann Leprince <yann.leprince@ylep.fr>.

This file is part of highres-cortex, a collection of software designed
to process high-resolution magnetic resonance images of the cerebral
cortex.

This software is governed by the CeCILL licence under French law and
abiding by the rules of distribution of free software. You can use,
modify and/or redistribute the software under the terms of the CeCILL
licence as circulated by CEA, CNRS and INRIA at the following URL:
<http://www.cecill.info/>.

As a counterpart to the access to the source code and rights to copy,
modify and redistribute granted by the licence, users are provided only
with a limited warranty and the software's author, the holder of the
economic rights, and the successive licensors have only limited
liability.

In this respect, the user's attention is drawn to the risks associated
with loading, using, modifying and/or developing or reproducing the
software by the user in light of its specific status of scientific
software, that may mean that it is complicated to manipulate, and that
also therefore means that it is reserved for developers and experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or
data to be ensured and, more generally, to use and operate it in the
same conditions as regards security.

The fact that you are presently reading this means that you have had
knowledge of the CeCILL licence and that you accept its terms.
*/

#include "upwinding.hh"

#include <algorithm>
#include <cmath>
#include <iostream>

#include <aims/connectivity/structuring_element.h>
#include <aims/vector/vector.h>

#include "front.hh"
#include "volume_util.hh"

#include "volume_util.tcc"

namespace
{

template<typename T>
inline T square(T value)
{
  return value * value;
}

inline int
upwind_direction(float field_left,
                 float field_centre,
                 float field_right)
{
  if(field_right < field_left && field_right < field_centre) {
    return 1;
  } else if(field_left < field_centre) {
    return -1;
  } else {
    return 0;
  }
}

/* depending on the compiler (gcc) version and standard (c++9x/ c++11),
   isnan may be a template or a regular function. I have found no other way
   than using a wrapping structure.
   (see also: http://stackoverflow.com/questions/17574242/how-to-use-isnan-as-a-predicate-function-to-stdfind-if-c11)
*/
struct float_isnan
{
  inline bool operator ()( float x ) const
  {
    return std::isnan( x );
  }
};

} // end of anonymous namespace

carto::VolumeRef<float>
yl::upwind_distance(const carto::VolumeRef<float>& upwind_field,
                    carto::VolumeRef<int16_t> domain,
                    const int16_t domain_label,
                    const int16_t origin_label,
                    const int16_t done_label,
                    const int16_t front_label)
{
  assert(domain_label != origin_label);
  assert(domain_label != front_label);
  assert(domain_label != done_label);
  assert(domain_label != front_label);
  assert(origin_label != front_label);
  assert(done_label != front_label);

  assert(yl::xyz_min_border(domain) >= 1);
  assert(yl::check_border_values(domain,
    std::bind1st(std::not_equal_to<int16_t>(), domain_label)));

  assert(yl::xyz_min_border(upwind_field) >= 1);
  // see above float_isnan struct for why not using std::isnan directly
  assert(yl::check_border_values(upwind_field, float_isnan()));

  std::auto_ptr<aims::strel::Connectivity> connectivity(
    aims::strel::ConnectivityFactory::create("6"));

  carto::VolumeRef<float> solution(domain.getSizeX(),
                                   domain.getSizeY(),
                                   domain.getSizeZ(),
                                   1, 1 /* border */);
  solution->copyHeaderFrom(domain.header());
  // Also fill the border
  solution->refVolume()->fill(std::numeric_limits<float>::quiet_NaN());

  std::clog << "initializing front..." << std::endl;
  // TODO better initialization (based on level set / halfway beween classes)
  SortedFront front(domain, front_label, done_label);
  {
    const int size_x = domain->getSizeX();
    const int size_y = domain->getSizeY();
    const int size_z = domain->getSizeZ();
    for(int z = 0 ; z < size_z ; ++z)
      for(int y = 0 ; y < size_y ; ++y)
        for(int x = 0 ; x < size_x ; ++x) {
          if(domain(x, y, z) == domain_label) {
            const Point3d point(x, y, z);
            bool add = false;

            for(aims::strel::Connectivity::const_iterator
                  neighbour_it = connectivity->begin() ;
                neighbour_it != connectivity->end() ;
                ++neighbour_it) {
              const Point3d& neighbour_offset = *neighbour_it;
              const Point3d neighbour_point = point + neighbour_offset;
              const int nx = neighbour_point[0];
              const int ny = neighbour_point[1];
              const int nz = neighbour_point[2];

              if(domain(nx, ny, nz) == origin_label) {
                solution(nx, ny, nz) = 0;
                add = true;
              }
            }
            if(add)
              front.add(point, upwind_field(x, y, z));
          }
        }
  }

  const std::vector<float> voxsize = domain->getVoxelSize();
  const float inv_voxsize_x = 1 / voxsize[0];
  const float inv_voxsize_y = 1 / voxsize[1];
  const float inv_voxsize_z = 1 / voxsize[2];

  unsigned int iter = 0;
  unsigned int num_nans = 0;

  std::clog << "upwinding..." << std::endl;

  while(!front.empty()) {
    ++iter;

    // the point with lowest value of upwind_field
    const Point3d point = front.pop();

    const int x = point[0];
    const int y = point[1];
    const int z = point[2];

    if(iter % 10000 == 0) { // verbose
      std::clog << "\r  iteration " << iter
                << ", current field " << upwind_field(x, y, z);
      if(front.constant_time_size)
        std::clog << ", " << front.size() << " voxels in front";
      std::clog << ".   " << std::flush;
    }

    float fx, fy, fz;
    float vx, vy, vz;
    int dir;
    dir = upwind_direction(upwind_field(x - 1, y, z),
                           upwind_field(x, y, z),
                           upwind_field(x + 1, y, z));
    if(dir) {
      fx = (upwind_field(x, y, z) - upwind_field(x + dir, y, z)) * inv_voxsize_x;
      vx = solution(x + dir, y, z);
    } else {
      fx = 0;
      vx = 0;
    }

    dir = upwind_direction(upwind_field(x, y - 1, z),
                           upwind_field(x, y, z),
                           upwind_field(x, y + 1, z));
    if(dir) {
      fy = (upwind_field(x, y, z) - upwind_field(x, y + dir, z)) * inv_voxsize_y;
      vy = solution(x, y + dir, z);
    } else {
      fy = 0;
      vy = 0;
    }

    dir = upwind_direction(upwind_field(x, y, z - 1),
                           upwind_field(x, y, z),
                           upwind_field(x, y, z + 1));
    if(dir) {
      fz = (upwind_field(x, y, z) - upwind_field(x, y, z + dir)) * inv_voxsize_z;
      vz = solution(x, y, z + dir);
    } else {
      fz = 0;
      vz = 0;
    }

    // Normalize the gradient and weigh the axes by their voxel sizes
    const float inv_fnorm = 1 / std::sqrt(square(fx) + square(fy) + square(fz));
    fx *= inv_fnorm * inv_voxsize_x;
    fy *= inv_fnorm * inv_voxsize_y;
    fz *= inv_fnorm * inv_voxsize_z;

    const float speed = 1;
    const float new_distance = (speed + (fx * vx + fy * vy + fz * vz))
      / (fx + fy + fz);

    if(!std::isnan(new_distance)) {
      solution(x, y, z) = new_distance;

      // add neighbours
      for(aims::strel::Connectivity::const_iterator
            neighbour_it = connectivity->begin() ;
          neighbour_it != connectivity->end() ;
          ++neighbour_it) {
        const Point3d& neighbour_offset = *neighbour_it;
        const Point3d neighbour_point = point + neighbour_offset;
        const int nx = neighbour_point[0];
        const int ny = neighbour_point[1];
        const int nz = neighbour_point[2];

        if(domain(nx, ny, nz) == domain_label) {
          front.add(neighbour_point, upwind_field(nx, ny, nz));
        }
      } // end for every neighbour
    } else {
      ++num_nans;
    }
  } // end of iteration on front

  // verbose
  std::clog << "\nfinished after " << iter << " iterations ("
            << num_nans << " NaN results)." << std::endl;

  return solution;
}
