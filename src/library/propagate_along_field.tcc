#include <yleprince/propagate_along_field.hh>

#include <cmath>
#include <iostream>
#include <stdexcept>

using carto::VolumeRef;
using std::clog;
using std::endl;

template<typename Tlabel>
Tlabel
yl::PropagateAlongField::
ascend_until_nonzero(const Point3df &start_point,
                     const VolumeRef<Tlabel> &seeds,
                     const Tlabel ignore_label) const
{
  if(debug_output >= 2 && m_verbose) {
    clog << "    ascend_until_nonzero at " << start_point << endl;
  }

  Point3df current_point = start_point;

  const int size_x = seeds.getSizeX();
  const int size_y = seeds.getSizeY();
  const int size_z = seeds.getSizeZ();

  unsigned int iter;
  for(iter = 0; iter < m_max_iter ; ++iter) {
    const float xp = current_point[0],
      yp = current_point[1], zp = current_point[2];

    const int ix = static_cast<int>((xp * m_invsize_x) + 0.5);
    if(ix < 0 || ix >= size_x) return 0;
    const int iy = static_cast<int>((yp * m_invsize_y) + 0.5);
    if(iy < 0 || iy >= size_y) return 0;
    const int iz = static_cast<int>((zp * m_invsize_z) + 0.5);
    if(iz < 0 || iz >= size_z) return 0;

    const Tlabel seed_value = seeds(ix, iy, iz);
    if(debug_output >= 2 && m_verbose) {
      clog << "      iteration " << iter << " at " << current_point
           << ", seed_value = " << seed_value << endl;
    }
    if(seed_value != 0 && seed_value != ignore_label) {
      return seed_value;
    }


    // Move along the field
    float gx = m_interp_fieldx.value(current_point);
    float gy = m_interp_fieldy.value(current_point);
    float gz = m_interp_fieldz.value(current_point);

    // Normalize the field, stop if too small or infinite on NaN.
    const float gn = std::sqrt(gx*gx + gy*gy + gz*gz);
    if(!std::isnormal(gn)) break;
    gx /= gn; gy /= gn; gz /= gn;

    current_point[0] = xp + m_step * gx;
    current_point[1] = yp + m_step * gy;
    current_point[2] = zp + m_step * gz;
  }

  if(m_verbose) {
    clog << "    ascension at " << start_point << " aborted after "
         << iter << " iterations" << endl;
  }
  return 0;
}

template<typename Tlabel>
VolumeRef<Tlabel>
yl::PropagateAlongField::
propagate_regions(const VolumeRef<Tlabel> &seeds,
                  const Tlabel target_label) const
{
  VolumeRef<Tlabel> regions(new carto::Volume<Tlabel>(*seeds));

  const int size_x = seeds.getSizeX();
  const int size_y = seeds.getSizeY();
  const int size_z = seeds.getSizeZ();

  if(m_verbose) {
    clog << "yl::PropagateAlongField::propagate_regions:\n"
         << "  maximum propagation distance: " << m_step * m_max_iter
         << " mm." << endl;
  }

  unsigned int n_propagated = 0, n_dead_end = 0, n_lost = 0;

  for(int z = 0; z < size_z; ++z)
  for(int y = 0; y < size_y; ++y)
  for(int x = 0; x < size_x; ++x)
  {
    if(m_verbose && x == 0 && y == 0) {
      clog << "  at slice " << z << " / " << size_z << ", "
           << n_propagated << " propagated, "
           << n_dead_end << " dead-end, "
           << n_lost << " lost."<< endl;
    }

    if(seeds(x, y, z) == target_label) {
      const Point3df point(x * m_voxel_size_x,
                           y * m_voxel_size_x,
                           z * m_voxel_size_x);
      Tlabel result = ascend_until_nonzero(point, seeds,
                                            target_label);
      if(result > 0) {
        regions(x, y, z) = result;
        ++n_propagated;
      } else if(result < 0) {
        ++n_dead_end;
      } else {  // result == 0
        ++n_lost;
      }
    }
  }

  return regions;
}
