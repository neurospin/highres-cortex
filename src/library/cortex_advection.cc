/*
Copyright CEA (2014).
Copyright Université Paris XI (2014).

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

#include "cortex_advection.hh"

#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <utility>

#include <cartodata/volume/volume.h>
#include <aims/utility/converter_volume.h>

#include "advection.hh"
#include "field.hh"

using carto::VolumeRef;
using namespace std;

// Anonymous namespace for file-local symbols
namespace {

const float no_value = std::numeric_limits<float>::quiet_NaN();

class TubeAdvection : public yl::Advection::Visitor
{
public:
  TubeAdvection(const yl::ScalarField& divergence_field,
                const yl::ScalarField& domain,
                const bool opposite_direction=false)
    : m_divergence_field(divergence_field),
      m_domain(domain),
      m_opposite_direction(opposite_direction),
      m_previous_point(),
      m_surface(1.f),
      m_volume(0.f),
      m_abort(false)
  {
  };

  void first(const Point3df& point)
  {
    m_previous_point = point;
  };

  void visit(const Point3df& point)
  {
    const float step = (point - m_previous_point).norm();
    const float old_surface = m_surface;
    float divergence_value;
    try {
      divergence_value = m_divergence_field.evaluate(point);
    } catch(const yl::Field::UndefinedField&) {
      m_abort = true;
      return;
    }

    if(m_opposite_direction)
      m_surface *= 1 - step * divergence_value;
    else
      m_surface *= 1 + step * divergence_value;

    if(m_surface < 0) {
      clog << "  Warning: TubeAdvection encountered negative surface, aborting" << endl;
      m_abort = true;
      return;
    }
    m_volume += step * (old_surface + m_surface) / 2;
    m_previous_point = point;
  };

  bool move_on(const Point3df& point) const
  {
    if(m_abort)
      return false;
    try {
      return m_domain.evaluate(point) >= 0.5f;
    } catch(const yl::Field::UndefinedField&) {
      return false;
    }
  };

  float surface() const
  {
    if(m_abort)
      return no_value;
    else
      return m_surface;
  };

  float volume() const
  {
    if(m_abort)
      return no_value;
    else
      return m_volume;
  };

private:
  const yl::ScalarField& m_divergence_field;
  const yl::ScalarField& m_domain;
  const bool m_opposite_direction;
  Point3df m_previous_point;
  float m_surface;
  float m_volume;
  bool m_abort;
};

class EuclideanAdvection : public yl::Advection::Visitor
{
public:
  EuclideanAdvection(const yl::ScalarField& domain)
    : m_domain(domain),
      m_previous_point(),
      m_length(0.f),
      m_abort(false)
  {
  };

  void first(const Point3df& point)
  {
    m_previous_point = point;
  };

  void visit(const Point3df& point)
  {
    const float step = (point - m_previous_point).norm();
    m_length += step;
    m_previous_point = point;
  };

  bool move_on(const Point3df& point) const
  {
    if(m_abort)
      return false;
    try {
      return m_domain.evaluate(point) >= 0.5f;
    } catch(const yl::Field::UndefinedField&) {
      return false;
    }
  };

  float length() const
  {
    if(m_abort)
      return no_value;
    else
      return m_length;
  };

private:
  const yl::ScalarField& m_domain;
  Point3df m_previous_point;
  float m_length;
  bool m_abort;
};

} // end of anonymous namespace

std::pair<VolumeRef<float>, VolumeRef<float> >
yl::advect_tubes(const yl::VectorField3d& advection_field,
                 const yl::ScalarField& divergence_field,
                 const VolumeRef<int16_t>& domain,
                 const float max_advection_distance,
                 const float step_size,
                 const int verbosity)
{
  assert(max_advection_distance > 0);

  const int size_x = domain.getSizeX();
  const int size_y = domain.getSizeY();
  const int size_z = domain.getSizeZ();

  const std::vector<float> voxel_size = domain->getVoxelSize();
  const float voxel_size_x = voxel_size[0];
  const float voxel_size_y = voxel_size[1];
  const float voxel_size_z = voxel_size[2];

  carto::VolumeRef<float> surface_result(size_x, size_y, size_z);
  surface_result->copyHeaderFrom(domain.header());
  surface_result.fill(no_value);
  carto::VolumeRef<float> volume_result(size_x, size_y, size_z);
  volume_result->copyHeaderFrom(domain.header());
  volume_result.fill(no_value);

  unsigned int n_success = 0, n_aborted = 0;

  yl::ConstantStepAdvection advection(advection_field, step_size);
  advection.set_max_iter(std::ceil(max_advection_distance
                                   / std::abs(step_size)));
  advection.set_verbose(verbosity - 1);

  // This could be more elegant: the domain is first converted as float, then
  // fed into a scalar field to ease interpolation.
  carto::Converter<VolumeRef<int16_t>, VolumeRef<float> > conv;
  carto::VolumeRef<float> float_domain(*conv(domain));
  yl::LinearlyInterpolatedScalarField domain_field(float_domain);

  const bool opposite_direction = step_size < 0;

  for(int z = 0; z < size_z; ++z)
  for(int y = 0; y < size_y; ++y)
  for(int x = 0; x < size_x; ++x)
  {
    if(verbosity && x == 0 && y == 0) {
      clog << "\r  at slice " << z << " / " << size_z << ", "
           << n_success << " successfully advected, "
           << n_aborted << " aborted." << flush;
    }

    if(domain(x, y, z)) {
      const Point3df point(x * voxel_size_x,
                           y * voxel_size_y,
                           z * voxel_size_z);

      TubeAdvection visitor(divergence_field, domain_field, opposite_direction);
      yl::Advection::Visitor& plain_visitor = visitor;
      const bool success = advection.visitor_advection(plain_visitor, point);

      if(success) {
        volume_result(x, y, z) = visitor.volume();
        surface_result(x, y, z) = visitor.surface();
        ++n_success;
      } else {
        ++n_aborted;
      }
   }
  }

  if(verbosity)
    clog << "\ryl::advect_unit_surface: "
         << n_success << " propagated, "
         << n_aborted << " aborted." << endl;

  return std::make_pair(volume_result, surface_result);
}

VolumeRef<float>
yl::advect_euclidean(const yl::VectorField3d& advection_field,
                     const VolumeRef<int16_t>& domain,
                     const float max_advection_distance,
                     const float step_size,
                     const int verbosity)
{
  assert(max_advection_distance > 0);

  const int size_x = domain.getSizeX();
  const int size_y = domain.getSizeY();
  const int size_z = domain.getSizeZ();

  const std::vector<float> voxel_size = domain->getVoxelSize();
  const float voxel_size_x = voxel_size[0];
  const float voxel_size_y = voxel_size[1];
  const float voxel_size_z = voxel_size[2];

  carto::VolumeRef<float> length_result(size_x, size_y, size_z);
  length_result->copyHeaderFrom(domain.header());
  length_result.fill(no_value);

  unsigned int n_success = 0, n_aborted = 0;

  yl::ConstantStepAdvection advection(advection_field, step_size);
  advection.set_max_iter(std::ceil(max_advection_distance
                                   / std::abs(step_size)));
  advection.set_verbose(verbosity - 1);

  // This could be more elegant: the domain is first converted as float, then
  // fed into a scalar field to ease interpolation.
  carto::Converter<VolumeRef<int16_t>, VolumeRef<float> > conv;
  carto::VolumeRef<float> float_domain(*conv(domain));
  yl::LinearlyInterpolatedScalarField domain_field(float_domain);

  for(int z = 0; z < size_z; ++z)
  for(int y = 0; y < size_y; ++y)
  for(int x = 0; x < size_x; ++x)
  {
    if(verbosity && x == 0 && y == 0) {
      clog << "\r  at slice " << z << " / " << size_z << ", "
           << n_success << " successfully advected, "
           << n_aborted << " aborted." << flush;
    }

    if(domain(x, y, z)) {
      const Point3df point(x * voxel_size_x,
                           y * voxel_size_y,
                           z * voxel_size_z);

      EuclideanAdvection visitor(domain_field);
      yl::Advection::Visitor& plain_visitor = visitor;
      const bool success = advection.visitor_advection(plain_visitor, point);

      if(success) {
        length_result(x, y, z) = visitor.length();
        ++n_success;
      } else {
        ++n_aborted;
      }
   }
  }

  if(verbosity)
    clog << "\ryl::advect_unit_surface: "
         << n_success << " propagated, "
         << n_aborted << " aborted." << endl;

  return length_result;
}
