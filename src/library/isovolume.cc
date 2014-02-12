#include <yleprince/isovolume.hh>

#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <utility>

#include <cartodata/volume/volume.h>
#include <aims/utility/converter_volume.h>

#include <yleprince/advection.hh>
#include <yleprince/field.hh>

using carto::VolumeRef;
using namespace std;

// Anonymous namespace for file-local symbols
namespace {

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

  static const float no_value;

private:
  const yl::ScalarField& m_divergence_field;
  const yl::ScalarField& m_domain;
  const bool m_opposite_direction;
  Point3df m_previous_point;
  float m_surface;
  float m_volume;
  bool m_abort;
};

const float TubeAdvection::no_value = std::numeric_limits<float>::quiet_NaN();

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

  const carto::Object& voxel_size = domain.header().getProperty("voxel_size");
  assert(voxel_size->isArray());
  const float voxel_size_x = voxel_size->getArrayItem(0)->value<float>();
  const float voxel_size_y = voxel_size->getArrayItem(1)->value<float>();
  const float voxel_size_z = voxel_size->getArrayItem(2)->value<float>();

  carto::VolumeRef<float> surface_result(size_x, size_y, size_z);
  surface_result.header() = domain.header();
  surface_result.fill(TubeAdvection::no_value);
  carto::VolumeRef<float> volume_result(size_x, size_y, size_z);
  volume_result.header() = domain.header();
  volume_result.fill(TubeAdvection::no_value);

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
           << n_success << " succesfully advected, "
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
