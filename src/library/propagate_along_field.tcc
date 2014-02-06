#include <yleprince/propagate_along_field.hh>

#include <cmath>
#include <iostream>
#include <stdexcept>

#include <boost/shared_ptr.hpp>
#include <boost/tuple/tuple.hpp>

using carto::VolumeRef;
using std::clog;
using std::endl;

namespace
{

template <typename Tlabel>
struct SeedValueAscensionResult
{
  typedef Tlabel result_type;
  Tlabel operator()(const Tlabel& label, const Point3df&, float) const
  {
    return label;
  };
  static const Tlabel& failure_result()
  {
    static const Tlabel failure_res = 0;
    return failure_res;
  };
};
template <typename Tlabel>
struct SeedAndPointAscensionResult
{
  typedef std::pair<Tlabel, Point3df> result_type;
  std::pair<Tlabel, Point3df>
  operator()(const Tlabel& label, const Point3df& point, float) const
  {
    return std::make_pair(label, point);
  };
  static const std::pair<Tlabel, Point3df>& failure_result()
  {
    static const std::pair<Tlabel, Point3df> failure_res
      = std::make_pair(0, Point3df(-1, -1, -1));
    return failure_res;
  }
};

template <typename Tlabel>
struct SeedPointDistAscensionResult
{
  typedef boost::tuple<Tlabel, Point3df, float> result_type;

  result_type
  operator()(const Tlabel& label, const Point3df& point, float distance) const
  {
    return boost::make_tuple(label, point, distance);
  };
  static const result_type& failure_result()
  {
    static const result_type failure_res
      = boost::make_tuple(0, Point3df(-1, -1, -1), 0.0f);
    return failure_res;
  }
};


} // end of anonymous namespace

template <template <typename> class ResultChooserTemplate, typename Tlabel>
typename ResultChooserTemplate<Tlabel>::result_type
inline
yl::PropagateAlongField::
internal_ascension(const Point3df &start_point,
                   const VolumeRef<Tlabel> &seeds,
                   const Tlabel ignore_label,
                   const ResultChooserTemplate<Tlabel>& result_chooser) const
{
  return internal_ascension<ResultChooserTemplate<Tlabel> >
           (start_point, seeds, ignore_label, result_chooser);
}

template <class ResultChooser, typename Tlabel>
typename ResultChooser::result_type
yl::PropagateAlongField::
internal_ascension(const Point3df &start_point,
                   const VolumeRef<Tlabel> &seeds,
                   const Tlabel ignore_label,
                   const ResultChooser& result_chooser) const
{
  if(debug_output >= 3 && m_verbose >= 3) {
    clog << "    ascension at " << start_point << endl;
  }

  const carto::Object &voxel_size = seeds.header().getProperty("voxel_size");
  assert(voxel_size->isArray());
  const float voxel_size_x = voxel_size->getArrayItem(0)->value<float>();
  const float voxel_size_y = voxel_size->getArrayItem(1)->value<float>();
  const float voxel_size_z = voxel_size->getArrayItem(2)->value<float>();
  const float invsize_x = 1.0f / voxel_size_x;
  const float invsize_y = 1.0f / voxel_size_y;
  const float invsize_z = 1.0f / voxel_size_z;
  if(!std::isnormal(voxel_size_x) ||
     !std::isnormal(voxel_size_y) ||
     !std::isnormal(voxel_size_z)) {
    throw std::runtime_error("inconsistent voxel_size value");
  }

  Point3df current_point = start_point;

  const int size_x = seeds.getSizeX();
  const int size_y = seeds.getSizeY();
  const int size_z = seeds.getSizeZ();

  unsigned int iter;
  for(iter = 0; iter < m_max_iter ; ++iter) {
    const float xp = current_point[0],
      yp = current_point[1], zp = current_point[2];

    const int ix = static_cast<int>((xp * invsize_x) + 0.5f);
    if(ix < 0 || ix >= size_x) break;
    const int iy = static_cast<int>((yp * invsize_y) + 0.5f);
    if(iy < 0 || iy >= size_y) break;
    const int iz = static_cast<int>((zp * invsize_z) + 0.5f);
    if(iz < 0 || iz >= size_z) break;

    const Tlabel seed_value = seeds(ix, iy, iz);
    if(debug_output >= 4 && m_verbose >= 4) {
      clog << "      iteration " << iter << " at " << current_point
           << ", seed_value = " << seed_value << endl;
    }
    if(seed_value != 0 && seed_value != ignore_label) {
      return result_chooser(seed_value, current_point, std::abs(m_step * iter));
    }

    // Move along the field
    Point3df local_field;
    float& gx = local_field[0], gy = local_field[1], gz = local_field[2];
    try {
      m_vector_field->evaluate(current_point, local_field);
    } catch(const VectorField::UndefinedField&) {
      if(m_verbose >= 2) {
        clog << "    ascension at " << start_point << " aborted after "
             << iter << " iterations: vector field is undefined at ("
             << gx << ", " << gy << ", " << gz << ")" << endl;
      }
      return ResultChooser::failure_result();
    }

    // Normalize the field, stop if too small or infinite or NaN.
    const float gn = std::sqrt(gx*gx + gy*gy + gz*gz);
    if(!std::isnormal(gn)) break;
    gx /= gn; gy /= gn; gz /= gn;

    current_point[0] = xp + m_step * gx;
    current_point[1] = yp + m_step * gy;
    current_point[2] = zp + m_step * gz;
  }

  if(m_verbose >= 2) {
    clog << "    ascension at " << start_point << " aborted after "
         << iter << " iterations" << endl;
  }
  return ResultChooser::failure_result();
}


class ScalarFieldSumVisitor
{
public:
  ScalarFieldSumVisitor(const boost::shared_ptr<yl::ScalarField>& field)
    : m_field(field), m_sum(0.f) {};

  void init(const Point3df&) {};
  void visit(const Point3df& point)
  {
    m_sum += m_field->evaluate(point);
  };

  float sum() const
  { return m_sum; };
private:
  boost::shared_ptr<yl::ScalarField> m_field;
  float m_sum;
};

class UnitSurfaceEvolver
{
public:
  UnitSurfaceEvolver(const boost::shared_ptr<yl::ScalarField>& divergence_field)
    : m_divergence_field(divergence_field), m_surface(1.f) {};

  void init(const Point3df& point)
  {
    m_previous_point = point;
  };

  void visit(const Point3df& point)
  {
    const float step = (point - m_previous_point).norm();
    m_surface *= (1.f + step * m_divergence_field->evaluate(point));
  };

  float end_surface() const
  { return m_surface; };
private:
  boost::shared_ptr<yl::ScalarField> m_divergence_field;
  float m_surface;
  Point3df m_previous_point;
};

template <typename Tlabel>
float yl::PropagateAlongField::
evolve_unit_surface(const Point3df& start_point,
                    const VolumeRef<Tlabel>& seeds,
                    const boost::shared_ptr<yl::ScalarField>& divergence_field,
                    const Tlabel ignore_label) const
{
  UnitSurfaceEvolver visitor(divergence_field);
  if(visitor_ascension(start_point, seeds, ignore_label, visitor)) {
    return visitor.end_surface();
  } else {
    throw AbortedAdvection();
  }
}

template <typename Tlabel>
float yl::PropagateAlongField::
integrate_field_along_advection(const Point3df& start_point,
                                const VolumeRef<Tlabel>& seeds,
                                const boost::shared_ptr<yl::ScalarField>& field,
                                const Tlabel ignore_label) const
{
  ScalarFieldSumVisitor visitor(field);
  if(visitor_ascension(start_point, seeds, ignore_label, visitor)) {
    return m_step * visitor.sum();
  } else {
    throw AbortedAdvection();
  }
}

template <class TVisitor, typename Tlabel>
bool
yl::PropagateAlongField::
visitor_ascension(const Point3df& start_point,
                  const VolumeRef<Tlabel>& seeds,
                  const Tlabel ignore_label,
                  TVisitor& visitor) const
{
  if(debug_output >= 3 && m_verbose >= 3) {
    clog << "    ascension at " << start_point << endl;
  }

  const carto::Object &voxel_size = seeds.header().getProperty("voxel_size");
  assert(voxel_size->isArray());
  const float voxel_size_x = voxel_size->getArrayItem(0)->value<float>();
  const float voxel_size_y = voxel_size->getArrayItem(1)->value<float>();
  const float voxel_size_z = voxel_size->getArrayItem(2)->value<float>();
  const float invsize_x = 1.0f / voxel_size_x;
  const float invsize_y = 1.0f / voxel_size_y;
  const float invsize_z = 1.0f / voxel_size_z;
  if(!std::isnormal(voxel_size_x) ||
     !std::isnormal(voxel_size_y) ||
     !std::isnormal(voxel_size_z)) {
    throw std::runtime_error("inconsistent voxel_size value");
  }

  Point3df current_point = start_point;
  visitor.init(start_point);

  const int size_x = seeds.getSizeX();
  const int size_y = seeds.getSizeY();
  const int size_z = seeds.getSizeZ();

  unsigned int iter;
  for(iter = 0; iter < m_max_iter ; ++iter) {
    const float xp = current_point[0],
      yp = current_point[1], zp = current_point[2];

    const int ix = static_cast<int>((xp * invsize_x) + 0.5f);
    if(ix < 0 || ix >= size_x) break;
    const int iy = static_cast<int>((yp * invsize_y) + 0.5f);
    if(iy < 0 || iy >= size_y) break;
    const int iz = static_cast<int>((zp * invsize_z) + 0.5f);
    if(iz < 0 || iz >= size_z) break;

    const Tlabel seed_value = seeds(ix, iy, iz);
    if(debug_output >= 4 && m_verbose >= 4) {
      clog << "      iteration " << iter << " at " << current_point
           << ", seed_value = " << seed_value << endl;
    }
    visitor.visit(current_point);
    if(seed_value != 0 && seed_value != ignore_label) {
      return true;
    }

    // Move along the field
    Point3df local_field;
    float& gx = local_field[0], gy = local_field[1], gz = local_field[2];
    try {
      m_vector_field->evaluate(current_point, local_field);
    } catch(const VectorField::UndefinedField&) {
      if(m_verbose >= 2) {
        clog << "    ascension at " << start_point << " aborted after "
             << iter << " iterations: vector field is undefined at ("
             << gx << ", " << gy << ", " << gz << ")" << endl;
      }
      return false;
    }

    // Normalize the field, stop if too small or infinite or NaN.
    const float gn = std::sqrt(gx*gx + gy*gy + gz*gz);
    if(!std::isnormal(gn)) break;
    gx /= gn; gy /= gn; gz /= gn;

    current_point[0] = xp + m_step * gx;
    current_point[1] = yp + m_step * gy;
    current_point[2] = zp + m_step * gz;
  }

  if(m_verbose >= 2) {
    clog << "    ascension at " << start_point << " aborted after "
         << iter << " iterations" << endl;
  }
  return false;
}


template <typename Tlabel>
Tlabel
yl::PropagateAlongField::
ascend_until_nonzero(const Point3df &start_point,
                     const VolumeRef<Tlabel> &seeds,
                     const Tlabel ignore_label) const
{
  return internal_ascension(start_point, seeds, ignore_label,
                            SeedValueAscensionResult<Tlabel>());
}


namespace {

template <typename Tlabel>
class SeedValueRecorder
{
public:
  SeedValueRecorder(const VolumeRef<Tlabel>& seeds)
    : m_seed_values(new carto::Volume<Tlabel>(*seeds)) {};

  void record(const int x, const int y, const int z,
              const std::pair<Tlabel, Point3df>& value)
  {
    m_seed_values(x, y, z) = value.first;
  };

  VolumeRef<Tlabel> result() const {return m_seed_values;};

private:
  VolumeRef<Tlabel> m_seed_values;
};

template <typename Tlabel>
class SeedValueAndPointRecorder
{
public:
  SeedValueAndPointRecorder(const VolumeRef<Tlabel>& seeds)
    : m_seed_values(new carto::Volume<Tlabel>(*seeds)),
      m_points(seeds.getSizeX(), seeds.getSizeY(), seeds.getSizeZ(), 3)
  {
    m_points.header() = seeds.header();
    m_points.fill(-1);
  };

  void record(const int x, const int y, const int z,
              const std::pair<Tlabel, Point3df>& value)
  {
    m_seed_values(x, y, z) = value.first;
    const Point3df& point = value.second;
    m_points(x, y, z, 0) = point[0];
    m_points(x, y, z, 1) = point[1];
    m_points(x, y, z, 2) = point[2];
  };

  std::pair<VolumeRef<Tlabel>, VolumeRef<float> >
  result() const {return std::make_pair(m_seed_values, m_points);};

private:
  VolumeRef<Tlabel> m_seed_values;
  VolumeRef<float> m_points;
};

} // end of anonymous namespace

template<typename Tlabel>
VolumeRef<Tlabel>
yl::PropagateAlongField::
propagate_regions(const VolumeRef<Tlabel> &seeds,
                  const Tlabel target_label) const
{
  SeedValueRecorder<Tlabel> result_recorder(seeds);
  internal_propagation(seeds, target_label, result_recorder);
  return result_recorder.result();
}

template<typename Tlabel>
std::pair<carto::VolumeRef<Tlabel>, carto::VolumeRef<float> >
yl::PropagateAlongField::
propagate_regions_keeping_dests(const VolumeRef<Tlabel> &seeds,
                                const Tlabel target_label) const
{
  SeedValueAndPointRecorder<Tlabel> result_recorder(seeds);
  internal_propagation(seeds, target_label, result_recorder);
  return result_recorder.result();
}

template<class ResultRecorder, typename Tlabel>
void
yl::PropagateAlongField::
internal_propagation(const VolumeRef<Tlabel> &seeds,
                     const Tlabel target_label,
                     ResultRecorder& result_recorder) const
{
  SeedAndPointAscensionResult<Tlabel> ascension_result_chooser;

  const int size_x = seeds.getSizeX();
  const int size_y = seeds.getSizeY();
  const int size_z = seeds.getSizeZ();

  const carto::Object &voxel_size = seeds.header().getProperty("voxel_size");
  assert(voxel_size->isArray());
  const float voxel_size_x = voxel_size->getArrayItem(0)->value<float>();
  const float voxel_size_y = voxel_size->getArrayItem(1)->value<float>();
  const float voxel_size_z = voxel_size->getArrayItem(2)->value<float>();

  if(m_verbose) {
    clog << "yl::PropagateAlongField::propagate_regions:\n"
            "  maximum propagation distance: " << m_step * m_max_iter
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
      const Point3df point(x * voxel_size_x,
                           y * voxel_size_y,
                           z * voxel_size_z);

      std::pair<Tlabel, Point3df> result
        = internal_ascension(point, seeds, target_label,
                             ascension_result_chooser);
      const Tlabel& result_label = result.first;

      if(result_label > 0) {
        result_recorder.record(x, y, z, result);
        ++n_propagated;
      } else if(result_label < 0) {
        ++n_dead_end;
      } else {  // result_label == 0
        ++n_lost;
      }
    }
  }

  if(m_verbose) {
    clog << "End of yl::PropagateAlongField::propagate_regions: "
         << n_propagated << " propagated, "
         << n_dead_end << " dead-end, "
         << n_lost << " lost."<< endl;
  }
}

template<typename Tlabel>
VolumeRef<float>
yl::PropagateAlongField::
evolve_unit_surface_from_region(const VolumeRef<Tlabel> &seeds,
                                const boost::shared_ptr<ScalarField>& divergence_field,
                                const Tlabel target_label) const
{
  const int size_x = seeds.getSizeX();
  const int size_y = seeds.getSizeY();
  const int size_z = seeds.getSizeZ();

  carto::VolumeRef<float> result_volume(size_x, size_y, size_z);
  result_volume.header() = seeds.header();
  static const float aborted_result = -1.f;
  static const float no_calculation_result = -2.f;
  result_volume.fill(no_calculation_result);

  const carto::Object &voxel_size = seeds.header().getProperty("voxel_size");
  assert(voxel_size->isArray());
  const float voxel_size_x = voxel_size->getArrayItem(0)->value<float>();
  const float voxel_size_y = voxel_size->getArrayItem(1)->value<float>();
  const float voxel_size_z = voxel_size->getArrayItem(2)->value<float>();

  if(m_verbose) {
    clog << "yl::PropagateAlongField::evolve_unit_surface_from_region:\n"
            "  maximum propagation distance: " << m_step * m_max_iter
         << " mm." << endl;
  }

  unsigned int n_propagated = 0, n_aborted = 0;

  for(int z = 0; z < size_z; ++z)
  for(int y = 0; y < size_y; ++y)
  for(int x = 0; x < size_x; ++x)
  {
    if(m_verbose && x == 0 && y == 0) {
      clog << "  at slice " << z << " / " << size_z << ", "
           << n_propagated << " propagated, "
           << n_aborted << " aborted." << endl;
    }

    if(seeds(x, y, z) == target_label) {
      const Point3df point(x * voxel_size_x,
                           y * voxel_size_y,
                           z * voxel_size_z);

      float local_result;
      try {
        local_result = evolve_unit_surface(point, seeds, divergence_field, target_label);
      } catch(const AbortedAdvection&) {
        ++n_aborted;
        result_volume(x, y, z) = aborted_result;
        continue;
      }

      result_volume(x, y, z) = local_result;
      ++n_propagated;
    }
  }

  if(m_verbose) {
    clog << "End of yl::PropagateAlongField::integrate_from_regions: "
         << n_propagated << " propagated, "
         << n_aborted << " aborted." << endl;
  }

  return result_volume;
}
