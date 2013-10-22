#include <yleprince/propagate_along_field.hh>

#include <cmath>
#include <iostream>
#include <stdexcept>

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

  Point3df current_point = start_point;

  const int size_x = seeds.getSizeX();
  const int size_y = seeds.getSizeY();
  const int size_z = seeds.getSizeZ();

  unsigned int iter;
  for(iter = 0; iter < m_max_iter ; ++iter) {
    const float xp = current_point[0],
      yp = current_point[1], zp = current_point[2];

    const int ix = static_cast<int>((xp * m_invsize_x) + 0.5f);
    if(ix < 0 || ix >= size_x) break;
    const int iy = static_cast<int>((yp * m_invsize_y) + 0.5f);
    if(iy < 0 || iy >= size_y) break;
    const int iz = static_cast<int>((zp * m_invsize_z) + 0.5f);
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
    float gx = m_interp_fieldx.value(current_point);
    float gy = m_interp_fieldy.value(current_point);
    float gz = m_interp_fieldz.value(current_point);

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
      const Point3df point(x * m_voxel_size_x,
                           y * m_voxel_size_y,
                           z * m_voxel_size_z);

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
