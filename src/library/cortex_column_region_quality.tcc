#include <yleprince/cortex_column_region_quality.hh>

#include <cmath>
#include <limits>
#include <algorithm>

namespace
{

static const float pi = 3.14159265358979323846f;

inline float square(float x)
{
  return x * x;
}

inline float region_size_criterion(float region_size)
{
  static const float canonical_region_size = 200;
  static const float region_size_variance = 200 * 200;

  return std::exp(-square(region_size - canonical_region_size)
                  / region_size_variance);
}

inline float proximity_to_border_criterion(float smallest_distance)
{
  return std::exp(-std::max(smallest_distance, 0.f));
}

inline float barycentre_dist_criterion(float rms_mean2_dist)
{
  return std::exp(-rms_mean2_dist / 5);
}

template <class InputIterator>
Point3df region_barycentre(const InputIterator& first_point,
                           const InputIterator& last_point)
{
  Point3df sum(0.f, 0.f, 0.f);
  std::size_t num_points = 0;
  for(InputIterator point_it = first_point;
      point_it != last_point;
      ++point_it)
  {
    const Point3d& point = *point_it;
    const Point3df pointf(point[0], point[1], point[2]);

    sum += pointf;
    ++num_points;
  }
  sum /= static_cast<float>(num_points);
  return sum;
}

template <class BaseIterator>
class ChainedIterator
  : public boost::iterator_adaptor<
  ChainedIterator<BaseIterator>,
  BaseIterator,
  boost::use_default,
  boost::forward_traversal_tag>
{
public:
  ChainedIterator()
    : ChainedIterator::iterator_adaptor_(0) {}

  ChainedIterator(const BaseIterator& begin_first,
                  const BaseIterator& end_first,
                  const BaseIterator& begin_second)
    : ChainedIterator::iterator_adaptor_(begin_first),
      m_end_first(end_first),
      m_begin_second(begin_second) {}

  ChainedIterator(const BaseIterator& begin)
    : ChainedIterator::iterator_adaptor_(begin),
      m_end_first(),
      m_begin_second() {}


private:
  const BaseIterator m_end_first, m_begin_second;

  friend class boost::iterator_core_access;
  void increment()
  {
    BaseIterator& base = this->base_reference();
    ++base;
    if(base == m_end_first && m_end_first != BaseIterator())
      base = m_begin_second;
  }
};

} // end of anonymous namespace


namespace yl
{

template <class PointIterator>
inline float CortexColumnRegionQuality::
evaluate(const PointIterator& point_it_begin,
         const PointIterator& point_it_end) const
{
  const Point3df barycentre = region_barycentre(point_it_begin,
                                                point_it_end);

  std::size_t region_size = 0;
  float sum_dist2_to_bary = 0;
  float min_dist_CSF = std::numeric_limits<float>::infinity(),
    min_dist_white = std::numeric_limits<float>::infinity();
  for(PointIterator point_it = point_it_begin;
      point_it != point_it_end;
      ++point_it)
  {
    const Point3d& point = *point_it;
    const Point3df pointf(point[0], point[1], point[2]);
    const int x = point[0];
    const int y = point[1];
    const int z = point[2];

    ++region_size;

    const float dist_CSF = m_distmap_from_CSF.at(x, y, z);
    const float dist_white = m_distmap_from_white.at(x, y, z);

    if(dist_CSF < min_dist_CSF)
      min_dist_CSF = dist_CSF;
    if(dist_white < min_dist_white)
      min_dist_white = dist_white;

    sum_dist2_to_bary += (static_cast<Point3df>(pointf) - barycentre).norm2();
  }

  assert(std::isnormal(min_dist_CSF));
  assert(std::isnormal(min_dist_white));

  return region_size_criterion(region_size)
    * proximity_to_border_criterion(min_dist_CSF)
    * proximity_to_border_criterion(min_dist_white)
    * barycentre_dist_criterion(sum_dist2_to_bary / region_size);
}

template <typename Tlabel>
float CortexColumnRegionQuality::
evaluate(const LabelVolume<Tlabel>& label_vol, const Tlabel label) const
{
  return evaluate(label_vol.region_begin(label),
                  label_vol.region_end(label));
}
template <typename Tlabel>
float CortexColumnRegionQuality::
evaluate(const LabelVolume<Tlabel>& label_vol,
         const Tlabel label1,
         const Tlabel label2) const
{
  typedef ChainedIterator<typename LabelVolume<Tlabel>::const_point_iterator>
    ChainedPointIterator;
  return evaluate(ChainedPointIterator(label_vol.region_begin(label1),
                                       label_vol.region_end(label1),
                                       label_vol.region_begin(label2)),
                  ChainedPointIterator(label_vol.region_end(label2)));
}

}; // namespace yl
