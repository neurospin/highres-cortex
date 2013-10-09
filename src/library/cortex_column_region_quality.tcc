#include <yleprince/cortex_column_region_quality.hh>

#include <cmath>

namespace
{

static const float pi = 3.14159265358979323846f;

inline float square(float x)
{
  return x * x;
}

}

namespace yl
{

template <typename Tlabel>
float CortexColumnRegionQuality::
evaluate(const LabelVolume<Tlabel>& label_vol,
         const Tlabel label1,
         const Tlabel label2) const
{
  static const std::size_t canonical_region_size = 100;
  static const float region_size_variance = 30 * 30;
  static const float normalization_factor
    = 1.f / std::sqrt(2 * pi * region_size_variance);

  std::size_t region1_size = label_vol.region_size(label1);
  std::size_t region2_size = label_vol.region_size(label2);
  return normalization_factor * std::exp(-square(region1_size + region2_size - canonical_region_size) / region_size_variance);

#if 0 // TODO
  for(label_vol::const_point_iterator point_it = label_vol.begin(),
        point_end = label_vol.end();
      point_it != point_end;
      ++point_it)
  {
    const Point3d point = *point_it;
    const int x = point.X();
    const int y = point.Y();
    const int z = point.Z();

    
  }
#endif
}

}; // namespace yl
