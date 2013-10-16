#ifndef CORTEX_COLUMN_REGION_QUALITY_HH_INCLUDED
#define CORTEX_COLUMN_REGION_QUALITY_HH_INCLUDED

#include <cartodata/volume/volume.h>
#include <yleprince/label_volume.hh>

namespace yl
{

class CortexColumnRegionQuality
{
public:
  CortexColumnRegionQuality(const carto::VolumeRef<float>& distmap_from_CSF,
                            const carto::VolumeRef<float>& distmap_from_white);

  void setBorderProxWeight(const float weight)
  { m_border_prox_weight = weight; };
  void setCompacityWeight(const float weight)
  { m_compacity_weight = weight; };
  void setSizeWeight(const float weight)
  { m_size_weight = weight; };

  template <typename Tlabel>
  float evaluate(const LabelVolume<Tlabel>&, Tlabel) const;

  template <typename Tlabel>
  float evaluate(const LabelVolume<Tlabel>&, Tlabel, Tlabel) const;

  template <class PointIterator>
  inline float evaluate(const PointIterator& point_it_begin,
                        const PointIterator& point_it_end) const;

  static float default_border_prox_weight() { return 100.f; };
  static float default_compacity_weight() { return 1.f; };
  static float default_size_weight() { return 0.1f; };

private:
  carto::VolumeRef<float> m_distmap_from_CSF;
  carto::VolumeRef<float> m_distmap_from_white;
  float m_border_prox_weight;
  float m_compacity_weight;
  float m_size_weight;
}; // class CortexColumnRegionQuality

}; // namespace yl

#endif // !defined(CORTEX_COLUMN_REGION_QUALITY_HH_INCLUDED)
