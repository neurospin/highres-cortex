#ifndef CORTEX_COLUMN_REGION_QUALITY_HH_INCLUDED
#define CORTEX_COLUMN_REGION_QUALITY_HH_INCLUDED

#include <cartodata/volume/volume.h>
#include <yleprince/label_volume.hh>

namespace yl
{

class CortexColumnRegionQuality
{
public:
  CortexColumnRegionQuality(const carto::VolumeRef<float>& CSF_projections,
                            const carto::VolumeRef<float>& white_projections);

  void setShapeParametres(float goal_diametre, float max_thickness);

  template <typename Tlabel>
  float evaluate(const LabelVolume<Tlabel>&, Tlabel) const;

  template <typename Tlabel>
  float evaluate(const LabelVolume<Tlabel>&, Tlabel, Tlabel) const;

  template <class PointIterator>
  inline float evaluate(const PointIterator& point_it_begin,
                        const PointIterator& point_it_end) const;

  static float default_goal_diametre();
  static float default_max_thickness();

private:
  carto::VolumeRef<float> m_CSF_projections;
  carto::VolumeRef<float> m_white_projections;
  float m_sorted_voxel_sizes[3];
  float m_pseudo_area_reliability_threshold;
  float m_pseudo_area_cutoff;
  float m_max_criterion;
}; // class CortexColumnRegionQuality

}; // namespace yl

#endif // !defined(CORTEX_COLUMN_REGION_QUALITY_HH_INCLUDED)
