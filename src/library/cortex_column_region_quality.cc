#include <yleprince/cortex_column_region_quality.hh>

using carto::VolumeRef;

yl::CortexColumnRegionQuality::
CortexColumnRegionQuality(const VolumeRef<float>& distmap_from_CSF,
                          const VolumeRef<float>& distmap_from_white)
  : m_distmap_from_CSF(distmap_from_CSF),
    m_distmap_from_white(distmap_from_white)
{
}

#include <yleprince/cortex_column_region_quality.tcc>

template float yl::CortexColumnRegionQuality::evaluate<int32_t>(
  const LabelVolume<int32_t>&, int32_t) const;

template float yl::CortexColumnRegionQuality::evaluate<int32_t>(
  const LabelVolume<int32_t>&, int32_t, int32_t) const;
