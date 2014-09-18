#include <yleprince/cortex_column_region_quality.hh>


#include <yleprince/cortex_column_region_quality.tcc>

template yl::CortexColumnRegionQuality::Cache
yl::CortexColumnRegionQuality::cache<int32_t>(
  const LabelVolume<int32_t>&, int32_t) const;
template float yl::CortexColumnRegionQuality::evaluate<int32_t>(
  const LabelVolume<int32_t>&, int32_t) const;
template float yl::CortexColumnRegionQuality::evaluate<int32_t>(
  const LabelVolume<int32_t>&, int32_t, int32_t) const;


#include <algorithm>


using carto::VolumeRef;

namespace
{

float diametre_to_pseudo_area(float diametre)
{
  return 0.125f * square(diametre);
}

} // end of anonymous namespace

yl::CortexColumnRegionQuality::
CortexColumnRegionQuality(const VolumeRef<float>& CSF_projections,
                          const VolumeRef<float>& white_projections)
  : m_CSF_projections(CSF_projections),
    m_white_projections(white_projections)
{
  std::vector<float> voxel_size = CSF_projections->getVoxelSize();
  voxel_size.resize(3);
  std::sort(voxel_size.begin(), voxel_size.end());
  std::copy(voxel_size.begin(), voxel_size.end(), m_sorted_voxel_sizes);

  m_pseudo_area_reliability_threshold
    = 0.5 * m_sorted_voxel_sizes[0] * m_sorted_voxel_sizes[1];
  setShapeParametres(default_goal_diametre(), default_max_thickness());
  assert(m_CSF_projections.getSizeT() == 3);
  assert(m_white_projections.getSizeT() == 3);
}

void
yl::CortexColumnRegionQuality::
setShapeParametres(float goal_diametre, float max_thickness)
{
  m_pseudo_area_cutoff = diametre_to_pseudo_area(goal_diametre);

  // TODO this could be improved for non-isotropic voxels
  float voxel_volume = m_sorted_voxel_sizes[0]
    * m_sorted_voxel_sizes[1]
    * m_sorted_voxel_sizes[2];
  m_max_criterion = 2 * pi / voxel_volume * max_thickness;
}

float yl::CortexColumnRegionQuality::default_goal_diametre()
{
  return 0.5f;
}

float yl::CortexColumnRegionQuality::default_max_thickness()
{
  return 2.f;
}
