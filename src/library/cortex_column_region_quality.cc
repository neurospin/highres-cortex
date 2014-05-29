#include <yleprince/cortex_column_region_quality.hh>


#include <yleprince/cortex_column_region_quality.tcc>

template yl::CortexColumnRegionQuality::Cache
yl::CortexColumnRegionQuality::cache<int32_t>(
  const LabelVolume<int32_t>&, int32_t) const;
template float yl::CortexColumnRegionQuality::evaluate<int32_t>(
  const LabelVolume<int32_t>&, int32_t) const;
template float yl::CortexColumnRegionQuality::evaluate<int32_t>(
  const LabelVolume<int32_t>&, int32_t, int32_t) const;


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
  const carto::Object& voxel_size = CSF_projections.header().getProperty("voxel_size");
  assert(voxel_size->isArray());
  m_sorted_voxel_sizes[0] = voxel_size->getArrayItem(0)->value<float>();
  m_sorted_voxel_sizes[1] = voxel_size->getArrayItem(1)->value<float>();
  m_sorted_voxel_sizes[2] = voxel_size->getArrayItem(2)->value<float>();
  std::sort(m_sorted_voxel_sizes, m_sorted_voxel_sizes + 3);
  m_pseudo_area_reliability_threshold
    = 0.5 * m_sorted_voxel_sizes[0] * m_sorted_voxel_sizes[1];
  setShapeParametres(default_goal_diametre(), default_max_thickness());
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
