#include <yleprince/propagate_along_field.hh>

#include <cmath>
#include <stdexcept>

using carto::VolumeRef;
using std::clog;
using std::endl;

namespace {
  const int debug_output = 0;
}

yl::PropagateAlongField::
PropagateAlongField(const VolumeRef<float> &fieldx,
                    const VolumeRef<float> &fieldy,
                    const VolumeRef<float> &fieldz)
  : m_interp_fieldx(fieldx), m_interp_fieldy(fieldy), m_interp_fieldz(fieldz),
    m_max_iter(default_max_iter), m_step(default_step), m_verbose(debug_output)
{
  const carto::Object &voxel_size = fieldx.header().getProperty("voxel_size");
  assert(voxel_size->isArray());
  m_voxel_size_x = voxel_size->getArrayItem(0)->value<float>();
  m_voxel_size_y = voxel_size->getArrayItem(1)->value<float>();
  m_voxel_size_z = voxel_size->getArrayItem(2)->value<float>();
  m_invsize_x = 1.0f / m_voxel_size_x;
  m_invsize_y = 1.0f / m_voxel_size_y;
  m_invsize_z = 1.0f / m_voxel_size_z;
  if(!std::isnormal(m_voxel_size_x) ||
     !std::isnormal(m_voxel_size_y) ||
     !std::isnormal(m_voxel_size_z)) {
    throw std::runtime_error("inconsistent voxel_size value");
  }
}

yl::PropagateAlongField::~PropagateAlongField()
{
}

void
yl::PropagateAlongField::setVerbose(bool verbose)
{
  m_verbose = verbose;
}

void
yl::PropagateAlongField::setStep(float step)
{
  m_step = step;
}

void
yl::PropagateAlongField::setMaxIter(unsigned int max_iter)
{
  m_max_iter = max_iter;
}


#include "propagate_along_field.tcc"

template
int16_t yl::PropagateAlongField::ascend_until_nonzero<int16_t>(
    const Point3df &start_point,
    const carto::VolumeRef<int16_t> &seeds,
    int16_t ignore_label
) const;
template
carto::VolumeRef<int16_t>
yl::PropagateAlongField::propagate_regions<int16_t>(
    const carto::VolumeRef<int16_t> &seeds,
    int16_t target_label=0
) const;

template
int32_t yl::PropagateAlongField::ascend_until_nonzero<int32_t>(
    const Point3df &start_point,
    const carto::VolumeRef<int32_t> &seeds,
    int32_t ignore_label
) const;
template
carto::VolumeRef<int32_t>
yl::PropagateAlongField::propagate_regions<int32_t>(
    const carto::VolumeRef<int32_t> &seeds,
    int32_t target_label=0
) const;
