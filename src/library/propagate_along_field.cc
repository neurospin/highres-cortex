#include <yleprince/propagate_along_field.hh>

#include <cmath>
#include <stdexcept>

#include <boost/make_shared.hpp>

using carto::VolumeRef;
using std::clog;
using std::endl;
using boost::shared_ptr;
using boost::make_shared;

yl::PropagateAlongField::
PropagateAlongField(const shared_ptr<VectorField3d>& vector_field)
  : m_vector_field(vector_field),
    m_max_iter(default_max_iter), m_step(default_step), m_verbose(debug_output)
{
}

yl::PropagateAlongField::
PropagateAlongField(const VolumeRef<float>& fieldx,
                    const VolumeRef<float>& fieldy,
                    const VolumeRef<float>& fieldz)
  : m_vector_field(make_shared<LinearlyInterpolatedVectorField3d>(fieldx,
                                                                  fieldy,
                                                                  fieldz)),
    m_max_iter(default_max_iter), m_step(default_step), m_verbose(debug_output)
{
}

yl::PropagateAlongField::~PropagateAlongField()
{
}

void
yl::PropagateAlongField::setVerbose(int verbosity)
{
  m_verbose = verbosity;
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
    int16_t target_label
) const;
template
std::pair<carto::VolumeRef<int16_t>, carto::VolumeRef<float> >
yl::PropagateAlongField::propagate_regions_keeping_dests<int16_t>(
    const carto::VolumeRef<int16_t> &seeds,
    int16_t target_label
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
    int32_t target_label
) const;
template
std::pair<carto::VolumeRef<int32_t>, carto::VolumeRef<float> >
yl::PropagateAlongField::propagate_regions_keeping_dests<int32_t>(
    const carto::VolumeRef<int32_t> &seeds,
    int32_t target_label
) const;
