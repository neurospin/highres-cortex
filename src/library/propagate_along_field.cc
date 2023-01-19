/*
Copyright CEA (2014).
Copyright Université Paris XI (2014).

Contributor: Yann Leprince <yann.leprince@ylep.fr>.

This file is part of highres-cortex, a collection of software designed
to process high-resolution magnetic resonance images of the cerebral
cortex.

This software is governed by the CeCILL licence under French law and
abiding by the rules of distribution of free software. You can use,
modify and/or redistribute the software under the terms of the CeCILL
licence as circulated by CEA, CNRS and INRIA at the following URL:
<http://www.cecill.info/>.

As a counterpart to the access to the source code and rights to copy,
modify and redistribute granted by the licence, users are provided only
with a limited warranty and the software's author, the holder of the
economic rights, and the successive licensors have only limited
liability.

In this respect, the user's attention is drawn to the risks associated
with loading, using, modifying and/or developing or reproducing the
software by the user in light of its specific status of scientific
software, that may mean that it is complicated to manipulate, and that
also therefore means that it is reserved for developers and experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or
data to be ensured and, more generally, to use and operate it in the
same conditions as regards security.

The fact that you are presently reading this means that you have had
knowledge of the CeCILL licence and that you accept its terms.
*/

#include "propagate_along_field.hh"

#include <cmath>
#include <stdexcept>

#include <boost/make_shared.hpp>

using carto::VolumeRef;
using std::clog;
using std::endl;
using boost::shared_ptr;
using boost::make_shared;


const float yl::PropagateAlongField::default_step = 0.03f;
const unsigned int yl::PropagateAlongField::default_max_iter = 1000;

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
