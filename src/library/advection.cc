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

#include "advection.hh"
#include <cmath>

yl::Advection::Advection(const VectorField3d& advection_field)
  : m_vector_field(advection_field), m_verbose(0), m_max_iter(default_max_iter)
{
}

yl::ConstantStepAdvection::
ConstantStepAdvection(const yl::VectorField3d& advection_field,
                      const float step)
  : Advection(advection_field), m_step(step)
{
}

void yl::ConstantStepAdvection::
move_one_step(Point3df& point,
              const Point3df& local_field) const
{
  float gx = local_field[0], gy = local_field[1], gz = local_field[2];
  // Normalize the field, stop if too small or infinite or NaN.
  const float gn = std::sqrt(gx*gx + gy*gy + gz*gz);
  if(!std::isnormal(gn))
    throw yl::Advection::AbnormalField();
  gx /= gn; gy /= gn; gz /= gn;

  point[0] += m_step * gx;
  point[1] += m_step * gy;
  point[2] += m_step * gz;
}

#include "advection.tcc"

template bool yl::Advection::visitor_advection<yl::Advection::Visitor>
(yl::Advection::Visitor& visitor, const Point3df& start_point) const;
