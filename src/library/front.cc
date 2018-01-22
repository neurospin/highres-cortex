/*
Copyright Télécom ParisTech (2015).

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

#include "front.hh"


yl::SortedFront::
SortedFront(carto::VolumeRef<int16_t>& domain,
            const int16_t front_label,
            const int16_t done_label)
    : m_domain(domain),
      m_front_label(front_label),
      m_done_label(done_label)
{
}

Point3d
yl::SortedFront::
pop()
{
  const Point3d point = get();
  // TODO decide where to put this
  m_domain(point[0], point[1], point[2]) = m_done_label;
  m_queue.pop();
  return point;
}

void
yl::SortedFront::
add(const Point3d& point, const float priority)
{
  assert(m_domain(point[0], point[1], point[2]) != m_front_label);

  m_queue.push(std::make_pair(priority, point));
  m_domain(point[0], point[1], point[2]) = m_front_label;
}
