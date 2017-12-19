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

#include "cortex_column_region_quality.hh"

#include <cmath>
#include <limits>
#include <algorithm>

#include <aims/data/fastAllocationData.h>
#include <aims/math/eigen.h>

#include "cortex.hh"

namespace
{

template <typename T>
inline T square(T x)
{
  return x * x;
}

// ----------------------------------------------------------------------------
// Numerical diagonalization of 3x3 matrcies
// Copyright (C) 2006  Joachim Kopp
// ----------------------------------------------------------------------------
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
// ----------------------------------------------------------------------------
void dsyevc3(const double A[3][3], double w[3])
// ----------------------------------------------------------------------------
// Calculates the eigenvalues of a symmetric 3x3 matrix A using Cardano's
// analytical algorithm.
// Only the diagonal and upper triangular parts of A are accessed. The access
// is read-only.
// ----------------------------------------------------------------------------
// Parameters:
//   A: The symmetric input matrix
//   w: Storage buffer for eigenvalues
// ----------------------------------------------------------------------------
{
  static const double SQRT3 = 1.73205080756887729352744634151;
  double m, c1, c0;

  // Determine coefficients of characteristic poynomial. We write
  //       | a   d   f  |
  //  A =  | d*  b   e  |
  //       | f*  e*  c  |
  double de = A[0][1] * A[1][2]; // d * e
  double dd = square(A[0][1]); // d^2
  double ee = square(A[1][2]); // e^2
  double ff = square(A[0][2]); // f^2
  m  = A[0][0] + A[1][1] + A[2][2];
  c1 = (A[0][0] * A[1][1] + A[0][0] * A[2][2] + A[1][1] * A[2][2])
          - (dd + ee + ff); // a*b + a*c + b*c - d^2 - e^2 - f^2
  c0 = A[2][2] * dd + A[0][0] * ee + A[1][1] * ff - A[0][0] * A[1][1] * A[2][2]
            - 2.0 * A[0][2]*de; // c*d^2 + a*e^2 + b*f^2 - a*b*c - 2*f*d*e)

  double p, sqrt_p, q, c, s, phi;
  p = square(m) - 3.0 * c1;
  q = m * (p - (3.0 / 2.0) * c1) - (27.0 / 2.0) * c0;
  sqrt_p = std::sqrt(std::fabs(p));

  phi = 27.0 * (0.25 * square(c1) * (p - c1) + c0 * (q + 27.0 / 4.0 * c0));
  phi = (1.0 / 3.0) * std::atan2(std::sqrt(std::fabs(phi)), q);

  c = sqrt_p * std::cos(phi);
  s = (1.0 / SQRT3) * sqrt_p * std::sin(phi);

  w[1]  = (1.0 / 3.0) * (m - c);
  w[2]  = w[1] + s;
  w[0]  = w[1] + c;
  w[1] -= s;
}

inline bool proj_is_valid(float xproj, float /*unused*/, float /*unused*/)
{
  // Invalid values are (-1, -1, -1), and no valid coordinate is supposed to be
  // negative
  return xproj >= 0;
}

float large_eigenvalue(const yl::MomentAccumulator& mom)
{
  if(mom.m_000 == 0)
    return 0.f;

  double centred_moment_matrix[3][3];

  // Only the lower triangular part is referenced
  centred_moment_matrix[0][0]
    = (mom.m_200 - mom.m_100 * mom.m_100 / mom.m_000) / mom.m_000;
  centred_moment_matrix[0][1]
    = (mom.m_110 - mom.m_100 * mom.m_010 / mom.m_000) / mom.m_000;
  centred_moment_matrix[1][1]
    = (mom.m_020 - mom.m_010 * mom.m_010 / mom.m_000) / mom.m_000;
  centred_moment_matrix[0][2]
    = (mom.m_101 - mom.m_100 * mom.m_001 / mom.m_000) / mom.m_000;
  centred_moment_matrix[1][2]
    = (mom.m_011 - mom.m_010 * mom.m_001 / mom.m_000) / mom.m_000;
  centred_moment_matrix[2][2]
    = (mom.m_002 - mom.m_001 * mom.m_001 / mom.m_000) / mom.m_000;

  double eigenvalues[3];
  // Compute the eigenvalues of the centred moment matrix
  dsyevc3(centred_moment_matrix, eigenvalues);

  return static_cast<float>(
    std::max(0., *std::max_element(eigenvalues, eigenvalues + 3)));
}

template <class BaseIterator>
class ChainedIterator
  : public boost::iterator_adaptor<
  ChainedIterator<BaseIterator>,
  BaseIterator,
  boost::use_default,
  boost::forward_traversal_tag>
{
public:
  ChainedIterator()
    : ChainedIterator::iterator_adaptor_(0) {}

  ChainedIterator(const BaseIterator& begin_first,
                  const BaseIterator& end_first,
                  const BaseIterator& begin_second)
    : ChainedIterator::iterator_adaptor_(begin_first),
      m_end_first(end_first),
      m_begin_second(begin_second) {}

  explicit ChainedIterator(const BaseIterator& begin)
    : ChainedIterator::iterator_adaptor_(begin),
      m_end_first(),
      m_begin_second() {}


private:
  const BaseIterator m_end_first, m_begin_second;

  friend class boost::iterator_core_access;
  void increment()
  {
    BaseIterator& base = this->base_reference();
    ++base;
    if(base == m_end_first && m_end_first != BaseIterator())
      base = m_begin_second;
  }
};

} // end of anonymous namespace


namespace yl
{

float CortexColumnRegionQuality::fusion_ordering(const Cache& cache) const
{
  const std::size_t region_size = cache.region_size();

  return region_size;
}

bool CortexColumnRegionQuality::want_fusion(const Cache& cache) const
{
  return pseudo_area(cache) < m_pseudo_area_cutoff;
}

float CortexColumnRegionQuality::
pseudo_area(const Cache& cache) const
{
  const yl::MomentAccumulator& CSF_moment = cache.CSF_moments();
  const yl::MomentAccumulator& white_moment = cache.white_moments();

  // The moments are calculated based on the projected point cloud, thus the
  // thicker regions lead to a denser point cloud and weigh more. Is this good?

  const float CSF_large_eigenval = large_eigenvalue(CSF_moment);
  const float white_large_eigenval = large_eigenvalue(white_moment);

  return CSF_large_eigenval + white_large_eigenval;
}

template <class PointIterator>
CortexColumnRegionQuality::Cache CortexColumnRegionQuality::
cache(const PointIterator& point_it_begin,
      const PointIterator& point_it_end) const
{
  Cache cache;
  yl::MomentAccumulator& CSF_moment = cache.CSF_moments();
  yl::MomentAccumulator& white_moment = cache.white_moments();
  std::size_t& region_size = cache.region_size();
  bool& touches_CSF = cache.touches_CSF();
  bool& touches_white = cache.touches_white();

  for(PointIterator point_it = point_it_begin;
      point_it != point_it_end;
      ++point_it)
  {
    const Point3d& point = *point_it;
    const int x = point[0];
    const int y = point[1];
    const int z = point[2];

    ++region_size;

    {
      const float xproj = m_CSF_projections.at(x, y, z, 0);
      const float yproj = m_CSF_projections.at(x, y, z, 1);
      const float zproj = m_CSF_projections.at(x, y, z, 2);
      if(proj_is_valid(xproj, yproj, zproj)) {
        CSF_moment.update(xproj, yproj, zproj);
      }
    }

    {
      const float xproj = m_white_projections.at(x, y, z, 0);
      const float yproj = m_white_projections.at(x, y, z, 1);
      const float zproj = m_white_projections.at(x, y, z, 2);
      if(proj_is_valid(xproj, yproj, zproj)) {
        white_moment.update(xproj, yproj, zproj);
      }
    }

    const int16_t voxel_classif = m_classif.at(x, y, z);
    if(voxel_classif == CSF_LABEL)
      touches_CSF = true;
    else if(voxel_classif == WHITE_LABEL)
      touches_white = true;
    else if(voxel_classif != CORTEX_LABEL)
      std::clog << "WARNING: CortexColumnRegionQuality: unknown classif label "
                << voxel_classif << std::endl;
  }

  return cache;
}

template <typename Tlabel>
CortexColumnRegionQuality::Cache CortexColumnRegionQuality::
cache(const LabelVolume<Tlabel>& label_vol, const Tlabel label) const
{
  return cache(label_vol.region_begin(label),
               label_vol.region_end(label));
}

template <class PointIterator>
float CortexColumnRegionQuality::
fusion_ordering(const PointIterator& point_it_begin,
                 const PointIterator& point_it_end) const
{
  return fusion_ordering(cache(point_it_begin, point_it_end));
}

template <typename Tlabel>
float CortexColumnRegionQuality::
fusion_ordering(const LabelVolume<Tlabel>& label_vol, const Tlabel label) const
{
  return fusion_ordering(label_vol.region_begin(label),
                          label_vol.region_end(label));
}

template <typename Tlabel>
float CortexColumnRegionQuality::
fusion_ordering(const LabelVolume<Tlabel>& label_vol,
                 const Tlabel label1,
                 const Tlabel label2) const
{
  typedef ChainedIterator<typename LabelVolume<Tlabel>::const_point_iterator>
    ChainedPointIterator;
  return fusion_ordering(
    ChainedPointIterator(label_vol.region_begin(label1),
                         label_vol.region_end(label1),
                         label_vol.region_begin(label2)),
    ChainedPointIterator(label_vol.region_end(label2)));
}

}; // namespace yl
