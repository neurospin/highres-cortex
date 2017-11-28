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

#include <gsl/gsl_eigen.h>

#include <aims/data/fastAllocationData.h>
#include <aims/math/eigen.h>

#include "cortex.hh"

using aims::AimsFastAllocationData;

namespace
{

inline float square(float x)
{
  return x * x;
}

inline bool proj_is_valid(float xproj, float /*unused*/, float /*unused*/)
{
  // Invalid values are (-1, -1, -1), and no valid coordinate is supposed to be
  // negative
  return xproj >= 0;
}

float large_eigenvalue(const yl::MomentAccumulator& mom,
                       gsl_eigen_symm_workspace* w)
{
  if(mom.m_000 == 0)
    return 0.f;

  double centred_moment_matrix[3 * 3];

  // Only the lower triangular part is referenced
  centred_moment_matrix[0 * 3 + 0]
    = (mom.m_200 - mom.m_100 * mom.m_100 / mom.m_000) / mom.m_000;
  centred_moment_matrix[1 * 3 + 0]
    = (mom.m_110 - mom.m_100 * mom.m_010 / mom.m_000) / mom.m_000;
  centred_moment_matrix[1 * 3 + 1]
    = (mom.m_020 - mom.m_010 * mom.m_010 / mom.m_000) / mom.m_000;
  centred_moment_matrix[2 * 3 + 0]
    = (mom.m_101 - mom.m_100 * mom.m_001 / mom.m_000) / mom.m_000;
  centred_moment_matrix[2 * 3 + 1]
    = (mom.m_011 - mom.m_010 * mom.m_001 / mom.m_000) / mom.m_000;
  centred_moment_matrix[2 * 3 + 2]
    = (mom.m_002 - mom.m_001 * mom.m_001 / mom.m_000) / mom.m_000;

  double eigenvalues[3];

  gsl_matrix_view mat_view
    = gsl_matrix_view_array(centred_moment_matrix, 3, 3);
  gsl_vector_view eval_view
    = gsl_vector_view_array(eigenvalues, 3);

  gsl_eigen_symm(&mat_view.matrix, &eval_view.vector, w);

  return std::max(0.f, static_cast<float>(gsl_vector_max(&eval_view.vector)));
};

float large_eigenvalue(const yl::MomentAccumulator& moments)
{
  struct WorkspaceHolder
  {
    WorkspaceHolder(gsl_eigen_symm_workspace* w_) : w(w_) {};
    ~WorkspaceHolder() {gsl_eigen_symm_free(w);};
    gsl_eigen_symm_workspace* const w;
  };
  WorkspaceHolder wh(gsl_eigen_symm_alloc(3));
  return large_eigenvalue(moments, wh.w);
};


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

  ChainedIterator(const BaseIterator& begin)
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
  region_size = 0;
  touches_CSF = false;
  touches_white = false;

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
