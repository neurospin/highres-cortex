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

#include <highres-cortex/cortex_column_region_quality.hh>


#include <highres-cortex/cortex_column_region_quality.tcc>

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

float diameter_to_pseudo_area(float diameter)
{
  return 0.125f * square(diameter);
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
  setShapeParametres(default_goal_diameter());
  assert(m_CSF_projections.getSizeT() == 3);
  assert(m_white_projections.getSizeT() == 3);
}

void
yl::CortexColumnRegionQuality::
setShapeParametres(float goal_diameter)
{
  m_pseudo_area_cutoff = diameter_to_pseudo_area(goal_diameter);
}

float yl::CortexColumnRegionQuality::default_goal_diameter()
{
  return 0.5f;
}
