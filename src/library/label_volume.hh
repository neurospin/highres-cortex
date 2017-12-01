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

#ifndef YL_LABEL_VOLUME_HH_INCLUDED
#define YL_LABEL_VOLUME_HH_INCLUDED

#include <cassert>
#include <functional>
#include <boost/iterator/transform_iterator.hpp>

#include <cartodata/volume/volume.h>
#include <aims/bucket/bucketMap.h>

namespace yl
{

/** Return the first element of a std::pair, like non-standard std::select1st */
template<typename TPair>
class select1st : public std::unary_function<TPair, typename TPair::first_type>
{
public:
  typename TPair::first_type&
  operator()(TPair& pair) const
  {
    return pair.first;
  }

  const typename TPair::first_type&
  operator()(const TPair& pair) const
  {
    return pair.first;
  }
};


template <typename Tlabel>
class LabelVolume
{
public:
  typedef aims::BucketMap<Void> BucketMap;
  typedef BucketMap::Bucket Bucket;

  typedef boost::transform_iterator<select1st<BucketMap::value_type>,
                                    BucketMap::const_iterator> const_regions_iterator;
  /** Iterator through a Bucket's points (yields Point3d) */
  typedef boost::transform_iterator<select1st<Bucket::value_type>,
                                    Bucket::const_iterator> const_point_iterator;

  explicit LabelVolume(const carto::VolumeRef<Tlabel> &vol, Tlabel background = 0);
  virtual ~LabelVolume() {};

  void merge_regions(Tlabel eating_label, Tlabel eaten_label);
  void discard_region(Tlabel label);

  Tlabel background_label() const
  {
    return m_background_label;
  }

  const carto::VolumeRef<Tlabel>& volume() const
  {
    return m_volume;
  }

  BucketMap::size_type n_regions() const
  {
    return m_bucketmap.size();
  }

  /** Iterate through all region labels */
  const_regions_iterator regions_begin() const
  {
    return const_regions_iterator(m_bucketmap.begin(),
                                  select1st<BucketMap::value_type>());
  };

  /** Iterate through all region labels */
  const_regions_iterator regions_end() const
  {
    return const_regions_iterator(m_bucketmap.end(),
                                  select1st<BucketMap::value_type>());
  };

  BucketMap::size_type region_size(Tlabel label) const
  {
    const BucketMap::const_iterator region_bucket_it = m_bucketmap.find(label);
    assert(region_bucket_it != m_bucketmap.end());
    return region_bucket_it->second.size();
  };

  /** Iterate through the coordinates of a region's voxels */
  const_point_iterator region_begin(Tlabel label) const
  {
    const BucketMap::const_iterator region_bucket_it = m_bucketmap.find(label);
    assert(region_bucket_it != m_bucketmap.end());
    return const_point_iterator(region_bucket_it->second.begin(),
                                select1st<Bucket::value_type>());
  };

  /** Iterate through the coordinates of a region's voxels */
  const_point_iterator region_end(Tlabel label) const
  {
    const BucketMap::const_iterator region_bucket_it = m_bucketmap.find(label);
    assert(region_bucket_it != m_bucketmap.end());
    return const_point_iterator(region_bucket_it->second.end(),
                                select1st<Bucket::value_type>());
  };

private:
  const Tlabel m_background_label;
  carto::VolumeRef<Tlabel> m_volume;
  aims::BucketMap<Void> m_bucketmap;
}; // class LabelVolume

}; // namespace yl

#endif // !defined(YL_LABEL_VOLUME_HH_INCLUDED)
