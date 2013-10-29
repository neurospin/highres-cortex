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

  LabelVolume(const carto::VolumeRef<Tlabel> &vol, Tlabel background = 0);
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
