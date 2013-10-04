#ifndef YL_LABEL_VOLUME_HH_INCLUDED
#define YL_LABEL_VOLUME_HH_INCLUDED

#include <functional>
#include <stdexcept>
#include <boost/iterator/transform_iterator.hpp>

#include <cartodata/volume/volume.h>
#include <aims/bucket/bucketMap.h>

namespace yl
{

/** Return the first element of a std::pair, like non-standard std::select1st */
template<typename TPair>
class select1st : public std::unary_function<TPair, typename TPair::first_type>
{
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

  /** Iterator through a Bucket's points (yields Point3d) */
  typedef boost::transform_iterator<select1st<Bucket::value_type>,
                                    Bucket::const_iterator> const_point_iterator;

  LabelVolume(const carto::Volume<Tlabel> &vol, Tlabel background = 0);
  virtual ~LabelVolume() {};

  void merge_regions(Tlabel eating_label, Tlabel eaten_label);
  void change_voxel_label(const Point3d &point, Tlabel new_label);

  /** Iterate through the coordinates of a region's voxels */
  const_point_iterator region_begin(Tlabel label) const
  {
    const BucketMap::const_iterator region_bucket_it = m_bucketmap.find(label);
    if(region_bucket_it == m_bucketmap.end())
      throw std::runtime_error("yl::LabelVolume::region_begin: "
                               "no such region");
    else
      return const_point_iterator(region_bucket_it->second.begin(),
                                  select1st<Bucket::value_type>());
  };

  /** Iterate through the coordinates of a region's voxels */
  const_point_iterator region_end(Tlabel label) const
  {
    const BucketMap::const_iterator region_bucket_it = m_bucketmap.find(label);
    if(region_bucket_it == m_bucketmap.end())
      throw std::runtime_error("yl::LabelVolume::region_end: "
                               "no such region");
    else
      return const_point_iterator(region_bucket_it->second.end(),
                                  select1st<Bucket::value_type>());
  };

private:
  Tlabel m_background;
  carto::Volume<Tlabel> m_volume;
  aims::BucketMap<Void> m_bucketmap;
}; // class LabelVolume

}; // namespace yl

#endif // !defined(YL_LABEL_VOLUME_HH_INCLUDED)
