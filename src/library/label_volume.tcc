namespace yl
{

template <typename Tlabel>
LabelVolume<Tlabel>::
LabelVolume(const carto::VolumeRef<Tlabel> &vol,
            Tlabel background)
  : m_background(background), m_volume(vol), m_bucketmap()
{
  const int size_x = m_volume.getSizeX();
  const int size_y = m_volume.getSizeY();
  const int size_z = m_volume.getSizeZ();

  for(int z = 0; z < size_z; ++z)
  for(int y = 0; y < size_y; ++y)
  for(int x = 0; x < size_x; ++x)
  {
    const Tlabel label = m_volume.at(x, y, z);
    if(label != background) {
      m_bucketmap[label].insert(Bucket::value_type(Point3d(x, y, z), Void()));
    }
  }
}

template <typename Tlabel>
void LabelVolume<Tlabel>::
merge_regions(Tlabel eating_label, Tlabel eaten_label)
{
  const BucketMap::iterator eating_bucket_it = m_bucketmap.find(eating_label);
  BucketMap::iterator eaten_bucket_it = m_bucketmap.find(eaten_label);

  assert(eating_bucket_it != m_bucketmap.end());
  assert(eaten_bucket_it != m_bucketmap.end());

  Bucket & eating_bucket = eating_bucket_it->second;
  Bucket & eaten_bucket = eaten_bucket_it->second;

  for(Bucket::iterator
        point_it = eaten_bucket.begin(),
        point_it_end = eaten_bucket.end();
      point_it != point_it_end;
      ++point_it) {
    const Point3d & point = point_it->first;
    const int & x = point[0];
    const int & y = point[1];
    const int & z = point[2];
    m_volume.at(x, y, z) = eating_label;
    eating_bucket.insert(*point_it);
  }

  m_bucketmap.erase(eaten_bucket_it);
}

} // namespace yl
