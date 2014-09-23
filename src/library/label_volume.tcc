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

namespace yl
{

template <typename Tlabel>
LabelVolume<Tlabel>::
LabelVolume(const carto::VolumeRef<Tlabel> &vol,
            Tlabel background)
  : m_background_label(background), m_volume(vol), m_bucketmap()
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

template <typename Tlabel>
void LabelVolume<Tlabel>::
discard_region(Tlabel label)
{
  BucketMap::iterator bucket_it = m_bucketmap.find(label);
  assert(bucket_it != m_bucketmap.end());

  Bucket & bucket = bucket_it->second;

  for(Bucket::iterator
        point_it = bucket.begin(),
        point_it_end = bucket.end();
      point_it != point_it_end;
      ++point_it) {
    const Point3d & point = point_it->first;
    const int & x = point[0];
    const int & y = point[1];
    const int & z = point[2];
    m_volume.at(x, y, z) = m_background_label;
  }

  m_bucketmap.erase(bucket_it);
}

} // namespace yl
