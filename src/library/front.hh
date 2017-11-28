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

#ifndef FRONT_HH_INCLUDED
#define FRONT_HH_INCLUDED

#include <utility>

#include <boost/heap/priority_queue.hpp>
#include <boost/heap/d_ary_heap.hpp>
#include <boost/unordered_map.hpp>

#include <cartodata/volume/volume.h>
#include <aims/vector/vector.h>


namespace yl
{

/** Helper class to use Point3d & friends in hash tables */
template <class Point>
struct PointHasher : std::unary_function<Point, std::size_t>
{
  std::size_t operator()(const Point& point) const
  {
    std::size_t seed = 0;
    for(int i = 0; i < point.size(); ++i)
      boost::hash_combine(seed, point[i]);
    return seed;
  }
};

enum {
  DEFAULT_FRONT_LABEL = -123,
  DEFAULT_DONE_LABEL = -10
};

class VariablePrioritySortedFront
{
public:
  VariablePrioritySortedFront(carto::VolumeRef<int16_t>& domain,
                              int16_t front_label,
                              int16_t done_label);

  const Point3d& get() const
  {
    return m_queue.top().second;
  };
  Point3d pop();

  void add(const Point3d&, float priority);
  void update_priority(const Point3d&, float new_priority);
  void increase_priority(const Point3d&, float new_priority);
  void decrease_priority(const Point3d&, float new_priority);
  void add_or_update(const Point3d&, float priority);

  bool empty() const
  {
    return m_queue.empty();
  };

  std::size_t size() const
  {
    return m_queue.size();
  };

private:
  struct VoxelOrdering {
    bool operator() (const std::pair<float, Point3d>& left,
                     const std::pair<float, Point3d>& right) const
    {
      return left.first > right.first;
    };
  };
#if 0 // lower asymptotic complexity, but slower in our tests
  typedef boost::heap::fibonacci_heap<std::pair<float, Point3d>,
                                      boost::heap::compare<VoxelOrdering> >
    FrontQueue;
#else
  typedef boost::heap::d_ary_heap<std::pair<float, Point3d>,
                                  boost::heap::arity<8>,
                                  boost::heap::mutable_<true>,
                                  boost::heap::compare<VoxelOrdering> >
    FrontQueue;
#endif
  typedef FrontQueue::handle_type QueueHandle;
  FrontQueue m_queue;

  boost::unordered_map<Point3d, QueueHandle, PointHasher<Point3d> >
    m_handle_map;

  carto::VolumeRef<int16_t> m_domain;
  const int16_t m_front_label;
  const int16_t m_done_label;

public:
  static const bool constant_time_size = FrontQueue::constant_time_size;
};


class SortedFront
{
public:
  SortedFront(carto::VolumeRef<int16_t>& domain,
              int16_t front_label,
              int16_t done_label);

  const Point3d& get() const
  {
    return m_queue.top().second;
  };
  Point3d pop();

  void add(const Point3d&, float priority);

  bool empty() const
  {
    return m_queue.empty();
  };

  std::size_t size() const
  {
    return m_queue.size();
  };

private:
  struct VoxelOrdering {
    bool operator() (const std::pair<float, Point3d>& left,
                     const std::pair<float, Point3d>& right) const
    {
      return left.first > right.first;
    };
  };

  typedef boost::heap::priority_queue<std::pair<float, Point3d>,
                                      boost::heap::compare<VoxelOrdering> >
    FrontQueue;
  FrontQueue m_queue;

  carto::VolumeRef<int16_t> m_domain;
  const int16_t m_front_label;
  const int16_t m_done_label;

public:
  static const bool constant_time_size = FrontQueue::constant_time_size;
};

} // namespace yl

#endif // !defined(FRONT_HH_INCLUDED)
