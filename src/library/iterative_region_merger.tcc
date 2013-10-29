#include <cassert>
#include <algorithm>
#include <functional>
#include <iostream>
#include <iterator>
#include <set>

#include <boost/heap/d_ary_heap.hpp>
#include <boost/iterator/indirect_iterator.hpp>

namespace
{

const int debug_output = 0;

inline bool approx_equal(float x, float y)
{
  return std::abs(x - y) < 0.1f;
}

template <typename Tlabel> class Region;
template <typename Tlabel>
std::ostream& operator << (std::ostream& stream, const Region<Tlabel>& region);

template <typename Tlabel>
class Region
{
  friend std::ostream& operator << <Tlabel> (std::ostream&, const Region<Tlabel>&);

public:
  typedef boost::indirect_iterator<typename std::set<Region*>::const_iterator>
  const_neighbour_iterator;
  typedef boost::indirect_iterator<typename std::set<Region*>::iterator>
  neighbour_iterator;

  Region(const Tlabel region_label, const float initial_quality)
    : m_label(region_label), m_neighbours(), m_quality(initial_quality)
  {
  }

  Region(const Region& other)
    : m_label(other.m_label),
      m_neighbours(other.m_neighbours),
      m_quality(other.m_quality)
  {
    // Ensure reciprocity of the neighbourhood relationship
    std::for_each(m_neighbours.begin(), m_neighbours.end(),
                  std::bind2nd(std::mem_fun(&Region::add_forward_neighbour), this));
  }

  virtual ~Region()
  {
    // Remove this disappearing node from the neighbourhood of its neighbours
    // (prevent broken edges in the region neighbourhood graph)
    std::for_each(m_neighbours.begin(), m_neighbours.end(),
                  std::bind2nd(std::mem_fun(&Region::remove_forward_neighbour), this));
  }

  Tlabel label() const
  {
    return m_label;
  }

  void add_neighbour(Region<Tlabel>& neighbour_region)
  {
    assert(&neighbour_region != this);
    m_neighbours.insert(&neighbour_region);
    neighbour_region.m_neighbours.insert(this);
  }

  float quality() const
  {
    return m_quality;
  }

  void update_quality(const float new_quality)
  {
    m_quality = new_quality;
  }

  const_neighbour_iterator neighbours_begin() const
  { return m_neighbours.begin(); };
  const_neighbour_iterator neighbours_end() const
  { return m_neighbours.end(); };
  neighbour_iterator neighbours_begin()
  { return m_neighbours.begin(); };
  neighbour_iterator neighbours_end()
  { return m_neighbours.end(); };

  bool operator == (const Region& other) const
  {
    return this == &other;
  }

  bool operator != (const Region& other) const
  {
    return this != &other;
  }

private:
  void add_forward_neighbour(Region<Tlabel>* neighbour_region)
  {
    m_neighbours.insert(neighbour_region);
  }

  void remove_forward_neighbour(Region<Tlabel>* neighbour_region)
  {
    m_neighbours.erase(neighbour_region);
  }

  const Tlabel m_label;
  std::set<Region*> m_neighbours;
  float m_quality;
};

template <typename Tlabel>
std::ostream& operator << (std::ostream& stream, const Region<Tlabel>& region)
{
  stream << "Region@" << &region << " (" << "label=" << region.label()
         << ", quality=" << region.quality()  << ", #neighbours="
         << region.m_neighbours.size();
  return stream;
}

template <typename Tlabel>
class RegionInQueue : public Region<Tlabel>
{
public:
  typedef boost::heap::d_ary_heap<RegionInQueue<Tlabel>,
                                  boost::heap::arity<2>,
                                  boost::heap::mutable_<true> > RegionQueue;
  typedef typename RegionQueue::handle_type Handle;

  RegionInQueue(const Tlabel region_label, const float initial_quality)
    : Region<Tlabel>(region_label, initial_quality), m_handle()
  {
  }

  RegionInQueue(const RegionInQueue<Tlabel>& other)
    : Region<Tlabel>(other), m_handle(other.m_handle)
  {
  }

  bool operator < (const Region<Tlabel>& other) const
  {
    return this->quality() > other.quality();
  }

  void set_handle(const Handle& handle)
  {
    m_handle = handle;
  }

  const Handle& handle()
  {
    return m_handle;
  }

private:
  Handle m_handle;
};

}; // end of anonymous namespace

namespace yl
{

template <typename Tlabel, class RegionQualityCriterion>
IterativeRegionMerger<Tlabel, RegionQualityCriterion>::
IterativeRegionMerger(const LabelVolume<Tlabel>& label_vol,
                      const RegionQualityCriterion& criterion,
                      const int verbosity)
  : m_label_volume(label_vol), m_criterion(criterion),
    m_max_region_size(std::numeric_limits<std::size_t>::max()),
    m_verbosity(verbosity)
{
}

template <typename Tlabel, class RegionQualityCriterion>
void
IterativeRegionMerger<Tlabel, RegionQualityCriterion>::
setVerbose(const int verbosity)
{
  m_verbosity = verbosity;
}

template <typename Tlabel, class RegionQualityCriterion>
void
IterativeRegionMerger<Tlabel, RegionQualityCriterion>::
merge_worst_regions_iteratively()
{
  if(m_verbosity) {
    std::clog << "yl::IterativeRegionMerger::merge_worst_regions_iteratively:\n"
              << "  computing initial region qualities..."
              << std::endl;
  }

  typedef boost::heap::d_ary_heap<RegionInQueue<Tlabel>,
                                  boost::heap::arity<2>,
                                  boost::heap::mutable_<true> > RegionQueue;
  typedef typename RegionQueue::handle_type Handle;

  // Hold the labels to retrieve them in order of increasing region quality
  RegionQueue queue;
  std::map<Tlabel, Handle> label_to_handle;

  unsigned int oversized_regions_ignored = 0;
  for(typename LabelVolume<Tlabel>::const_regions_iterator
        labels_it = m_label_volume.regions_begin(),
        labels_end = m_label_volume.regions_end();
      labels_it != labels_end;
      ++labels_it)
  {
    const Tlabel label = *labels_it;

    if(m_label_volume.region_size(label) > m_max_region_size) {
      ++oversized_regions_ignored;
      continue;
    }

    Handle handle = queue.push(
      RegionInQueue<Tlabel>(label,
                            m_criterion.evaluate(m_label_volume, label)));
    (*handle).set_handle(handle);
    assert(label_to_handle.find(label) == label_to_handle.end());
    label_to_handle[label] = handle;
  }

  if(m_verbosity) {
    std::clog << "  ignoring " << oversized_regions_ignored
              << " regions larger than " << m_max_region_size << " voxels.\n"
              << "  " << queue.size() << " regions will be processed.\n"
              << "  filling in neighbourhoods..." << std::endl;
  }

  {
    const Tlabel background_label = m_label_volume.background_label();
    const int size_x = m_label_volume.volume().getSizeX();
    const int size_y = m_label_volume.volume().getSizeY();
    const int size_z = m_label_volume.volume().getSizeZ();
    for(int z = 0; z < size_z - 1; ++z)
    for(int y = 0; y < size_y - 1; ++y)
    for(int x = 0; x < size_x - 1; ++x)
    {
      const Tlabel self_label = m_label_volume.volume().at(x, y, z);

      // If the region was not included into the queue
      if(label_to_handle.count(self_label) == 0)
        continue;

      if(self_label != background_label) {
        Region<Tlabel>& self_region = *label_to_handle[self_label];
        const Tlabel xplus_label = m_label_volume.volume().at(x + 1, y, z);
        const Tlabel yplus_label = m_label_volume.volume().at(x, y + 1, z);
        const Tlabel zplus_label = m_label_volume.volume().at(x, y, z + 1);

        if(label_to_handle.count(xplus_label) != 0
           && xplus_label != background_label && self_label != xplus_label) {
          Region<Tlabel>& xplus_region = *label_to_handle[xplus_label];
          self_region.add_neighbour(xplus_region);
        }

        if(label_to_handle.count(yplus_label) != 0
           && yplus_label != background_label && self_label != yplus_label) {
          Region<Tlabel>& yplus_region = *label_to_handle[yplus_label];
          self_region.add_neighbour(yplus_region);
        }

        if(label_to_handle.count(zplus_label) != 0
           && zplus_label != background_label && self_label != zplus_label) {
          Region<Tlabel>& zplus_region = *label_to_handle[zplus_label];
          self_region.add_neighbour(zplus_region);
        }
      }
    }
  }

  if(m_verbosity) {
    std::clog << "  iteratively merging regions..." << std::endl;
  }

  // Iteratively merge the worst region with one of its neighbours, until no
  // region can be improved further by merging a neighbour.
  while(!queue.empty())
  {
    const RegionInQueue<Tlabel>& worst_region = queue.top();
    const Tlabel worst_label = worst_region.label();

    // const_cast is ok because the "best neighbour region" is used only if it
    // is actually a neighbour region. The pointer is left alone if it retains
    // its original value (see "if" below).
    Region<Tlabel>* best_neighbour_region_p
      = const_cast<RegionInQueue<Tlabel>*>(&worst_region);
    float best_quality = worst_region.quality();

    if(m_verbosity >= 2) {
      std::clog << "  ";
      if(RegionQueue::constant_time_size) {
        std::clog << queue.size() << " to go, ";
      }
      std::clog << m_label_volume.n_regions() << " regions, q = "
                << best_quality << "\r" << std::flush;
    }

    // assert(approx_equal(m_criterion.evaluate(m_label_volume, worst_label),
    //                     best_quality));

    for(typename Region<Tlabel>::neighbour_iterator
          neighbour_it = worst_region.neighbours_begin(),
          neighbour_end = worst_region.neighbours_end();
        neighbour_it != neighbour_end;
        ++neighbour_it)
    {
      Region<Tlabel>& neighbour_region = *neighbour_it;
      const Tlabel neighbour_label = neighbour_region.label();
      const float conjunction_quality =
        m_criterion.evaluate(m_label_volume, worst_label, neighbour_label);
      if(conjunction_quality > best_quality) {
        best_quality = conjunction_quality;
        best_neighbour_region_p = &neighbour_region;
      }
    }

    RegionInQueue<Tlabel>& best_neighbour_region
      = *dynamic_cast<RegionInQueue<Tlabel>*>(best_neighbour_region_p);
    const Tlabel best_neighbour_label = best_neighbour_region.label();

    if(best_neighbour_region != worst_region) {
      if(m_verbosity >= 3) {
        std::clog << "\n    merging with best neighbour " << best_neighbour_region
                  << " (new q=" << best_quality << ")" << std::endl;
      }

      // Region with worst_label is eaten by its neighbour_label
      m_label_volume.merge_regions(best_neighbour_label, worst_label);

      // Update new region's neighbourhood
      for(typename Region<Tlabel>::neighbour_iterator
            neighbour_it = worst_region.neighbours_begin(),
            neighbour_end = worst_region.neighbours_end();
          neighbour_it != neighbour_end;
          ++neighbour_it)
      {
        Region<Tlabel>& neighbour_region = *neighbour_it;
        if(neighbour_region != best_neighbour_region) {
          best_neighbour_region.add_neighbour(neighbour_region);
        }
      }

      // Get rid of worst_region. This needs to be done *before* updating the
      // heap!
      queue.pop();

      best_neighbour_region.update_quality(best_quality);
      queue.update(best_neighbour_region.handle());
    } else {
      if(m_verbosity >= 3) {
        std::clog << "\n    region " << worst_label << " (quality=" << best_quality
                  << ") cannot be improved by merging a neighbour" << std::endl;
      }
      queue.pop();
    }
  }
  if(m_verbosity >= 1) {
    std::clog << "end: " << m_label_volume.n_regions() << " regions." << std::endl;
  }
}

};
