#ifndef ITERATIVE_REGION_MERGER_HH_INCLUDED
#define ITERATIVE_REGION_MERGER_HH_INCLUDED

#include <cartodata/volume/volume.h>

#include <yleprince/label_volume.hh>

namespace yl
{

template <typename Tlabel, class RegionQualityCriterion>
class IterativeRegionMerger
{
public:
  IterativeRegionMerger(
    const LabelVolume<Tlabel>& label_vol,
    const RegionQualityCriterion& criterion=RegionQualityCriterion(),
    int verbosity=0);
  virtual ~IterativeRegionMerger() {};

  void setVerbose(int verbosity=1);

  void set_max_region_size(const std::size_t max_region_size)
  { m_max_region_size = max_region_size; };
  std::size_t max_region_size() const { return m_max_region_size; };

  void merge_worst_regions_iteratively();

  carto::VolumeRef<Tlabel> volume() const
  {
    return m_label_volume.volume();
  }

private:
  LabelVolume<Tlabel> m_label_volume;
  RegionQualityCriterion m_criterion;
  std::size_t m_max_region_size;
  int m_verbosity;
}; // class IterativeRegionMerger

}; // namespace yl

#endif // !defined(ITERATIVE_REGION_MERGER_HH_INCLUDED)
