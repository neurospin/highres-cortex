#ifndef CORTEX_COLUMN_REGION_QUALITY_HH_INCLUDED
#define CORTEX_COLUMN_REGION_QUALITY_HH_INCLUDED

#include <cartodata/volume/volume.h>
#include <yleprince/label_volume.hh>

namespace yl
{

class CortexColumnRegionQuality
{
public:
  CortexColumnRegionQuality(const carto::VolumeRef<float>& distmap_from_CSF,
                            const carto::VolumeRef<float>& distmap_from_white);

  template <typename Tlabel>
  float evaluate(const LabelVolume<Tlabel>&, Tlabel) const;

  template <typename Tlabel>
  float evaluate(const LabelVolume<Tlabel>&, Tlabel, Tlabel) const;
private:
  carto::VolumeRef<float> m_distmap_from_CSF;
  carto::VolumeRef<float> m_distmap_from_white;
}; // class CortexColumnRegionQuality


template <typename Tlabel>
inline float CortexColumnRegionQuality::
evaluate(const LabelVolume<Tlabel>& label_vol, const Tlabel label) const
{
  return evaluate(label_vol, label, label);
}


}; // namespace yl

#endif // !defined(CORTEX_COLUMN_REGION_QUALITY_HH_INCLUDED)
