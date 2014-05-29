#ifndef CORTEX_COLUMN_REGION_QUALITY_HH_INCLUDED
#define CORTEX_COLUMN_REGION_QUALITY_HH_INCLUDED

#include <cartodata/volume/volume.h>
#include <yleprince/label_volume.hh>

namespace yl
{

struct MomentAccumulator
{
public:
  MomentAccumulator()
    : m_000(0.f),
      m_100(0.f), m_010(0.f), m_001(0.f),
      m_200(0.f),
      m_110(0.f), m_020(0.f),
      m_101(0.f), m_011(0.f), m_002(0.f) {}

  MomentAccumulator(float m000,
                    float m100, float m010, float m001,
                    float m200,
                    float m110, float m020,
                    float m101, float m011, float m002)
    : m_000(m000),
      m_100(m100), m_010(m010), m_001(m001),
      m_200(m200),
      m_110(m110), m_020(m020),
      m_101(m101), m_011(m011), m_002(m002) {}

  void update(float x, float y, float z)
  {
    m_000 += 1;
    m_100 += x; m_010 += y; m_001 += z;
    m_200 += x * x;
    m_110 += y * x; m_020 += y * y;
    m_101 += z * x; m_011 += z * y; m_002 += z * z;
  };

  MomentAccumulator operator + (const MomentAccumulator& other) const
  {
    return MomentAccumulator(m_000 + other.m_000,
                             m_100 + other.m_100,
                             m_010 + other.m_010,
                             m_001 + other.m_001,
                             m_200 + other.m_200,
                             m_110 + other.m_110,
                             m_020 + other.m_020,
                             m_101 + other.m_101,
                             m_011 + other.m_011,
                             m_002 + other.m_002);
  };

  std::size_t m_000;
  float m_100, m_010, m_001;
  float m_200;
  float m_110, m_020;
  float m_101, m_011, m_002;
};


class CortexColumnRegionQuality
{
public:
  CortexColumnRegionQuality(const carto::VolumeRef<float>& CSF_projections,
                            const carto::VolumeRef<float>& white_projections);

  void setShapeParametres(float goal_diametre, float max_thickness);

  template <typename Tlabel>
  float evaluate(const LabelVolume<Tlabel>&, Tlabel) const;

  template <typename Tlabel>
  float evaluate(const LabelVolume<Tlabel>&, Tlabel, Tlabel) const;

  template <class PointIterator>
  inline float evaluate(const PointIterator& point_it_begin,
                        const PointIterator& point_it_end) const;

  class Cache
  {
  public:
    Cache() : m_CSF_moments(), m_white_moments() {};
    Cache(const MomentAccumulator& CSF_moments,
          const MomentAccumulator& white_moments)
      : m_CSF_moments(CSF_moments), m_white_moments(white_moments) {};

    Cache operator + (const Cache& other) const
    {
      return Cache(CSF_moments() + other.CSF_moments(),
                   white_moments() + other.white_moments());
    };

    const MomentAccumulator& CSF_moments() const {return m_CSF_moments;};
    MomentAccumulator& CSF_moments() {return m_CSF_moments;};
    const MomentAccumulator& white_moments() const {return m_white_moments;};
    MomentAccumulator& white_moments() {return m_white_moments;};
    const std::size_t& region_size() const {return m_region_size;};
    std::size_t& region_size() {return m_region_size;};
  private:
    MomentAccumulator m_CSF_moments;
    MomentAccumulator m_white_moments;
    std::size_t m_region_size;
  };
  float evaluate(const Cache&) const;
  float evaluate_without_size_penalty(const Cache&) const;
  template <class PointIterator>
  Cache cache(const PointIterator& point_it_begin,
              const PointIterator& point_it_end) const;
  template <typename Tlabel>
  Cache cache(const LabelVolume<Tlabel>&, Tlabel) const;

  static float default_goal_diametre();
  static float default_max_thickness();

private:
  carto::VolumeRef<float> m_CSF_projections;
  carto::VolumeRef<float> m_white_projections;
  float m_sorted_voxel_sizes[3];
  float m_pseudo_area_reliability_threshold;
  float m_pseudo_area_cutoff;
  float m_max_criterion;
}; // class CortexColumnRegionQuality

}; // namespace yl

#endif // !defined(CORTEX_COLUMN_REGION_QUALITY_HH_INCLUDED)
