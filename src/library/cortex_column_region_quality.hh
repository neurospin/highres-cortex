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

// TODO rename this to cortical_traverse_quality

#ifndef CORTEX_COLUMN_REGION_QUALITY_HH_INCLUDED
#define CORTEX_COLUMN_REGION_QUALITY_HH_INCLUDED

#include <cartodata/volume/volume.h>
#include "label_volume.hh"

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
                            const carto::VolumeRef<float>& white_projections,
                            const carto::VolumeRef<int16_t>& classif);

  void setShapeParametres(float goal_diameter);

  template <typename Tlabel>
  float fusion_ordering(const LabelVolume<Tlabel>&, Tlabel) const;

  template <typename Tlabel>
  float fusion_ordering(const LabelVolume<Tlabel>&, Tlabel, Tlabel) const;

  template <class PointIterator>
  inline float fusion_ordering(const PointIterator& point_it_begin,
                               const PointIterator& point_it_end) const;

  class Cache
  {
  public:
    Cache() : m_CSF_moments(), m_white_moments() {};
    Cache(const MomentAccumulator& CSF_moments,
          const MomentAccumulator& white_moments,
          std::size_t region_size,
          bool touches_CSF,
          bool touches_white)
      : m_CSF_moments(CSF_moments),
        m_white_moments(white_moments),
        m_region_size(region_size),
        m_touches_CSF(touches_CSF),
        m_touches_white(touches_white)
    {};

    Cache operator + (const Cache& other) const
    {
      return Cache(CSF_moments() + other.CSF_moments(),
                   white_moments() + other.white_moments(),
                   region_size() + other.region_size(),
                   touches_CSF() || other.touches_CSF(),
                   touches_white() || other.touches_white());
    };

    const MomentAccumulator& CSF_moments() const {return m_CSF_moments;};
    MomentAccumulator& CSF_moments() {return m_CSF_moments;};
    const MomentAccumulator& white_moments() const {return m_white_moments;};
    MomentAccumulator& white_moments() {return m_white_moments;};
    std::size_t region_size() const {return m_region_size;};
    std::size_t& region_size() {return m_region_size;};
    bool touches_CSF() const {return m_touches_CSF;};
    bool& touches_CSF() {return m_touches_CSF;};
    bool touches_white() const {return m_touches_white;};
    bool& touches_white() {return m_touches_white;};
  private:
    MomentAccumulator m_CSF_moments;
    MomentAccumulator m_white_moments;
    std::size_t m_region_size;
    bool m_touches_CSF;
    bool m_touches_white;
  };
  float fusion_ordering(const Cache&) const;
  bool want_fusion (const Cache&) const;
  float pseudo_area(const Cache&) const;
  template <class PointIterator>
  Cache cache(const PointIterator& point_it_begin,
              const PointIterator& point_it_end) const;
  template <typename Tlabel>
  Cache cache(const LabelVolume<Tlabel>&, Tlabel) const;

  static float default_goal_diameter();

private:
  carto::VolumeRef<float> m_CSF_projections;
  carto::VolumeRef<float> m_white_projections;
  carto::VolumeRef<int16_t> m_classif;
  float m_pseudo_area_cutoff;
}; // class CortexColumnRegionQuality

}; // namespace yl

#endif // !defined(CORTEX_COLUMN_REGION_QUALITY_HH_INCLUDED)
