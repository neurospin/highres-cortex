/*
Copyright Forschungszentrum Jülich GmbH (2017).
Copyright CEA (2014, 2017).
Copyright Université Paris XI (2014).

Contributor: Yann Leprince <yann.leprince@ylep.fr>.
Contributor: Denis Rivière <denis.riviere@cea.fr>.

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

#include "cortex_advection.hh"

#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <utility>

#include <cartodata/volume/volume.h>
#include <aims/utility/converter_volume.h>

#include "advection.hh"
#include "field.hh"

using carto::VolumeRef;
using namespace std;

// Anonymous namespace for file-local symbols
namespace {

const float no_value = std::numeric_limits<float>::quiet_NaN();

class TubeAdvection : public yl::Advection::Visitor
{
public:
  TubeAdvection(const yl::ScalarField& divergence_field,
                const yl::ScalarField& domain,
                VolumeRef<float> & volume_result,
                VolumeRef<float> & surface_result,
                const bool opposite_direction=false)
    : m_divergence_field(divergence_field),
      m_domain(domain),
      m_opposite_direction(opposite_direction),
      m_previous_point(),
      m_surface(1.f),
      m_volume(0.f),
      m_volume_result(volume_result),
      m_surface_result(surface_result),
      m_abort(false),
      m_voxel_size(volume_result->getVoxelSize())
  {
  };

  void first(const Point3df& point)
  {
    m_previous_point = point;
  };

  void visit(const Point3df& point)
  {
    const float step = (point - m_previous_point).norm();
    const float old_surface = m_surface;
    float divergence_value;
    try {
      divergence_value = m_divergence_field.evaluate(point);
    } catch(const yl::Field::UndefinedField&) {
      m_abort = true;
      return;
    }

    if(m_opposite_direction)
      m_surface *= 1 - step * divergence_value;
    else
      m_surface *= 1 + step * divergence_value;

    // The test is written inverted so that NaN gives false
    if(!(m_surface > 0)) {
      clog << "  Warning: TubeAdvection encountered non-positive surface" << endl;
      m_abort = true;
      return;
    }
    m_volume += step * (old_surface + m_surface) / 2;
    m_previous_point = point;
  };

  bool move_on(const Point3df& point) const
  {
    if(m_abort)
      return false;
    try {
      return m_domain.evaluate(point) >= 0.5f;
    } catch(const yl::Field::UndefinedField&) {
      return false;
    }
  };

  void finished(const Point3df& start_point)
  {
    Point3d sp_int = Point3d(
      lrint(start_point[0] / m_voxel_size[0]),
      lrint(start_point[1] / m_voxel_size[1]),
      lrint(start_point[2] / m_voxel_size[2]));
    m_volume_result(sp_int) = volume();
    m_surface_result(sp_int) = surface();
  }

  float surface() const
  {
    if(m_abort)
      return no_value;
    return m_surface;
  };

  float volume() const
  {
    if(m_abort)
      return no_value;
    return m_volume;
  };

  bool aborted() const
  {
    return m_abort;
  };

  VolumeRef<float>& volume_result() const
  {
    return m_volume_result;
  }

  VolumeRef<float>& surface_result() const
  {
    return m_surface_result;
  }

private:
  const yl::ScalarField& m_divergence_field;
  const yl::ScalarField& m_domain;
  const bool m_opposite_direction;
  Point3df m_previous_point;
  float m_surface;
  float m_volume;
  VolumeRef<float>& m_volume_result;
  VolumeRef<float>& m_surface_result;
  bool m_abort;
  std::vector<float> m_voxel_size;
};

class EuclideanAdvection : public yl::Advection::Visitor
{
public:
  EuclideanAdvection(const yl::ScalarField& domain,
                     VolumeRef<float>& length_result)
    : m_domain(domain),
      m_previous_point(),
      m_length(0.f),
      m_abort(false),
      m_length_result(length_result),
      m_voxel_size(length_result->getVoxelSize())
  {
  };

  void first(const Point3df& point)
  {
    m_previous_point = point;
  };

  void visit(const Point3df& point)
  {
    const float step = (point - m_previous_point).norm();
    m_length += step;
    m_previous_point = point;
  };

  bool move_on(const Point3df& point) const
  {
    if(m_abort)
      return false;
    try {
      return m_domain.evaluate(point) >= 0.5f;
    } catch(const yl::Field::UndefinedField&) {
      return false;
    }
  };

  void finished(const Point3df& start_point)
  {
    Point3d sp_int = Point3d(
      lrint(start_point[0] / m_voxel_size[0]),
      lrint(start_point[1] / m_voxel_size[1]),
      lrint(start_point[2] / m_voxel_size[2]));
    m_length_result(sp_int) = length();
  }

  float length() const
  {
    if(m_abort)
      return no_value;
    return m_length;
  };

  bool aborted() const {
    return m_abort;
  };

  VolumeRef<float>& length_result() const
  {
    return m_length_result;
  }

private:
  const yl::ScalarField& m_domain;
  Point3df m_previous_point;
  float m_length;
  bool m_abort;
  VolumeRef<float>& m_length_result;
  std::vector<float> m_voxel_size;
};


template <typename T>
class ValueAdvection : public yl::Advection::Visitor
{
public:
  ValueAdvection(const yl::ScalarField& domain,
                 const VolumeRef<T>& value_seed,
                 VolumeRef<T>& value_result)
    : m_domain(domain),
      m_previous_point(),
      m_length(0.f),
      m_abort(false),
      m_value_seed(value_seed),
      m_value_result(value_result),
      m_voxel_size(value_seed->getVoxelSize())
  {
  };

  void first(const Point3df& point)
  {
    m_previous_point = point;
  };

  void visit(const Point3df& point)
  {
    m_previous_point = point;
  };

  bool move_on(const Point3df& point) const
  {
    if(m_abort)
      return false;
    try {
      return m_domain.evaluate(point) >= 0.5f;
    } catch(const yl::Field::UndefinedField&) {
      return false;
    }
  };

  void finished(const Point3df& start_point)
  {
    Point3d end_pos = Point3d(
      lrint(m_previous_point[0] / m_voxel_size[0]),
      lrint(m_previous_point[1] / m_voxel_size[1]),
      lrint(m_previous_point[2] / m_voxel_size[2]));
    Point3d pos = Point3d(
      lrint(start_point[0] / m_voxel_size[0]),
      lrint(start_point[1] / m_voxel_size[1]),
      lrint(start_point[2] / m_voxel_size[2]));
    m_value_result(pos) = m_value_seed(end_pos);
  }

  bool aborted() const {
    return m_abort;
  };

  VolumeRef<float>& value_result() const
  {
    return m_value_result;
  }

private:
  const yl::ScalarField& m_domain;
  Point3df m_previous_point;
  float m_length;
  bool m_abort;
  const VolumeRef<T>& m_value_seed;
  VolumeRef<T>& m_value_result;
  std::vector<float> m_voxel_size;
};

/** PathAdvection: record a wireframe mesh of advection "tracts"
 */
class PathAdvection : public yl::Advection::Visitor
{
public:
  PathAdvection(const yl::ScalarField& domain,
                AimsSurface<2>& path_result)
    : m_domain(domain),
      m_abort(false),
      m_first(true),
      m_path_result(path_result)
  {
  };

  void first(const Point3df& /*point*/)
  {
    m_first = true;
  };

  void visit(const Point3df& point)
  {
    m_path_result.vertex().push_back(point);
    if( m_first )
      m_first = false;
    else
    {
      uint32_t n = uint32_t(m_path_result.vertex().size()) - 1;
      m_path_result.polygon().push_back(AimsVector<uint32_t, 2>(n-1, n));
    }
  };

  bool move_on(const Point3df& point) const
  {
    if(m_abort)
      return false;
    try {
      return m_domain.evaluate(point) >= 0.5f;
    } catch(const yl::Field::UndefinedField&) {
      return false;
    }
  };

  void finished(const Point3df& /*start_point*/)
  {
  }

  bool aborted() const {
    return m_abort;
  };

  void abort()
  {
    std::cout << "path aborted.\n";
  }

  AimsSurface<2>& path_result() const
  {
    return m_path_result;
  }

private:
  const yl::ScalarField& m_domain;
  bool m_abort;
  bool m_first;
  AimsSurface<2>& m_path_result;
};



/** VisitorTraits manages instantiation of Visitor specializations, and of
    input / output data in the advect() function.
    By default it does basically nothing (no input/output data) and
    instantiates the visitor type with the domain field as argument (which will
    probably not be enough in many cases).
    It should be specialized.
*/
template <typename TVisitor>
class VisitorTraits
{
public:
  typedef Void ResultType;
  typedef Void InputType;
  static ResultType init_result(const VolumeRef<int16_t>& /*result_like*/,
                                const InputType& /*inputs*/) {}
  static inline TVisitor build_visitor(const yl::ScalarField& domain_field,
                                       const InputType& /*inputs*/,
                                       ResultType& /*result*/)
  { return TVisitor(domain_field); }
};

template <>
class VisitorTraits<EuclideanAdvection>
{
public:
  typedef VolumeRef<float> ResultType;
  typedef Void InputType;

  static ResultType init_result(const VolumeRef<int16_t>& result_like,
                                const InputType& /*inputs*/)
  {
    const int size_x = result_like.getSizeX();
    const int size_y = result_like.getSizeY();
    const int size_z = result_like.getSizeZ();
    VolumeRef<float> length_result(size_x, size_y, size_z);
    length_result->copyHeaderFrom(result_like.header());
    length_result.fill(no_value);
    return length_result;
  }

  static inline EuclideanAdvection build_visitor(
    const yl::ScalarField& domain_field,
    const InputType& /*inputs*/,
    ResultType& result)
  {
    return EuclideanAdvection(domain_field, result);
  }
};


template <>
class VisitorTraits<TubeAdvection>
{
public:
  typedef std::pair<VolumeRef<float>, VolumeRef<float> > ResultType;
  typedef std::pair<const yl::ScalarField &, bool> InputType;

  static ResultType init_result(const VolumeRef<int16_t>& result_like,
                                const InputType& /*inputs*/)
  {
    const int size_x = result_like.getSizeX();
    const int size_y = result_like.getSizeY();
    const int size_z = result_like.getSizeZ();
    carto::VolumeRef<float> surface_result(size_x, size_y, size_z);
    surface_result->copyHeaderFrom(result_like.header());
    surface_result.fill(no_value);
    carto::VolumeRef<float> volume_result(size_x, size_y, size_z);
    volume_result->copyHeaderFrom(result_like.header());
    volume_result.fill(no_value);
    return std::make_pair(volume_result, surface_result);
  }

  static inline TubeAdvection build_visitor(
    const yl::ScalarField& domain_field,
    const InputType& inputs,
    ResultType& result)
  {
    return TubeAdvection(inputs.first, domain_field, result.first,
                         result.second, inputs.second);
  }
};


template <typename T>
class VisitorTraits<ValueAdvection<T> >
{
public:
  typedef VolumeRef<T> ResultType;
  typedef VolumeRef<T> InputType;

  static ResultType init_result(const VolumeRef<int16_t>& result_like,
                                const InputType& inputs)
  {
    const int size_x = result_like.getSizeX();
    const int size_y = result_like.getSizeY();
    const int size_z = result_like.getSizeZ();
    VolumeRef<T> value_result(size_x, size_y, size_z);
    value_result->copyHeaderFrom(result_like.header());
    *value_result = *inputs;
    return value_result;
  }

  static inline ValueAdvection<T> build_visitor(
    const yl::ScalarField& domain_field,
    const InputType& inputs,
    ResultType& result)
  {
    return ValueAdvection<T>(domain_field, inputs, result);
  }
};


template <>
class VisitorTraits<PathAdvection>
{
public:
  typedef AimsSurface<2> ResultType;
  typedef Void InputType;

  static ResultType init_result(const VolumeRef<int16_t>& /*result_like*/,
                                const InputType& /*inputs*/)
  {
    AimsSurface<2> path_result;
    return path_result;
  }

  static inline PathAdvection build_visitor(
    const yl::ScalarField& domain_field,
    const InputType& /*inputs*/,
    ResultType& result)
  {
    return PathAdvection(domain_field, result);
  }
};


template <typename TDomainField>
class DomainFieldTraits
{
public:
  static TDomainField build_field(const VolumeRef<int16_t>& domain);
};


template <>
class DomainFieldTraits<yl::LinearlyInterpolatedScalarField>
{
public:
  static yl::LinearlyInterpolatedScalarField build_field(
    const VolumeRef<int16_t>& domain)
  {
    // This could be more elegant: the domain is first converted as float, then
    // fed into a scalar field to ease interpolation.
    carto::Converter<VolumeRef<int16_t>, VolumeRef<float> > conv;
    carto::VolumeRef<float> float_domain(*conv(domain));
    return yl::LinearlyInterpolatedScalarField(float_domain);
  }
};


template <>
class DomainFieldTraits<yl::BooleanScalarField>
{
public:
  static yl::BooleanScalarField build_field(
    const VolumeRef<int16_t>& domain)
  {
    return yl::BooleanScalarField(domain);
  }
};

/** Perform an advection for all seeds in a volume

    \warning TVisitor must be thread-safe, because the advection from multiple
    seeds will be performed in parallel if OpenMP is enabled. You should pay
    particular attention to this when writing the advection result in the
    visitor's finished() method.
 */
template <class TVisitor, class TAdvection=yl::ConstantStepAdvection>
typename VisitorTraits<TVisitor>::ResultType
advect(const yl::VectorField3d& advection_field,
       const VolumeRef<int16_t>& domain,
       const float max_advection_distance,
       const float step_size,
       const int verbosity,
       const typename VisitorTraits<TVisitor>::InputType & inputs,
       const yl::ScalarField & domain_field,
       const VolumeRef<int16_t>& advect_seeds_domain = VolumeRef<int16_t>())
{
  assert(max_advection_distance > 0);

  const VolumeRef<int16_t> & advect_seeds_domain2
    = advect_seeds_domain.getSizeX() <= 1 ? domain : advect_seeds_domain;

  const int size_x = domain.getSizeX();
  const int size_y = domain.getSizeY();
  const int size_z = domain.getSizeZ();

  const std::vector<float> voxel_size = domain->getVoxelSize();
  const float voxel_size_x = voxel_size[0];
  const float voxel_size_y = voxel_size[1];
  const float voxel_size_z = voxel_size[2];

  typename VisitorTraits<TVisitor>::ResultType result
    = VisitorTraits<TVisitor>::init_result(domain, inputs);

  unsigned int n_success = 0, n_aborted = 0;

  TAdvection advection(advection_field, step_size);
  advection.set_max_iter(std::ceil(max_advection_distance
                                  / std::abs(step_size)));
  advection.set_verbose(verbosity - 1);

  int slices_done = 0;
  #pragma omp parallel for schedule(dynamic)
  for(int z = 0; z < size_z; ++z)
  {
    for(int y = 0; y < size_y; ++y)
    for(int x = 0; x < size_x; ++x)
    {
      if(advect_seeds_domain2(x, y, z)) {
        const Point3df point(x * voxel_size_x,
                             y * voxel_size_y,
                             z * voxel_size_z);

        TVisitor visitor
          = VisitorTraits<TVisitor>::build_visitor(domain_field, inputs,
                                                   result);
        yl::Advection::Visitor& plain_visitor = visitor;
        const bool success = advection.visitor_advection(plain_visitor,
                                                         point);

        if(success) {
          // Each thread writes to different array elements, as a result no
          // synchronization should be needed. However, while this is safe on
          // most platforms, it does not seem to be guaranteed by the OpenMP
          // specification (OpenMP API v3.1, July 2011, p. 14, l. 16).
          #pragma omp atomic
          ++n_success;
        } else {
          #pragma omp atomic
          ++n_aborted;
        }
      }
    }

    #pragma omp atomic
    ++slices_done;

    if(verbosity) {
      #pragma omp critical(print_stderr)
      clog << "\r  " << slices_done << " / " << size_z << " slices processed. "
          << n_success << " voxels successfully advected, "
          << n_aborted << " aborted..." << flush;
    }
  }

  if(verbosity)
    clog << "\nyl::advect: "
        << n_success << " voxels successfully advected, "
        << n_aborted << " aborted." << endl;

  return result;
}


template <class TVisitor, class TAdvection=yl::ConstantStepAdvection,
          class TDomainField>
inline typename VisitorTraits<TVisitor>::ResultType
advect(const yl::VectorField3d& advection_field,
       const VolumeRef<int16_t>& domain,
       const float max_advection_distance,
       const float step_size,
       const int verbosity,
       const typename VisitorTraits<TVisitor>::InputType & inputs,
       const VolumeRef<int16_t>& advect_seeds_domain = VolumeRef<int16_t>())
{
  return advect<TVisitor, TAdvection>(
    advection_field, domain, max_advection_distance, step_size,
    verbosity, inputs,
    DomainFieldTraits<TDomainField>::build_field(domain),
    advect_seeds_domain);
}

} // end of anonymous namespace

namespace yl
{

template <class TDomainField>
std::pair<VolumeRef<float>, VolumeRef<float> >
advect_tubes(const yl::VectorField3d& advection_field,
             const yl::ScalarField& divergence_field,
             const VolumeRef<int16_t>& domain,
             const float max_advection_distance,
             const float step_size,
             const int verbosity,
             const VolumeRef<int16_t>& advect_seeds_domain)
{
  bool opposite_direction = step_size < 0;
  return advect<TubeAdvection, yl::ConstantStepAdvection, TDomainField>(
    advection_field, domain,
    max_advection_distance, step_size, verbosity,
    std::pair<const yl::ScalarField&, bool>(
      divergence_field, opposite_direction),
    advect_seeds_domain);
}

std::pair<VolumeRef<float>, VolumeRef<float> >
advect_tubes(const yl::VectorField3d& advection_field,
             const yl::ScalarField& divergence_field,
             const VolumeRef<int16_t>& domain,
             const float max_advection_distance,
             const float step_size,
             const yl::ScalarField & domain_field,
             const int verbosity,
             const VolumeRef<int16_t>& advect_seeds_domain)
{
  bool opposite_direction = step_size < 0;
  return advect<TubeAdvection, yl::ConstantStepAdvection>(
    advection_field, domain,
    max_advection_distance, step_size, verbosity,
    std::pair<const yl::ScalarField&, bool>(
      divergence_field, opposite_direction),
    domain_field, advect_seeds_domain);
}


template <class TDomainField>
VolumeRef<float>
advect_euclidean(const yl::VectorField3d& advection_field,
                 const VolumeRef<int16_t>& domain,
                 const float max_advection_distance,
                 const float step_size,
                 const int verbosity,
                 const VolumeRef<int16_t>& advect_seeds_domain)
{
  return advect<EuclideanAdvection, yl::ConstantStepAdvection,
                TDomainField>(
    advection_field, domain,
    max_advection_distance,
    step_size, verbosity, Void(),
    advect_seeds_domain);
}

VolumeRef<float>
advect_euclidean(const yl::VectorField3d& advection_field,
                 const VolumeRef<int16_t>& domain,
                 const float max_advection_distance,
                 const float step_size,
                 const yl::ScalarField & domain_field,
                 const int verbosity,
                 const VolumeRef<int16_t>& advect_seeds_domain)
{
  return advect<EuclideanAdvection, yl::ConstantStepAdvection>(
    advection_field, domain,
    max_advection_distance,
    step_size, verbosity, Void(),
    domain_field, advect_seeds_domain);
}

template <typename T, class TDomainField>
VolumeRef<T>
advect_value(const yl::VectorField3d& advection_field,
             const VolumeRef<T> & value_seeds,
             const VolumeRef<int16_t>& domain,
             const float max_advection_distance,
             const float step_size,
             const int verbosity,
             const VolumeRef<int16_t>& advect_seeds_domain)
{
  return advect<ValueAdvection<T>, yl::ConstantStepAdvection,
                TDomainField>(
    advection_field, domain,
    max_advection_distance,
    step_size, verbosity, value_seeds,
    advect_seeds_domain);
}

template <typename T>
VolumeRef<T>
advect_value(const yl::VectorField3d& advection_field,
             const VolumeRef<T> & value_seeds,
             const VolumeRef<int16_t>& domain,
             const float max_advection_distance,
             const float step_size,
             const yl::ScalarField & domain_field,
             const int verbosity,
             const VolumeRef<int16_t>& advect_seeds_domain)
{
  return advect<ValueAdvection<T>, yl::ConstantStepAdvection>(
    advection_field, domain,
    max_advection_distance,
    step_size, verbosity, value_seeds,
    domain_field, advect_seeds_domain);
}


template <class TDomainField>
AimsSurface<2>
advect_path(const yl::VectorField3d& advection_field,
            const carto::VolumeRef<int16_t>& domain,
            float max_advection_distance,
            float step_size,
            int verbosity,
            const carto::VolumeRef<int16_t>& advect_seeds_domain)
{
  return advect<PathAdvection, yl::ConstantStepAdvection, TDomainField>(
    advection_field, domain,
    max_advection_distance,
    step_size, verbosity,
    advect_seeds_domain);
}


AimsSurface<2>
advect_path(const yl::VectorField3d& advection_field,
            const carto::VolumeRef<int16_t>& domain,
            float max_advection_distance,
            float step_size,
            const yl::ScalarField & domain_field,
            int verbosity,
            const carto::VolumeRef<int16_t>& advect_seeds_domain)
{
  return advect<PathAdvection, yl::ConstantStepAdvection>(
    advection_field, domain,
    max_advection_distance,
    step_size, verbosity, Void(),
    domain_field, advect_seeds_domain);
}


template <class TDomainField>
yl::ScalarField*
create_domain_field(const carto::VolumeRef<int16_t>& domain)
{
  return
    new TDomainField(DomainFieldTraits<TDomainField>::build_field(domain));
}


template
VolumeRef<int16_t>
advect_value(const yl::VectorField3d& advection_field,
             const VolumeRef<int16_t> & value_seeds,
             const VolumeRef<int16_t>& domain,
             const float max_advection_distance,
             const float step_size,
             const int verbosity,
             const VolumeRef<int16_t>& advect_seeds_domain);

template
VolumeRef<int16_t>
advect_value<int16_t, yl::BooleanScalarField>(
             const yl::VectorField3d& advection_field,
             const VolumeRef<int16_t> & value_seeds,
             const VolumeRef<int16_t>& domain,
             const float max_advection_distance,
             const float step_size,
             const int verbosity,
             const VolumeRef<int16_t>& advect_seeds_domain);
template
VolumeRef<int16_t>
advect_value(const yl::VectorField3d& advection_field,
             const VolumeRef<int16_t> & value_seeds,
             const VolumeRef<int16_t>& domain,
             const float max_advection_distance,
             const float step_size,
             const yl::ScalarField & domain_field,
             const int verbosity,
             const VolumeRef<int16_t>& advect_seeds_domain);

template
VolumeRef<float>
advect_value(const yl::VectorField3d& advection_field,
             const VolumeRef<float> & value_seeds,
             const VolumeRef<int16_t>& domain,
             const float max_advection_distance,
             const float step_size,
             const int verbosity,
             const VolumeRef<int16_t>& advect_seeds_domain);

template
VolumeRef<float>
advect_value<float, yl::BooleanScalarField>(
             const yl::VectorField3d& advection_field,
             const VolumeRef<float> & value_seeds,
             const VolumeRef<int16_t>& domain,
             const float max_advection_distance,
             const float step_size,
             const int verbosity,
             const VolumeRef<int16_t>& advect_seeds_domain);
template
VolumeRef<float>
advect_value(const yl::VectorField3d& advection_field,
             const VolumeRef<float> & value_seeds,
             const VolumeRef<int16_t>& domain,
             const float max_advection_distance,
             const float step_size,
             const yl::ScalarField & domain_field,
             const int verbosity,
             const VolumeRef<int16_t>& advect_seeds_domain);

template
yl::ScalarField*
create_domain_field(const carto::VolumeRef<int16_t>& domain);

template
yl::ScalarField*
create_domain_field<yl::BooleanScalarField>(
  const carto::VolumeRef<int16_t>& domain);

} // namespace yl
