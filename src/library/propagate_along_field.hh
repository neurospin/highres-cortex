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

#ifndef YL_PROPAGATE_ALONG_FIELD_HH
#define YL_PROPAGATE_ALONG_FIELD_HH

#include <boost/shared_ptr.hpp>

#include <cartodata/volume/volume.h>
#include <aims/vector/vector.h>
#include "field.hh"

namespace yl
{

/** Propagate labels along a vector field.

    From a given starting point, the field is followed in fixed-size steps
    until a seed is reached, or the maximum number of iterations is exceeded.
 */
class PropagateAlongField
{
public:
  PropagateAlongField(const boost::shared_ptr<VectorField3d>& vector_field);
  /** Provide the field as separate scalar components.

      The field will be linearly interpolated.
   */
  PropagateAlongField(const carto::VolumeRef<float>& fieldx,
                      const carto::VolumeRef<float>& fieldy,
                      const carto::VolumeRef<float>& fieldz);
  virtual ~PropagateAlongField();

  /** Indicate progress on stderr.

      - 0: be silent
      - 1: show progress and statistics
      - 2: show aborted field ascensions
      - more: for debugging, see source code and #debug_output.
   */
  void setVerbose(int=1);

  /** Move in steps of the specified size (millimetres).

      The default step size is #default_step.
   */
  void setStep(float);

  /** Abort following the gradient after a number of iterations.

      This is to prevent infinite looping when the propagation falls in a local
      minimum or loops over itself. By default this is #default_max_iter.
   */
  void setMaxIter(unsigned int);

  /** Move along the field until a nonzero label is found.

      \arg \a start_point is to be specified in millimetre coordinates.
      \arg \a seeds is interpreted as for propagate_regions().
      \arg \a ignore_label counts as propagation region.

      Return the value of the seed encountered, or zero in the following cases:
      - the path leads outside of the image;
      - the maximum number of iterations is exceeded (e.g. the path loops);
      - a vector field value is encountered that is too close to zero, or
        infinite, or NaN.
   */
  template<typename Tlabel>
  Tlabel
  ascend_until_nonzero(const Point3df &start_point,
                       const carto::VolumeRef<Tlabel> &seeds,
                       Tlabel ignore_label) const;

  /** Propagate labelled regions down a vector field.

      The labels are interpreted as follows:
      - zero: propagation zone,
      - positive: seed, propagation is stopped and the target voxel is given
        the seed's value (except if the label is \a target_label),
      - negative: forbidden region, propagation is stopped and the target voxel
        remains unset.

        An optional \a target_label can be given, in which case the propagation
        is started only from thus labelled voxels.
   */
  template<typename Tlabel>
  carto::VolumeRef<Tlabel>
  propagate_regions(const carto::VolumeRef<Tlabel> &seeds,
                    Tlabel target_label=0) const;

  /** Propagate regions, keeping track of the provenance of each voxel.

      \see propagate_regions()

      The second output volume is the same size as the input \a seeds volume, except
      is has 3 time points which represent x, y, and z coordinates of the point
      of origin, respectively.
   */
  template<typename Tlabel>
  std::pair<carto::VolumeRef<Tlabel>, carto::VolumeRef<float> >
  propagate_regions_keeping_dests(const carto::VolumeRef<Tlabel> &seeds,
                                  Tlabel target_label=0) const;


  static const float default_step = 0.1f;
  static const unsigned int default_max_iter = 1000;

  /** Visit each point along the advection path

      \retval true if the advection succeeds (reaches a goal region).
      \retval false if the advection it is stopped early and fails.
   */
  template <class TVisitor, typename Tlabel>
  bool
  visitor_ascension(const Point3df &start_point,
                    const carto::VolumeRef<Tlabel> &seeds,
                    Tlabel ignore_label,
                    TVisitor& visitor) const;

  /** Return the integral of a scalar field along an advection path.

      \exception AbortedAdvection is thrown if the advection cannot be
      completed.
   */
  template <typename Tlabel>
  float
  integrate_field_along_advection(const Point3df &start_point,
                                  const carto::VolumeRef<Tlabel> &seeds,
                                  const boost::shared_ptr<ScalarField>& field,
                                  Tlabel ignore_label) const;


  /** Exception raised by certain members if an advection cannot be completed */
  class AbortedAdvection
  {
  };

  template<typename Tlabel>
  carto::VolumeRef<float>
  integrate_from_region(const carto::VolumeRef<Tlabel> &seeds,
                        const boost::shared_ptr<ScalarField>& field,
                        Tlabel target_label=0) const;

private:
  boost::shared_ptr<VectorField3d> m_vector_field;
  float m_voxel_size_x, m_voxel_size_y, m_voxel_size_z;
  float m_invsize_x, m_invsize_y, m_invsize_z;
  unsigned int m_max_iter;
  float m_step;
  int m_verbose;

  static const int debug_output = 0;

  template <template <typename> class ResultChooserTemplate, typename Tlabel>
  typename ResultChooserTemplate<Tlabel>::result_type
  internal_ascension(const Point3df &start_point,
                     const carto::VolumeRef<Tlabel> &seeds,
                     Tlabel ignore_label,
                     const ResultChooserTemplate<Tlabel>& result_chooser) const;

  template <class ResultChooser, typename Tlabel>
  typename ResultChooser::result_type
  internal_ascension(const Point3df &start_point,
                     const carto::VolumeRef<Tlabel> &seeds,
                     Tlabel ignore_label,
                     const ResultChooser& result_chooser) const;

  template<class ResultRecorder, typename Tlabel>
  void
  internal_propagation(const carto::VolumeRef<Tlabel> &seeds,
                       const Tlabel target_label,
                       ResultRecorder& result_recorder) const;
};

};

#endif // !defined(YL_PROPAGATE_ALONG_FIELD_HH)
