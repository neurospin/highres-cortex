#ifndef YL_PROPAGATE_ALONG_FIELD_HH
#define YL_PROPAGATE_ALONG_FIELD_HH

#include <cartodata/volume/volume.h>
#include <aims/vector/vector.h>
#include <aims/resampling/linearInterpolator.h>

namespace yl
{

  /** Propagate labels along a vector field.

      From a given starting point, the field is followed in fixed-size steps
      until a seed is reached, or the maximum number of iterations is exceeded.
   */
class PropagateAlongField
{
public:
  /** Provide the field as separate scalar components.
   */
  PropagateAlongField(const carto::VolumeRef<float> &fieldx,
                      const carto::VolumeRef<float> &fieldy,
                      const carto::VolumeRef<float> &fieldz);
  virtual ~PropagateAlongField();

  /** Indicate progress on stderr.

      - 0: be silent
      - 1: show progress and statistics
      - 2: show aborted field ascensions
      - more: for debugging, see source code and #debug_output.
   */
  void setVerbose(int = 1);

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

      Return the value of the seed encontered, or zero in the following cases:
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

  static const float default_step = 0.1f;
  static const unsigned int default_max_iter = 1000;

private:
  aims::LinearInterpolator<float> m_interp_fieldx;
  aims::LinearInterpolator<float> m_interp_fieldy;
  aims::LinearInterpolator<float> m_interp_fieldz;
  float m_voxel_size_x, m_voxel_size_y, m_voxel_size_z;
  float m_invsize_x, m_invsize_y, m_invsize_z;
  unsigned int m_max_iter;
  float m_step;
  int m_verbose;

  static const int debug_output = 0;
};

};

#endif // !defined(YL_PROPAGATE_ALONG_FIELD_HH)
