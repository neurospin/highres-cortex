#ifndef YL_VECTOR_FIELD_HH_INCLUDED
#define YL_VECTOR_FIELD_HH_INCLUDED

#include <aims/vector/vector.h>
#include <aims/resampling/linearInterpolator.h>

namespace yl
{

/** Provide access to a vector field at possibly non-integer coordinates */
class VectorField
{
public:
  /** Exception thrown when the field cannot be evaluated */
  class UndefinedField
  {
  };


  /** Evaluate the vector field at possibly non-integer coordinates

      \exception UndefinedVectorField if the field cannot be evaluated at the
      given position.
   */
  Point3df evaluate(const Point3df& pos) const
  {
    Point3df ret;
    evaluate(pos, ret);
    return ret;
  };

  /** Evaluate the vector field at possibly non-integer coordinates

      \exception UndefinedVectorField if the field cannot be evaluated at the
      given position.
  */
  virtual void evaluate(const Point3df& pos, Point3df& output) const = 0;
};

/** Access a vector field stored as three volumes

    The components are linearly interpolated between integer coordinates.
 */
class LinearlyInterpolatedVectorField3D : public VectorField
{
public:
  LinearlyInterpolatedVectorField3D(const carto::VolumeRef<float>& fieldx,
                                    const carto::VolumeRef<float>& fieldy,
                                    const carto::VolumeRef<float>& fieldz);

  virtual void evaluate(const Point3df& pos, Point3df& output) const;
private:
  aims::LinearInterpolator<float> m_interp_fieldx;
  aims::LinearInterpolator<float> m_interp_fieldy;
  aims::LinearInterpolator<float> m_interp_fieldz;
};

}

#endif // !defined(YL_VECTOR_FIELD_HH_INCLUDED)
