#ifndef YL_FIELD_HH_INCLUDED
#define YL_FIELD_HH_INCLUDED

#include <aims/vector/vector.h>
#include <aims/resampling/linearInterpolator.h>

namespace yl
{

/** Store and access a field (scalar or vector) defined over a 3D domain */
class Field
{
public:
  /** Exception thrown when the field has no defined value */
  class UndefinedField
  {
  };
};

/** Store a vector field and access it at any coordinates */
class VectorField3d : public Field
{
public:
  /** Evaluate the field's value at possibly non-integer coordinates

      \exception UndefinedField if the field cannot be evaluated at the
      given position.
   */
  Point3df evaluate(const Point3df& pos) const
  {
    Point3df ret;
    evaluate(pos, ret);
    return ret;
  };

  /** Evaluate the field's value at possibly non-integer coordinates

      \exception UndefinedField if the field cannot be evaluated at the
      given position.
  */
  virtual void evaluate(const Point3df& pos, Point3df& output) const = 0;
};

/** Store a scalar field and access it at any coordinates */
class ScalarField : public Field
{
public:
  /** Evaluate the field's value at possibly non-integer coordinates

      \exception UndefinedField if the field cannot be evaluated at the
      given position.
   */
  virtual float evaluate(const Point3df& pos) const = 0;

};

/** Access a vector field stored as three volumes

    The components are linearly interpolated between integer coordinates.
 */
class LinearlyInterpolatedVectorField3d : public VectorField3d
{
public:
  LinearlyInterpolatedVectorField3d(const carto::VolumeRef<float>& fieldx,
                                    const carto::VolumeRef<float>& fieldy,
                                    const carto::VolumeRef<float>& fieldz);

  virtual void evaluate(const Point3df& pos, Point3df& output) const;
private:
  aims::LinearInterpolator<float> m_interp_fieldx;
  aims::LinearInterpolator<float> m_interp_fieldy;
  aims::LinearInterpolator<float> m_interp_fieldz;
};

/** Access a scalar field stored in a volume

    The field's value is linearly interpolated between integer coordinates.
 */
class LinearlyInterpolatedScalarField : public ScalarField
{
public:
  LinearlyInterpolatedScalarField(const carto::VolumeRef<float>& field_volume);

  virtual float evaluate(const Point3df& pos) const;
private:
  aims::LinearInterpolator<float> m_interp_field;
};

}

#endif // !defined(YL_FIELD_HH_INCLUDED)
