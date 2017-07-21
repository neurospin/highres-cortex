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

class LinearlyInterpolatedScalarFieldGradient : public VectorField3d
{
public:
  LinearlyInterpolatedScalarFieldGradient(
    const carto::VolumeRef<float>& scalar_field);
  virtual void evaluate(const Point3df& pos, Point3df& output) const;
private:
  carto::rc_ptr<aims::Interpolator> m_interp_gradx;
  carto::rc_ptr<aims::Interpolator> m_interp_grady;
  carto::rc_ptr<aims::Interpolator> m_interp_gradz;
  float m_xoffset, m_yoffset, m_zoffset;
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

/** Access a scalar field stored in a volume

    The field's value is linearly interpolated between integer coordinates.
 */
class BooleanScalarField : public ScalarField
{
public:
  BooleanScalarField(const carto::VolumeRef<int16_t>& field_volume);

  virtual float evaluate(const Point3df& pos) const;
private:
  const carto::VolumeRef<int16_t>& m_field;
  Point3df m_voxel_size;
};

}

#endif // !defined(YL_FIELD_HH_INCLUDED)
