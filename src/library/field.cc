#include <yleprince/field.hh>

yl::LinearlyInterpolatedVectorField3d::
LinearlyInterpolatedVectorField3d(const carto::VolumeRef<float>& fieldx,
                                  const carto::VolumeRef<float>& fieldy,
                                  const carto::VolumeRef<float>& fieldz)
  : m_interp_fieldx(fieldx), m_interp_fieldy(fieldy), m_interp_fieldz(fieldz)
{
}

void
yl::LinearlyInterpolatedVectorField3d::
evaluate(const Point3df& pos, Point3df& output) const
{
  if(! m_interp_fieldx.isValid(pos[0], pos[1], pos[2])
     || ! m_interp_fieldy.isValid(pos[0], pos[1], pos[2])
     || ! m_interp_fieldz.isValid(pos[0], pos[1], pos[2]))
  {
    throw UndefinedField();
  }
  output[0] = m_interp_fieldx.value(pos);
  output[1] = m_interp_fieldy.value(pos);
  output[2] = m_interp_fieldz.value(pos);
}

yl::LinearlyInterpolatedScalarField::
LinearlyInterpolatedScalarField(const carto::VolumeRef<float>& field_volume)
  : m_interp_field(field_volume)
{
}

float
yl::LinearlyInterpolatedScalarField::
evaluate(const Point3df& pos) const
{
  if(! m_interp_field.isValid(pos[0], pos[1], pos[2]))
  {
    throw UndefinedField();
  }
  return m_interp_field.value(pos);
}
