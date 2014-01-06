#include <yleprince/vector_field.hh>

yl::LinearlyInterpolatedVectorField3D::
LinearlyInterpolatedVectorField3D(const carto::VolumeRef<float>& fieldx,
                                  const carto::VolumeRef<float>& fieldy,
                                  const carto::VolumeRef<float>& fieldz)
  : m_interp_fieldx(fieldx), m_interp_fieldy(fieldy), m_interp_fieldz(fieldz)
{
}

void
yl::LinearlyInterpolatedVectorField3D::
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
