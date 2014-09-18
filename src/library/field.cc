#include <yleprince/field.hh>

#include <aims/math/gradient.h>

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

yl::LinearlyInterpolatedScalarFieldGradient::
LinearlyInterpolatedScalarFieldGradient(const carto::VolumeRef<float>& scalar_field)
  : m_interp_gradx(), m_interp_grady(), m_interp_gradz()
{
  AimsGradient<float> gradient(AIMS_GRADIENT_DPLUS);
  AimsData<float> gradx = gradient.X(scalar_field);
  m_interp_gradx = getLinearInterpolator(gradx);
  AimsData<float> grady = gradient.Y(scalar_field);
  m_interp_grady = getLinearInterpolator(grady);
  AimsData<float> gradz = gradient.Z(scalar_field);
  m_interp_gradz = getLinearInterpolator(gradz);

  // Problem: the last line of each gradient is filled with zeros, which
  // means that the interpolation is meaningless near these borders.
  // This is worked around non-elegantly in the evaluate method.

  const std::vector<float> voxel_size = scalar_field->getVoxelSize();
  m_xoffset = -0.5f * voxel_size[0];
  m_yoffset = -0.5f * voxel_size[1];
  m_zoffset = -0.5f * voxel_size[2];
}

void
yl::LinearlyInterpolatedScalarFieldGradient::
evaluate(const Point3df& pos, Point3df& output) const
{
  if(! m_interp_gradx->isValid(pos[0] + m_xoffset, pos[1], pos[2])
     || ! m_interp_grady->isValid(pos[0], pos[1] + m_yoffset, pos[2])
     || ! m_interp_gradz->isValid(pos[0], pos[1], pos[2] + m_zoffset))
  {
    throw UndefinedField();
  }

  // Workaround for problem described in the constructor: test if we are near
  // the border. Should work, but not elegant...
  if(! m_interp_gradx->isValid(pos[0] - m_xoffset, pos[1], pos[2])
     || ! m_interp_grady->isValid(pos[0], pos[1] - m_yoffset, pos[2])
     || ! m_interp_gradz->isValid(pos[0], pos[1], pos[2] - m_zoffset))
  {
    throw UndefinedField();
  }

  output[0] = m_interp_gradx->value(pos[0] + m_xoffset, pos[1], pos[2]);
  output[1] = m_interp_grady->value(pos[0], pos[1] + m_yoffset, pos[2]);
  output[2] = m_interp_gradz->value(pos[0], pos[1], pos[2] + m_zoffset);
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
