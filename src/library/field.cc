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

#include "field.hh"

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
