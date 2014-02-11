#include <yleprince/advection.hh>
#include <cmath>

yl::Advection::Advection(const VectorField3d& advection_field)
  : m_vector_field(advection_field), m_verbose(0), m_max_iter(default_max_iter)
{
}

yl::ConstantStepAdvection::
ConstantStepAdvection(const yl::VectorField3d& advection_field,
                      const float step)
  : Advection(advection_field), m_step(step)
{
}

void yl::ConstantStepAdvection::
move_one_step(Point3df& point,
              const Point3df& local_field) const
{
  float gx = local_field[0], gy = local_field[1], gz = local_field[2];
  // Normalize the field, stop if too small or infinite or NaN.
  const float gn = std::sqrt(gx*gx + gy*gy + gz*gz);
  if(!std::isnormal(gn))
    throw yl::Advection::AbnormalField();
  gx /= gn; gy /= gn; gz /= gn;

  point[0] += m_step * gx;
  point[1] += m_step * gy;
  point[2] += m_step * gz;
}

#include <yleprince/advection.tcc>

template bool yl::Advection::visitor_advection<yl::Advection::Visitor>
(yl::Advection::Visitor& visitor, Point3df starst_point) const;
