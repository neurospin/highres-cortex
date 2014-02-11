#include <iostream>

using std::clog;
using std::endl;

namespace
{
int debug_output = 0;
}

template <class TVisitor>
bool
yl::Advection::
visitor_advection(TVisitor& visitor,
                  Point3df point) const
{
  if(m_verbose >= 2) {
    clog << "yl::Advection::visitor_advection starting at " << point << endl;
  }

  unsigned int iter = 0;
  visitor.first(point);

  while(visitor.move_on(point)) {
    if(iter >= m_max_iter) {
      if(m_verbose >= 1) {
        clog << "  advection aborted after " << iter
             << " iterations: too many iterations." << endl;
      }
      visitor.abort();
      return false;
    }

    if(debug_output >= 3 && m_verbose >= 3) {
      clog << "  iteration " << iter << " at " << point << endl;
    }

    ++iter;
    visitor.visit(point);

    // Move along the field
    Point3df local_field;
    try {
      m_vector_field.evaluate(point, local_field);
    } catch(const Field::UndefinedField&) {
      if(m_verbose >= 1) {
        clog << "  advection aborted after " << iter
             << " iterations: vector field is undefined at " << point << endl;
      }

      visitor.abort();
      return false;
    }

    try {
      move_one_step(point, local_field);
    } catch(const AbnormalField&) {
      if(m_verbose >= 1) {
        clog << "  advection aborted after " << iter
             << " iterations: vector field is abnormal"
                " (too small, infinite, or NaN) at "
             << point << endl;
      }

      visitor.abort();
      return false;
    }
  }

  visitor.finished();

  if(m_verbose >= 2) {
    clog << "  advection finished after " << iter << " iterations" << endl;
  }
  return true;
}
