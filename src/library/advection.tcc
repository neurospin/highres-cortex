/*
Copyright CEA (2014, 2017).
Copyright Université Paris XI (2014).

Contributor: Yann Leprince <yann.leprince@ylep.fr>.
Contributor: Denis Rivière <denis.riviere@cea.fr>.

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
  Point3df start_point = point;
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

  visitor.finished(start_point);

  if(m_verbose >= 2) {
    clog << "  advection finished after " << iter << " iterations" << endl;
  }
  return !visitor.aborted();
}
