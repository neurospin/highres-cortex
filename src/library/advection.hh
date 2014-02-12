#ifndef YL_ADVECTION_HH_INCLUDED
#define YL_ADVECTION_HH_INCLUDED

#include <yleprince/field.hh>

namespace yl
{

/** Advect a visitor along a vector field */
class Advection
{
public:
  Advection(const VectorField3d& advection_field);
  virtual ~Advection() {};

  /** Advect a visitor along a vector field

      \arg visitor object being taken through all the advection path, see
      Visitor. Can be derived from Visitor, or must have all the same methods.
      \arg start_point the start of the advection path

      \retval true if the advection is stopped by the visitor
      \retval false if the advection stops for another reason
   */
  template <class TVisitor>
  bool
  visitor_advection(TVisitor& visitor,
                    Point3df start_point) const;

  /** Abstract base class for advection visitors */
  class Visitor
  {
  public:
    Visitor() {};
    virtual ~Visitor() {};

    /** Called on the first point of the advection path */
    virtual void first(const Point3df& point) = 0;
    /** Called for every point of the advection path, except the first */
    virtual void visit(const Point3df& point) = 0;
    /** Predicate that decides if the advection stops */
    virtual bool move_on(const Point3df& point) const = 0;
    /** Called after the advection is stopped by move_on() */
    virtual void finished() {};
    /** Called when the advection cannot finish successfully */
    virtual void abort() {};
  };

  /** Default iteration limit */
  static const unsigned int default_max_iter = 1000;
  /** Set the maximum number of iterations */
  void set_max_iter(unsigned int max_iter)
  {
    m_max_iter = max_iter;
  };

  /** Set the verbosity level (output to stderr)

      0: silent (default on instance creation)
      1: errors
      2: informational messages
      >=3: debug (needs setting debug_output in advection.tcc)
   */
  void set_verbose(const int verbosity)
  {
    m_verbose = verbosity;
  };

  /** Thrown by move_one_step() to abort the advection */
  class AbnormalField
  {
  };

private:
  virtual void move_one_step(Point3df& point,
                             const Point3df& local_field) const = 0;

  const yl::VectorField3d& m_vector_field;
  int m_verbose;
  unsigned int m_max_iter;
};

/** Forward advection using a constant step size

    The step size can be negative, in which case the advection is done in the
    opposite direction.
 */
class ConstantStepAdvection : public Advection
{
public:
  ConstantStepAdvection(const yl::VectorField3d& advection_field, float step);
  virtual ~ConstantStepAdvection() {};

  /** Change the step length */
  void set_step(float step)
  {
    m_step = step;
  };

private:
  virtual void move_one_step(Point3df& point,
                             const Point3df& local_field) const;

  float m_step;
};

} // namespace yl

#endif // !defined(YL_ADVECTION_HH_INCLUDED)
