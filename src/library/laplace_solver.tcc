/*
Copyright Télécom ParisTech (2015).

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

#include "laplace_solver.hh"

#include <iostream>

#include <boost/format.hpp>
#include <boost/unordered_set.hpp>
#include <boost/functional/hash.hpp>

// Fix obscure compilation error, we do not use float128 anyway
#define BOOST_MATH_DISABLE_FLOAT128
#include <boost/math/constants/constants.hpp>

#include <aims/utility/converter_volume.h>
#include <aims/connectivity/structuring_element.h>

#include "cortex.hh"
#include "front.hh"

namespace {

template <typename T>
inline T square(T x)
{
  return x * x;
}

/** Check if there is a local extremum at the given coordinates.

    The check always accesses all 6-neighbours of the given point: therefore,
    all coordinates must be between 1 and (size - 2).
 */
template <typename Real>
bool
is_local_extremum(const carto::VolumeRef<Real>& solution,
                  const int x, const int y, const int z)
{
  unsigned int num_smaller = 0;
  unsigned int num_larger = 0;

  if(solution(x, y, z) < solution(x - 1, y, z))
    ++num_smaller;
  else if(solution(x, y, z) > solution(x - 1, y, z))
    ++num_larger;

  if(solution(x, y, z) < solution(x + 1, y, z))
    ++num_smaller;
  else if(solution(x, y, z) > solution(x + 1, y, z))
    ++num_larger;

  if(solution(x, y, z) < solution(x, y - 1, z))
    ++num_smaller;
  else if(solution(x, y, z) > solution(x, y - 1, z))
    ++num_larger;

  if(solution(x, y, z) < solution(x, y + 1, z))
    ++num_smaller;
  else if(solution(x, y, z) > solution(x, y + 1, z))
    ++num_larger;

  if(solution(x, y, z) < solution(x, y, z - 1))
    ++num_smaller;
  else if(solution(x, y, z) > solution(x, y, z - 1))
    ++num_larger;

  if(solution(x, y, z) < solution(x, y, z + 1))
    ++num_smaller;
  else if(solution(x, y, z) > solution(x, y, z + 1))
    ++num_larger;

  return (num_smaller == 6 || num_larger == 6);
}

} // end of anonymous namespace


template <typename Real>
yl::LaplaceSolver<Real>::
LaplaceSolver(const carto::VolumeRef<int16_t>& classif)
  : m_classif(classif),
    m_solution(classif->getSizeX(),
               classif->getSizeY(),
               classif->getSizeZ(),
               1, s_border_width),
    m_verbosity(0)
{
  m_solution->copyHeaderFrom(classif->header());
}

template <typename Real>
void
yl::LaplaceSolver<Real>::
initialize_solution()
{
  carto::RescalerInfo rescaler_info;
  rescaler_info.vmin = yl::CSF_LABEL;  // 0
  rescaler_info.vmax = yl::WHITE_LABEL;  // 200
  rescaler_info.omin = 0;
  rescaler_info.omax = 1;

  carto::Converter<carto::VolumeRef<int16_t>, carto::VolumeRef<Real> >
    converter(true, rescaler_info);
  converter.convert(m_classif, m_solution);
}

// This is a Successive Over-Relaxation method using a Gauss-Seidel iteration
// with the largest time step ensuring stability, and even-odd traversal with
// Chebyshev acceleration. See W.H. Press, W.H Teukolsky, W.T. Vetterling and
// B.P. Flannery; Numerical Recipes in C++: The Art of Scientific Computing;
// (Cambridge University Press, 2002); pages 868--871.
template <typename Real>
void
yl::LaplaceSolver<Real>::
SOR(const Real absolute_precision,
    const float typical_cortical_thickness)
{
  static const Real pi = boost::math::constants::pi<Real>();

  assert(absolute_precision > 0);
  assert(typical_cortical_thickness > 0);

  const int size_x = m_classif->getSizeX();
  const int size_y = m_classif->getSizeY();
  const int size_z = m_classif->getSizeZ();

  const std::vector<float> voxsize = m_classif->getVoxelSize();
  const Real voxsize_x = voxsize[0];
  const Real voxsize_y = voxsize[1];
  const Real voxsize_z = voxsize[2];
  const Real average_voxsize = (voxsize_x + voxsize_y + voxsize_z) / 3;

  const Real factx = 0.5f / (1 + square(voxsize_x / voxsize_y) +
                             square(voxsize_x / voxsize_z));
  const Real facty = 0.5f / (square(voxsize_y / voxsize_x) + 1 +
                             square(voxsize_y / voxsize_z));
  const Real factz = 0.5f / (square(voxsize_z / voxsize_x) +
                             square(voxsize_z / voxsize_y) + 1);

  // Using the average of the voxel size in the estimation of cortical
  // thickness in voxels gave the best result in some tests with 3:1
  // anisotropic voxels.
  const Real cortical_voxels = typical_cortical_thickness / average_voxsize;
  const Real estimated_jacobi_spectral_radius =
    (std::cos(pi / cortical_voxels) + 2) / 3;
  const Real omega_update_factor =
    0.25f * square(estimated_jacobi_spectral_radius);

  if(m_verbosity) {
    std::clog << boost::format(
      "solving Laplace equation using Successive Over-Relaxation\n"
      "(requested absolute precision %1$.1e)...") % absolute_precision
              << std::endl;
  }

  Real max_residual;
  Real omega = 1;
  unsigned int iter = 1;
  unsigned int max_iter = 0;  // initially no limit on iteration count
  do {
    max_residual = 0;
    for(int even_odd = 0; even_odd <= 1; ++even_odd) {
      #if _OPENMP >= 201107
      // OpenMP 3.1 or later (201107) is needed for the reduction(max) clause
      #pragma omp parallel for schedule(dynamic) reduction(max:max_residual)
      #endif
      for(int z = 0 ; z < size_z ; ++z)
        for(int y = 0 ; y < size_y ; ++y)
          for(int x = (y + z + even_odd) % 2 ; x < size_x ; x += 2) {
            if(m_classif(x, y, z) == yl::CORTEX_LABEL) {
              const Real old_value = m_solution(x, y, z);
              const Real gauss_seidel_new_value =
                factx * (m_solution(x - 1, y, z) + m_solution(x + 1, y, z)) +
                facty * (m_solution(x, y - 1, z) + m_solution(x, y + 1, z)) +
                factz * (m_solution(x, y, z - 1) + m_solution(x, y, z + 1));
              m_solution(x, y, z) = omega * gauss_seidel_new_value
                + (1 - omega) * old_value;
              const Real abs_residual = std::abs(gauss_seidel_new_value
                                                 - old_value);
              if(abs_residual > max_residual) {
                max_residual = abs_residual;
              }
            }
          }

      // Chebyshev acceleration
      if(iter == 1 && even_odd == 0)
        omega = 1 / (1 - 0.5f * square(estimated_jacobi_spectral_radius));
      else
        omega = 1 / (1 - omega * omega_update_factor);
      assert(omega > 0);
      assert(omega < 2);
    }

    // The residual cannot be guaranteed to decrease below this threshold
    // because of numerical errors. The value 11 * epsilon comes from a crude
    // analysis of the propagation of numerical errors (potential infinite loop
    // BUG if this has been underestimated...)
    //
    // When we detect this, we can predict the number of iterations necessary
    // to reach machine precision, assuming geometric convergence rate.
    //
    // Possible refinement: set omega = 1 for a few last steps (how many?):
    // that should decrease numerical errors by up to a factor omega.
    static const Real numerical_error_threshold
      = 11 * std::numeric_limits<Real>::epsilon();
    if(max_residual <= numerical_error_threshold && max_iter == 0) {
      static const Real epsilon = std::numeric_limits<Real>::epsilon();
      static const Real initial_error = 0.5;

      max_iter = static_cast<unsigned int>(
        iter * (std::log(epsilon / initial_error)
                / std::log(numerical_error_threshold / initial_error)) + 1);
      // Prevent infinite looping if an error in the above calculation leads to
      // max_iter <= iter
      max_iter = std::max(max_iter, iter + 1);
      if(m_verbosity) {
        std::clog << boost::format("\nmachine precision (%1$.1e) should"
                                   " be attained by iteration %2$u")
          % epsilon % max_iter << std::endl;
      }
    }

    if(m_verbosity) {
      std::clog << boost::format("\riteration %1$u, max residual = %2$.1e")
        % iter % max_residual << std::flush;
    }

    ++iter;
  } while(max_residual > absolute_precision && iter != max_iter);

  if(m_verbosity)
    std::clog << std::endl;
}

template <typename Real>
void
yl::LaplaceSolver<Real>::
eliminate_extrema()
{
  const int size_x = m_classif->getSizeX();
  const int size_y = m_classif->getSizeY();
  const int size_z = m_classif->getSizeZ();

  const std::vector<float> voxsize = m_classif->getVoxelSize();
  const Real voxsize_x = voxsize[0];
  const Real voxsize_y = voxsize[1];
  const Real voxsize_z = voxsize[2];

  const Real factx = 0.5f / (1 + square(voxsize_x / voxsize_y) +
                              square(voxsize_x / voxsize_z));
  const Real facty = 0.5f / (square(voxsize_y / voxsize_x) + 1 +
                              square(voxsize_y / voxsize_z));
  const Real factz = 0.5f / (square(voxsize_z / voxsize_x) +
                              square(voxsize_z / voxsize_y) + 1);

  if(m_verbosity) {
    std::clog << "searching for remaining local extrema... "
              << std::flush;
  }

  typedef boost::unordered_set<Point3d, yl::PointHasher<Point3d> >
    PointSetType;
  PointSetType extremum_points;

  #pragma omp parallel for schedule(dynamic)
  for(int z = 1 ; z < size_z - 1 ; ++z)
    for(int y = 1 ; y < size_y - 1 ; ++y)
      for(int x = 1 ; x < size_x - 1 ; ++x) {
        if(m_classif(x, y, z) == yl::CORTEX_LABEL) {
          if(is_local_extremum(m_solution, x, y, z)) {
            #pragma omp critical(extremum_point_set)
            extremum_points.insert(Point3d(x, y, z));
          }
        }
      }

  if(m_verbosity)
    std::clog << extremum_points.size() << " found." << std::endl;

  if(extremum_points.empty())
    return;

  std::auto_ptr<aims::strel::Connectivity> connectivity(
    aims::strel::ConnectivityFactory::create("6"));

  while(!extremum_points.empty()) {
    PointSetType::const_iterator point_it = extremum_points.begin();
    const Point3d point = *point_it;
    extremum_points.erase(point_it);

    const int x = point[0];
    const int y = point[1];
    const int z = point[2];

    const Real new_value =
      factx * (m_solution(x - 1, y, z) + m_solution(x + 1, y, z)) +
      facty * (m_solution(x, y - 1, z) + m_solution(x, y + 1, z)) +
      factz * (m_solution(x, y, z - 1) + m_solution(x, y, z + 1));
    m_solution(x, y, z) = new_value;

    if(m_verbosity && is_local_extremum(m_solution, x, y, z)) {
      std::clog << "\nwarning: local extremum at " << point
                << " could not be eliminated" << std::endl;
    }

    for(aims::strel::Connectivity::const_iterator
          neighbour_it = connectivity->begin() ;
        neighbour_it != connectivity->end() ;
        ++neighbour_it) {
      const Point3d& neighbour_offset = Point3d((*neighbour_it)[0],
                                                (*neighbour_it)[1],
                                                (*neighbour_it)[2]);
      const Point3d neighbour_point = point + neighbour_offset;
      const int nx = neighbour_point[0];
      const int ny = neighbour_point[1];
      const int nz = neighbour_point[2];

      if(!(nx >= 1 && nx < size_x - 1
           && ny >= 1 && ny < size_y - 1
           && nz >= 1 && nz < size_z - 1)) {
        // We do not check for local extrema on the edge of the volume.
        continue;
      }

      if(m_classif(nx, ny, nz) == yl::CORTEX_LABEL
         && is_local_extremum(m_solution, nx, ny, nz)) {
        extremum_points.insert(neighbour_point);
      }
    }

    if(m_verbosity) {
      std::clog << "\reliminating local extrema... "
                << extremum_points.size() << " remaining. "
                << std::flush;
    }
  }

  if(m_verbosity)
    std::clog << std::endl;
}

template <typename Real>
void
yl::LaplaceSolver<Real>::
clamp_to_range(const Real min, const Real max)
{
  assert(max >= min);

  const int size_x = m_solution->getSizeX();
  const int size_y = m_solution->getSizeY();
  const int size_z = m_solution->getSizeZ();

  #pragma omp parallel for
  for(int z = 0 ; z < size_z ; ++z)
    for(int y = 0 ; y < size_y ; ++y)
      for(int x = 0 ; x < size_x ; ++x) {
        if(m_solution(x, y, z) < min)
          m_solution(x, y, z) = min;
        else if(m_solution(x, y, z) > max)
          m_solution(x, y, z) = max;
      }
}

template <typename Real>
carto::VolumeRef<Real>
yl::LaplaceSolver<Real>::
solution() const
{
  return m_solution;
}
