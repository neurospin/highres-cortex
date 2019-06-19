/*
Copyright Forschungszentrum Jülich GmbH (2017).
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

#ifndef YL_ISOVOLUME_HH_INCLUDED
#define YL_ISOVOLUME_HH_INCLUDED

#include <cartodata/volume/volume.h>
#include <aims/mesh/surface.h>

namespace yl
{

class VectorField3d;
class ScalarField;
class LinearlyInterpolatedScalarField;

/** Advect a tube along a field, starting with unit surface

    This non-template function differs from the other (template) one in that it
    takes the domain field as an additional, dynamic argument.

    \arg \p seeds : mask of advection starting points
    \arg \p domain : the advection domain with zero outside, one inside
    \arg \p advection_field : vector field to advect along
    \arg \p divergence_field : the divergence of the normalized advection field
    \arg \p max_advection_distance : the maximum length of the advection path
    \arg \p step_size : the constant length of an advection step. Can be
      negative to advect in the opposite direction.
    \arg \p verbosity : verbosity to stderr, (verbosity - 1) is passed to
      Advection::set_verbose() and to Advection::Visitor::set_verbose()

    \return pair of (tube's volume, tube's end surface)
 */
std::pair<carto::VolumeRef<float>, carto::VolumeRef<float> >
advect_tubes(const carto::VolumeRef<int16_t>& seeds,
             const yl::ScalarField& domain,
             const yl::VectorField3d& advection_field,
             const yl::ScalarField& divergence_field,
             float max_advection_distance,
             float step_size,
             int verbosity=0);

/** Advect a point along a field, keeping track of the distance

    This non-template function differs from the other (template) one in that it
    takes the domain field as an additional, dynamic argument.

    \arg \p seeds : mask of advection starting points mask
    \arg \p domain : the advection domain with zero outside, one inside
    \arg \p advection_field : vector field to advect along
    \arg \p max_advection_distance : the maximum length of the advection path
    \arg \p step_size : the constant length of an advection step. Can be
      negative to advect in the opposite direction.
    \arg \p verbosity : verbosity to stderr, (verbosity - 1) is passed to
      Advection::set_verbose() and to Advection::Visitor::set_verbose()

    \return Euclidean length of the advection path
 */
carto::VolumeRef<float>
advect_euclidean(const carto::VolumeRef<int16_t>& seeds,
                 const yl::ScalarField& domain,
                 const yl::VectorField3d& advection_field,
                 float max_advection_distance,
                 float step_size,
                 int verbosity=0);

/** Advect a point along a field, and propagate end points values to the
    starting point

    This function differs from the other one (with 2 template parameters) in
    that it takes the domain field as an additional, dynamic argument.

    \p T : values image voxel type

    \arg \p seeds : mask of advection starting points
    \arg \p domain : the advection domain with zero outside, one inside
    \arg \p advection_field : vector field to advect along
    \arg \p values : values to be propagated
    \arg \p max_advection_distance : the maximum length of the advection path
    \arg \p step_size : the constant length of an advection step. Can be
      negative to advect in the opposite direction.
    \arg \p verbosity : verbosity to stderr, (verbosity - 1) is passed to
      Advection::set_verbose() and to Advection::Visitor::set_verbose()

    \return advected values
 */
template <typename T>
carto::VolumeRef<T>
advect_value(const carto::VolumeRef<int16_t>& seeds,
             const yl::ScalarField& domain,
             const yl::VectorField3d& advection_field,
             const carto::VolumeRef<T>& values,
             float max_advection_distance,
             float step_size,
             int verbosity=0);

/** Advect a point along a field, recording advection tracts in a wireframe
    mesh

    This non-template function differs from the other (template) one in that it
    takes the domain field as an additional, dynamic argument.

    \arg \p seeds : mask of advection starting points
    \arg \p domain : the advection domain with zero outside, one inside
    \arg \p advection_field : vector field to advect along
    \arg \p max_advection_distance : the maximum length of the advection path
    \arg \p step_size : the constant length of an advection step. Can be
      negative to advect in the opposite direction.
    \arg \p verbosity : verbosity to stderr, (verbosity - 1) is passed to
      Advection::set_verbose() and to Advection::Visitor::set_verbose()

    \return Euclidean length of the advection path
 */
AimsSurface<2>
advect_path(const carto::VolumeRef<int16_t>& seeds,
            const yl::ScalarField& domain,
            const yl::VectorField3d& advection_field,
            float max_advection_distance,
            float step_size,
            int verbosity=0);

/** Build a ScalarField of the given type from an int16_t volume
 */
template <class TDomainField>
yl::ScalarField*
create_domain_field(const carto::VolumeRef<int16_t>& domain);

} // namespace yl

#endif // !defined(YL_ISOVOLUME_HH_INCLUDED)
