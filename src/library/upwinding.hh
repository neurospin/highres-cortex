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

#ifndef UPWINDING_HH_INCLUDED
#define UPWINDING_HH_INCLUDED

#include <cartodata/volume/volume.h>

namespace yl
{

enum {
  DEFAULT_FRONT_LABEL = -123,
  DEFAULT_DONE_LABEL = -10
};

/** Compute the distance along the gradient of a scalar field.

  In the object defined by domain_label, the distance \f$d\f$ to origin_label
  is computed along the field lines of the gradient of upwind_field (\f$T\f$).
  This is done by integrating the following equation, in a single sweep through
  the domain, using an upwinding condition on \f$T\f$:

  \f[
  \nabla d \cdot \frac{\nabla T}{\|\nabla T\|} = 1
  \f]

  The domain image will be modified in-place for keeping track of the front
  (front_label) and visited voxels (done_label).

  Preconditions:
   - upwind_field must have at least a 1-voxel border, which must be
     filled with NaN
   - domain must have at least a 1-voxel border, which will must not contain
     domain_label
   - domain_label, origin_label, done_label, and front_label must all have
     different values
*/
carto::VolumeRef<float>
upwind_distance(const carto::VolumeRef<float> upwind_field,
                carto::VolumeRef<int16_t> domain,
                int16_t domain_label,
                int16_t origin_label,
                int16_t done_label = DEFAULT_DONE_LABEL,
                int16_t front_label = DEFAULT_FRONT_LABEL);

} // namespace yl

#endif // !defined(UPWINDING_HH_INCLUDED)
