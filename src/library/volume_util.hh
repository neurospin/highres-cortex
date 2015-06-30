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

#ifndef VOLUME_UTIL_HH_INCLUDED
#define VOLUME_UTIL_HH_INCLUDED

#include <vector>

#include <cartodata/volume/volume.h>

namespace yl
{

/** Minimum border width along X, Y, and Z dimensions
 *
 * Call this function like this, or use the overload xyz_min_border(const carto::VolumeRef<T>&):
 * \code
 * yl::xyz_min_border(volume->getBorders())
 * \endcode
 */
int xyz_min_border(const std::vector<int>& borders);

/** Minimum border width along X, Y, and Z dimensions */
template <typename T>
int xyz_min_border(const carto::VolumeRef<T>& volume);

template <typename T, class Predicate>
bool
check_border_values(const carto::VolumeRef<T>& volume,
                    const Predicate& predicate);

} // namespace yl

#endif // !defined(VOLUME_UTIL_HH_INCLUDED)
