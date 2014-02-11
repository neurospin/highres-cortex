#ifndef YL_ISOVOLUME_HH_INCLUDED
#define YL_ISOVOLUME_HH_INCLUDED

#include <cartodata/volume/volume.h>

namespace yl
{

class VectorField3d;
class ScalarField;

/** Advect a tube along a field, starting with unit surface

    \arg advection_field vector field to advect along
    \arg divergence_field the divergence of the normalized advection field
    \arg domain the advection domain with zero outside, one inside
    \arg max_advection_distance the maximum length of the advection path
    \arg step_size the constant length of an advection step
    \arg verbosity verbosity to stderr, (verbosity - 1) is passed to
    Advection::set_verbose()

    \return pair of (tube's volume, tube's end surface)
 */
std::pair<carto::VolumeRef<float>, carto::VolumeRef<float> >
advect_tubes(const yl::VectorField3d& advection_field,
             const yl::ScalarField& divergence_field,
             const carto::VolumeRef<int16_t>& domain,
             float max_advection_distance,
             float step_size,
             int verbosity=0);

} // namespace yl

#endif // !defined(YL_ISOVOLUME_HH_INCLUDED)
