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

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <sstream>

#include <soma-io/allocator/allocator.h>
#include <cartobase/config/verbose.h>
#include <cartodata/volume/volume.h>
#include <aims/getopt/getopt2.h>
#include <aims/io/reader.h>
#include <aims/io/writer.h>

#include <highres-cortex/upwinding.hh>
#include <highres-cortex/cortex.hh>
#include <highres-cortex/volume_util.hh>

using std::clog;
using std::endl;
using carto::VolumeRef;
using carto::verbose;
using soma::AllocatorStrategy;
using soma::AllocatorContext;


// Anonymous namespace for file-local symbols
namespace
{
const int EXIT_USAGE_ERROR = 2;
std::string program_name;
} // end of anonymous namespace


int main(const int argc, const char **argv)
{
  // Initialize command-line option parsing
  aims::Reader<AimsData<float> > upwind_field_reader;
  aims::Reader<VolumeRef<float> > speed_map_reader;
  aims::Reader<AimsData<int16_t> > domain_reader;
  aims::Writer<VolumeRef<float> > output_writer;
  int16_t domain_label = yl::CORTEX_LABEL;
  int16_t origin_label = yl::CSF_LABEL;
  bool invert_field = false, opposite_speed = false;

  program_name = argv[0];
  aims::AimsApplication app(argc, argv,
"Compute the distance to a boundary along the gradient of a scalar field."
);
  app.addOption(domain_reader, "--domain",
                "label image defining the computation domain");
  app.addOption(upwind_field_reader, "--field",
                "scalar field whose gradient is used as the"
                " integration direction");
  app.addOption(invert_field, "--invert",
                "work on inverted field (downfield instead of upfield)",
                true);
  app.addOption(speed_map_reader, "--speed-map",
                "speed map of the field evolution (UNIMPLEMENTED)", true);
  app.addOption(opposite_speed, "--opposite-speed",
                "use -1 * the speed map", true);
  app.addOption(output_writer, "--output",
                "output volume containing the distance");
  {
    std::ostringstream help_str;
    help_str << "label of the propagation domain (default: "
             << domain_label << ")";
    app.addOption(domain_label, "--domain-label", help_str.str(), true);
  }
  {
    std::ostringstream help_str;
    help_str << "label of the origin object (default: "
             << origin_label << ")";
    app.addOption(origin_label, "--origin-label", help_str.str(), true);
  }

  // Process command-line options
  try
  {
    app.initialize();
  }
  catch(const carto::user_interruption &)
  {
    // Exit after printing e.g. help
    return EXIT_SUCCESS;
  }
  catch(const std::runtime_error &e)
  {
    clog << program_name << ": error processing command-line options: "
         << e.what() << endl;
    return EXIT_USAGE_ERROR;
  }

  if(domain_label == origin_label) {
    clog << program_name << ": domain-label and origin-label must be different"
         << std::endl;
    return EXIT_USAGE_ERROR;
  }

  // find internal labels that are different from domain and origin
  int16_t front_label = yl::DEFAULT_FRONT_LABEL;
  while(front_label == domain_label || front_label == origin_label)
    front_label--;

  int16_t done_label = yl::DEFAULT_DONE_LABEL;
  while(done_label == domain_label || done_label == origin_label
        || done_label == front_label)
    done_label--;

  const bool using_speed_map = !speed_map_reader.fileName().empty();
  if(opposite_speed && !using_speed_map) {
    clog << program_name << ": cannot have --opposite-speed without a"
      " speed map" << std::endl;
    return EXIT_USAGE_ERROR;
  }


  if(verbose != 0) clog << program_name << ": reading field..." << endl;
  AimsData<float> upwind_field;
  upwind_field_reader.setAllocatorContext(
    AllocatorContext(AllocatorStrategy::InternalModif)); // for border
  bool success = upwind_field_reader.read(upwind_field, 1 /* border */);
  if(!success) {
    clog << program_name << ": error reading file '"
         << upwind_field_reader.fileName()
         << "'specified as --upwind-field, aborting" << endl;
    return EXIT_FAILURE;
  }
  upwind_field.fillBorder(std::numeric_limits<float>::quiet_NaN());

  if(invert_field) {
    int x, y, z;
    ForEach3d(upwind_field, x, y, z)
      upwind_field(x, y, z) = -upwind_field(x, y, z);
  }

  AimsData<int16_t> domain_volume;
  if(verbose != 0) clog << program_name << ": reading domain volume..." << endl;
  domain_reader.setAllocatorContext(
    AllocatorContext(AllocatorStrategy::InternalModif));
  if(!domain_reader.read(domain_volume, 1 /* border */))
  {
    clog << program_name << ": cannot read domain volume" << endl;
    return EXIT_FAILURE;
  }
  domain_volume.fillBorder(done_label);

  carto::VolumeRef<float> speed_map;
  if(using_speed_map) {
    if(verbose != 0) clog << program_name << ": reading speed map..." << endl;
    speed_map_reader.setAllocatorContext(
      AllocatorContext(AllocatorStrategy::ReadOnly));
    bool success = speed_map_reader.read(speed_map);
    if(!success) {
      clog << program_name << ": error reading file '"
           << speed_map_reader.fileName()
           << "'specified as --speed-map, aborting" << endl;
      return EXIT_FAILURE;
    }
    if(opposite_speed) {
      for(int z = 0; z < speed_map->getSizeZ(); ++z)
        for(int y = 0; y < speed_map->getSizeY(); ++y)
          for(int x = 0; x < speed_map->getSizeX(); ++x)
            speed_map(x, y, z) = -speed_map(x, y, z);
    }
  }

  if(verbose != 0) clog << program_name << ": computing upwind distance..." << endl;
  VolumeRef<float> result_distance;
  if(using_speed_map) {
#if 0
    result_distance =
      yl::upwind_distance_with_speed(upwind_field, speed_map,
                                     domain_volume, domain_label, origin_label);
#else
    std::clog << "Not implemented yet." << std::endl;
    return EXIT_FAILURE;
#endif
  } else {
    result_distance =
      yl::upwind_distance(upwind_field.volume(),
                          domain_volume.volume(), domain_label, origin_label);
  }


  if(verbose != 0) clog << program_name << ": writing output..." << endl;
  success = output_writer.write(result_distance);
  if(!success) {
    clog << program_name << ": cannot write output volume" << endl;
  }

  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
