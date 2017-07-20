/*
Copyright CEA (2014).
Copyright Université Paris XI (2014).

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
#include <iostream>
#include <sstream>

#include <boost/make_shared.hpp>

#include <soma-io/allocator/allocator.h>
#include <cartobase/config/verbose.h>
#include <cartodata/volume/volume.h>
#include <aims/getopt/getopt2.h>
#include <aims/io/reader.h>
#include <aims/io/writer.h>

#include <highres-cortex/field.hh>
#include <highres-cortex/cortex_advection.hh>

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
}


int main(const int argc, const char **argv)
{
  // Initialize command-line option parsing
  aims::Reader<VolumeRef<float> > fieldx_reader, fieldy_reader, fieldz_reader;
  aims::Reader<VolumeRef<float> > grad_field_reader;
  aims::Reader<VolumeRef<int16_t> > domain_reader;
  aims::Reader<VolumeRef<int16_t> > advection_domain_reader;
  float step = 0.03f;
  float max_advection_distance = 6.f;
  aims::Writer<VolumeRef<float> > length_output_writer;

  program_name = argv[0];
  aims::AimsApplication app(argc, argv,
"Advect a line from each voxel, keeping track of its length."
);
  app.addOption(domain_reader, "--domain",
                "mask of the calculation domain: one inside, zero outside");
  app.addOption(advection_domain_reader, "--advect-domain",
                "mask of the advection seeds domain: one inside, zero "
                "outside - default: same as domain", true);
  app.addOption(grad_field_reader, "--grad-field",
                "use the gradient of this scalar field", true);
  app.addOption(fieldx_reader, "--fieldx",
                "x component of vector field", true);
  app.addOption(fieldy_reader, "--fieldy",
                "y component of vector field", true);
  app.addOption(fieldz_reader, "--fieldz",
                "z component of vector field", true);
  app.addOption(length_output_writer, "--output-length",
                "output volume containing the advection distance");
  {
    std::ostringstream help_str;
    help_str << "move in steps this big (millimetres) [default: "
             << step << "]";
    app.addOption(step, "--step", help_str.str(), true);
  }
  {
    std::ostringstream help_str;
    help_str << "maximum advection distance (millimetres) [default: "
             << max_advection_distance << "]";
    app.addOption(max_advection_distance, "--max-dist", help_str.str(), true);
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

  boost::shared_ptr<yl::VectorField3d> advection_field;
  bool success = true;

  bool fieldx_provided = !fieldx_reader.fileName().empty();
  bool fieldy_provided = !fieldy_reader.fileName().empty();
  bool fieldz_provided = !fieldz_reader.fileName().empty();
  bool grad_field_provided = !grad_field_reader.fileName().empty();

  if(!grad_field_provided && !(fieldx_provided
                               && fieldy_provided
                               && fieldz_provided)) {
    clog << program_name
         << ": must provide either --grad-field or --field{x,y,z}"
         << endl;
    return EXIT_USAGE_ERROR;
  }

  if(grad_field_provided && (fieldx_provided
                             || fieldy_provided
                             || fieldz_provided)) {
    clog << program_name
         << ": must provide either --grad-field or --field{x,y,z}, not both"
         << endl;
    return EXIT_USAGE_ERROR;
  }

  if(grad_field_provided) {
    // --grad-field provided
    if(verbose) clog << program_name << ": reading field..." << endl;
    VolumeRef<float> grad_field;
    grad_field_reader.setAllocatorContext(
      AllocatorContext(AllocatorStrategy::ReadOnly));
    success = grad_field_reader.read(grad_field);
    if(!success) {
      clog << program_name << ": error reading file '"
           << fieldx_reader.fileName()
           << "'specified as --grad-field, aborting" << endl;
      return EXIT_FAILURE;
    }
    advection_field = boost::make_shared<yl::LinearlyInterpolatedScalarFieldGradient>
      (grad_field);
  } else {
    // --fieldx, --fieldy, --fieldz provided
    if(verbose) clog << program_name << ": reading field..." << endl;
    VolumeRef<float> fieldx, fieldy, fieldz;
    fieldx_reader.setAllocatorContext(
      AllocatorContext(AllocatorStrategy::ReadOnly));
    success = fieldx_reader.read(fieldx);
    if(!success) {
      clog << program_name << ": error reading file '"
           << fieldx_reader.fileName() << "'specified as --fieldx, aborting"
           << endl;
      return EXIT_FAILURE;
    }
    fieldy_reader.setAllocatorContext(
      AllocatorContext(AllocatorStrategy::ReadOnly));
    success = fieldy_reader.read(fieldy);
    if(!success) {
      clog << program_name << ": error reading file '"
           << fieldy_reader.fileName() << "'specified as --fieldy, aborting"
           << endl;
      return EXIT_FAILURE;
    }
    fieldz_reader.setAllocatorContext(
      AllocatorContext(AllocatorStrategy::ReadOnly));
    success = fieldz_reader.read(fieldz);
    if(!success) {
      clog << program_name << ": error reading file '"
           << fieldz_reader.fileName() << "'specified as --fieldz, aborting"
           << endl;
      return EXIT_FAILURE;
    }
    if(fieldx.getSizeX() != fieldy.getSizeX() ||
       fieldx.getSizeX() != fieldz.getSizeX() ||
       fieldx.getSizeY() != fieldy.getSizeY() ||
       fieldx.getSizeY() != fieldz.getSizeY() ||
       fieldx.getSizeZ() != fieldy.getSizeZ() ||
       fieldx.getSizeZ() != fieldz.getSizeZ()) {
      clog << program_name << ": the sizes of the field volumes do not match"
           << endl;
      return EXIT_FAILURE;
    }
    advection_field = boost::make_shared<yl::LinearlyInterpolatedVectorField3d>
      (fieldx, fieldy, fieldz);
  }

  VolumeRef<int16_t> domain_volume;
  if(verbose) clog << program_name << ": reading domain volume..." << endl;
  domain_reader.setAllocatorContext(
    AllocatorContext(AllocatorStrategy::ReadOnly));
  if(!domain_reader.read(domain_volume))
  {
    clog << program_name << ": cannot read domain volume" << endl;
    return EXIT_FAILURE;
  }

  VolumeRef<int16_t> advection_domain_volume;
  if( !advection_domain_reader.fileName().empty() )
  {
    if(verbose) clog << program_name << ": reading advection domain volume..."
      << endl;
    advection_domain_reader.setAllocatorContext(
      AllocatorContext(AllocatorStrategy::ReadOnly));
    if(!advection_domain_reader.read(advection_domain_volume))
    {
      clog << program_name << ": cannot read advecton domain volume" << endl;
      return EXIT_FAILURE;
    }
  }

  VolumeRef<float> result_distance =
    yl::advect_euclidean(*advection_field, domain_volume,
                         max_advection_distance, step, verbose,
                         advection_domain_volume);

  success = length_output_writer.write(result_distance);
  if(!success) {
    clog << program_name << ": cannot write output volume" << endl;
  }

  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
