#include <cstdlib>
#include <iostream>
#include <sstream>

#include <cartobase/config/verbose.h>
#include <cartodata/volume/volume.h>
#include <aims/getopt/getopt2.h>
#include <aims/io/reader.h>
#include <aims/io/writer.h>

#include <yleprince/field.hh>
#include <yleprince/isovolume.hh>

using std::clog;
using std::endl;
using carto::VolumeRef;
using carto::verbose;


// Anonymous namespace for file-local symbols
namespace
{
std::string program_name;
}


int main(const int argc, const char **argv)
{
  // Initialize command-line option parsing
  aims::Reader<VolumeRef<float> > fieldx_reader, fieldy_reader, fieldz_reader;
  aims::Reader<VolumeRef<float> > divergence_field_reader;
  aims::Reader<VolumeRef<int16_t> > domain_reader;
  float step = 0.03f;
  float max_advection_distance = 6.f;
  aims::Writer<VolumeRef<float> > volume_output_writer, surface_output_writer;

  program_name = argv[0];
  aims::AimsApplication app(argc, argv,
"Advect a tube from each voxel, keeping track of its volume and end surface."
);
  app.addOption(domain_reader, "--domain",
                "mask of the calculation domain: one inside, zero outside");
  app.addOption(fieldx_reader, "--fieldx", "x component of vector field");
  app.addOption(fieldy_reader, "--fieldy", "y component of vector field");
  app.addOption(fieldz_reader, "--fieldz", "z component of vector field");
  app.addOption(divergence_field_reader, "--divergence",
                "divergence of the normalized vector field");
  app.addOption(volume_output_writer, "--output-volumes",
                "output volume containing the tubes' volume");
  app.addOption(surface_output_writer, "--output-surfaces",
                "output volume containing the tubes' end surface");
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
    return EXIT_FAILURE;
  }

  if(verbose) clog << program_name << ": reading field..." << endl;
  VolumeRef<float> fieldx, fieldy, fieldz;
  fieldx_reader.read(fieldx);
  fieldy_reader.read(fieldy);
  fieldz_reader.read(fieldz);
  if(fieldx.getSizeX() != fieldy.getSizeX() ||
     fieldx.getSizeX() != fieldz.getSizeX() ||
     fieldx.getSizeY() != fieldy.getSizeY() ||
     fieldx.getSizeY() != fieldz.getSizeY() ||
     fieldx.getSizeZ() != fieldy.getSizeZ() ||
     fieldx.getSizeZ() != fieldz.getSizeZ()) {
    clog << program_name << ": the sizes of the field volumes do not match"
         << endl;
    return false;
  }

  yl::LinearlyInterpolatedVectorField3d advection_field(fieldx, fieldy, fieldz);

  VolumeRef<int16_t> domain_volume;
  if(verbose) clog << program_name << ": reading domain volume..." << endl;
  if(!domain_reader.read(domain_volume))
  {
    clog << program_name << ": cannot read domain volume" << endl;
    return EXIT_FAILURE;
  }

  VolumeRef<float> divergence_field_volume;
  if(!divergence_field_reader.read(divergence_field_volume)) {
    return false; // Failure
  }
  yl::LinearlyInterpolatedScalarField divergence_field(divergence_field_volume);

  // yl::PropagateAlongField propagator(fieldx, fieldy, fieldz);
  // propagator.setVerbose(verbose);
  // propagator.setStep(step);
  // propagator.setMaxIter(max_iter);

  std::pair<VolumeRef<float>, VolumeRef<float> > results =
    yl::advect_tubes(advection_field, divergence_field, domain_volume,
                     max_advection_distance, step, verbose);

  bool success = true;
  {
    bool write_success = volume_output_writer.write(results.first);
    if(!write_success) {
      clog << program_name << ": cannot write output volume" << endl;
    }
    success = success && write_success;
  }
  {
    bool write_success = surface_output_writer.write(results.second);
    if(!write_success) {
      clog << program_name << ": cannot write output surface" << endl;
    }
    success = success && write_success;
  }

  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
