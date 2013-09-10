#include <cstdlib>
#include <iostream>
#include <sstream>

#include <cartobase/config/verbose.h>
#include <cartodata/volume/volume.h>
#include <aims/getopt/getopt2.h>
#include <aims/io/reader.h>
#include <aims/io/writer.h>

#include <yleprince/propagate_along_field.hh>

using carto::VolumeRef;


int main(const int argc, const char **argv)
{
  // Initialize command-line option parsing
  aims::Reader<VolumeRef<int16_t> > seeds_reader;
  aims::Reader<VolumeRef<float> > fieldx_reader, fieldy_reader, fieldz_reader;
  aims::Writer<VolumeRef<int16_t> > regions_writer;
  float step = yl::PropagateAlongField::default_step;
  unsigned int max_iter = yl::PropagateAlongField::default_max_iter;
  int16_t target_label = 0;

  aims::AimsApplication app(argc, argv,
                            "Propagate regions down a vector field\n"
"\n"
"From all voxels of the seed volume that have the target label, the field\n"
"is followed in fixed-size steps until a seed is reached, or the maximum\n"
"number of iterations is exceeded.");
  app.addOption(seeds_reader, "--seeds",
                "volume of seeds: \n"
                "- positive labels are seeds, \n"
                "- zero is the region of propagation, \n"
                "- negative labels are forbidden regions.");
  app.addOption(fieldx_reader, "--fieldx", "x component of vector field");
  app.addOption(fieldy_reader, "--fieldy", "y component of vector field");
  app.addOption(fieldz_reader, "--fieldz", "z component of vector field");
  app.addOption(target_label, "--target-label",
                "only voxels with this label are used as starting points"
                "[default: 0]", true);
  {
    std::ostringstream help_str;
    help_str << "move in steps this big (millimetres) [default: "
             << step << "]";
    app.addOption(step, "--step", help_str.str(), true);
  }
  {
    std::ostringstream help_str;
    help_str << "abort after so many iterations [default: "
             << max_iter << "]";
    app.addOption(max_iter, "--max-iter", help_str.str(), true);
  }
  app.addOption(regions_writer, "--output", "output the propagated regions");


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
    std::clog << argv[0] << ": error processing command-line options: "
              << e.what() << std::endl;
    return EXIT_FAILURE;
  }

  VolumeRef<float> fieldx, fieldy, fieldz;
  fieldx_reader.read(fieldx);
  fieldy_reader.read(fieldy);
  fieldz_reader.read(fieldz);

  VolumeRef<int16_t> seeds;
  seeds_reader.read(seeds);

  if(fieldx.getSizeX() != fieldy.getSizeX() ||
     fieldx.getSizeX() != fieldz.getSizeX() ||
     fieldx.getSizeY() != fieldy.getSizeY() ||
     fieldx.getSizeY() != fieldz.getSizeY() ||
     fieldx.getSizeZ() != fieldy.getSizeZ() ||
     fieldx.getSizeZ() != fieldz.getSizeZ()) {
    std::clog << argv[0] << ": the sizes of the field volumes do not match"
              << std::endl;
    return EXIT_FAILURE;
  }

  if(seeds.getSizeX() != fieldx.getSizeX() ||
     seeds.getSizeY() != fieldx.getSizeY() ||
     seeds.getSizeZ() != fieldx.getSizeZ()) {
    std::clog << argv[0] << ": the size of the seed volume does not match "
      "the size of the field volumes" << std::endl;
    return EXIT_FAILURE;
  }

  yl::PropagateAlongField propagator(fieldx, fieldy, fieldz);
  VolumeRef<int16_t> regions = propagator.propagate_regions(seeds);

  regions_writer.write(regions);

  return EXIT_SUCCESS;
}
