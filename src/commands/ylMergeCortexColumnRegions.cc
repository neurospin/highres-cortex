#include <cstdlib>
#include <iostream>
#include <sstream>

#include <cartobase/config/verbose.h>
#include <cartodata/volume/volume.h>
#include <aims/getopt/getopt2.h>
#include <aims/io/reader.h>
#include <aims/io/writer.h>

#include <yleprince/iterative_region_merger.hh>
#include <yleprince/cortex_column_region_quality.hh>

using std::clog;
using std::endl;
using carto::VolumeRef;
using carto::verbose;


int main(const int argc, const char **argv)
{
  // Initialize command-line option parsing
  aims::Reader<VolumeRef<int32_t> > input_reader;
  aims::Reader<VolumeRef<float> > distCSF_reader;
  aims::Reader<VolumeRef<float> > distwhite_reader;
  aims::Writer<VolumeRef<int32_t> > output_writer;
  aims::AimsApplication app(argc, argv,
    "TODO");
  app.addOption(input_reader, "--input", "input label volume");
  app.addOption(distCSF_reader, "--dist-csf",
                "distance map to CSF");
  app.addOption(distwhite_reader, "--dist-white",
                "distance map to white matter");
  app.addOption(output_writer, "--output", "output label volume");
  app.alias("-i", "--input");
  app.alias("-o", "--output");


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
    clog << argv[0] << ": error processing command-line options: "
         << e.what() << endl;
    return EXIT_FAILURE;
  }

  VolumeRef<int32_t> input_regions;
  input_reader.read(input_regions);
  VolumeRef<float> distCSF;
  distCSF_reader.read(distCSF);
  VolumeRef<float> distwhite;
  distwhite_reader.read(distwhite);

  yl::IterativeRegionMerger<int32_t, yl::CortexColumnRegionQuality>
    region_merger(input_regions,
                  yl::CortexColumnRegionQuality(distCSF, distwhite),
                  verbose);

  region_merger.merge_worst_regions_iteratively();

  VolumeRef<int32_t> output_volume = region_merger.volume();

  bool success = output_writer.write(output_volume);
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
