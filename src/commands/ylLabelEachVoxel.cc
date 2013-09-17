#include <cstdlib>
#include <iostream>
#include <sstream>

#include <cartobase/config/verbose.h>
#include <cartodata/volume/volume.h>
#include <aims/getopt/getopt2.h>
#include <aims/io/reader.h>
#include <aims/io/writer.h>

using std::clog;
using std::endl;
using carto::VolumeRef;
using carto::verbose;


int main(const int argc, const char **argv)
{
  // Initialize command-line option parsing
  aims::Reader<VolumeRef<int16_t> > input_reader;
  aims::Writer<VolumeRef<int32_t> > output_writer;
  int32_t first_label = 1;
  aims::AimsApplication app(argc, argv,
    "Assign a unique label to each voxel of a mask.\n"
    "\n"
    "Non-zero voxels of the mask are each assigned a positive unique label.");
  app.addOption(input_reader, "--input", "input mask");
  app.addOption(output_writer, "--output",
                "output label volume with S32 datatype");
  {
    std::ostringstream help_str;
    help_str << "assign labels starting with this value [default: "
             << first_label << "]";
    app.addOption(first_label, "--first", help_str.str(), true);
  }
  app.alias("-i", "--input");
  app.alias("-o", "--output");
  app.alias("--first-label", "--first");


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

  VolumeRef<int16_t> input_mask;
  input_reader.read(input_mask);
  const int size_x = input_mask.getSizeX();
  const int size_y = input_mask.getSizeY();
  const int size_z = input_mask.getSizeZ();
  VolumeRef<int32_t> labelled_volume(size_x, size_y, size_z);
  labelled_volume->header() = input_mask->header();
  labelled_volume->fill(0);

  int32_t next_label = first_label;
  for(int z = 0; z < size_z; ++z)
  for(int y = 0; y < size_y; ++y)
  for(int x = 0; x < size_x; ++x)
  {
    if(input_mask(x, y, z)) {
      labelled_volume(x, y, z) = next_label;
      next_label++;
    }
  }

  if(verbose) {
    clog << argv[0] << ": assigned labels between " << first_label
         << " and " << next_label - 1 << "." << endl;
  }

  bool success = output_writer.write(labelled_volume);
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
