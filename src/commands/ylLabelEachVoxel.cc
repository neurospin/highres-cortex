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

#include <cartobase/config/verbose.h>
#include <cartodata/volume/volume.h>
#include <aims/getopt/getopt2.h>
#include <aims/io/reader.h>
#include <aims/io/writer.h>

using std::clog;
using std::endl;
using carto::VolumeRef;
using carto::verbose;
using soma::AllocatorStrategy;
using soma::AllocatorContext;

namespace
{
const int EXIT_USAGE_ERROR = 2;
std::string program_name;
}

int main(const int argc, const char **argv)
{
  // Initialize command-line option parsing
  aims::Reader<VolumeRef<int16_t> > input_reader;
  aims::Writer<VolumeRef<int32_t> > output_writer;
  int32_t first_label = 1;

  program_name = argv[0];
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
    clog << program_name << ": error processing command-line options: "
         << e.what() << endl;
    return EXIT_USAGE_ERROR;
  }

  VolumeRef<int16_t> input_mask;
  input_reader.setAllocatorContext(
    AllocatorContext(AllocatorStrategy::ReadOnly));
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
    if(input_mask(x, y, z) != 0) {
      labelled_volume(x, y, z) = next_label;
      next_label++;
    }
  }

  if(verbose != 0) {
    clog << program_name << ": assigned labels between " << first_label
         << " and " << next_label - 1 << "." << endl;
  }

  bool success = output_writer.write(labelled_volume);
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
