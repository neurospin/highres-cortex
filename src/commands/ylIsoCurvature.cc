/*
Copyright Forschungszentrum Jülich GmbH (2016).

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
#include <cmath>

#include <boost/format.hpp>

#include <cartobase/config/verbose.h>
#include <cartodata/volume/volume.h>
#include <aims/getopt/getopt2.h>
#include <aims/io/reader.h>
#include <aims/io/writer.h>


using std::clog;
using std::endl;
using std::flush;
using carto::VolumeRef;
using carto::verbose;
using soma::AllocatorStrategy;
using soma::AllocatorContext;


// Anonymous namespace for file-local symbols
namespace
{
const int EXIT_USAGE_ERROR = 2;
std::string program_name;

// The sum of principal curvatures is computed according to the method
// described in this paper:
//
// Jean-Philippe Thirion, and Alexis Gourdon. ‘Computing the Differential
// Characteristics of Isointensity Surfaces’. Computer Vision and Image
// Understanding 61, no. 2 (March 1995): 190–202. doi:10.1006/cviu.1995.1015.
template<typename Real>
VolumeRef<Real> sum_curvatures(const VolumeRef<Real>& input)
{
  const int size_x = input->getSizeX();
  const int size_y = input->getSizeY();
  const int size_z = input->getSizeZ();

  VolumeRef<Real> result(size_x, size_y, size_z);
  result->copyHeaderFrom(input->header());
  result->fill(std::numeric_limits<Real>::quiet_NaN());

  /* fact[xyz] is the multiplicative factor of the finite difference to obtain
     the derivative, taking the voxel size into account
     (1 / (2 * voxel spacing in direction [xyz])) */
  const std::vector<float> voxsize = input->getVoxelSize();
  const Real factx = 0.5f / voxsize[0];
  const Real facty = 0.5f / voxsize[1];
  const Real factz = 0.5f / voxsize[2];

  int slices_done = 0;
  #pragma omp parallel for schedule(dynamic)
  for(int z = 1 ; z < size_z - 1 ; ++z) {
    for(int y = 1 ; y < size_y - 1 ; ++y)
      for(int x = 1 ; x < size_x - 1 ; ++x) {
        const Real fx = factx * (input(x + 1, y, z) - input(x - 1, y, z));
        const Real fy = facty * (input(x, y + 1, z) - input(x, y - 1, z));
        const Real fz = factz * (input(x, y, z + 1) - input(x, y, z - 1));
        const Real norm2 = fx * fx + fy * fy + fz * fz;

        if(norm2 != 0) {
          /* Calculate the second derivative using only first-order
             neighbours. This allows the kernel to fit in a 26-neighbourhood
             (3x3x3) */
          const Real fxx = 4 * factx * factx
            * (input(x + 1, y, z) - 2 * input(x, y, z) + input(x - 1, y, z));
          const Real fyy = 4 * facty * facty
            * (input(x, y + 1, z) - 2 * input(x, y, z) + input(x, y - 1, z));
          const Real fzz = 4 * factz * factz
            * (input(x, y, z + 1) - 2 * input(x, y, z) + input(x, y, z - 1));
          const Real fxy = factx * facty
            * (input(x + 1, y + 1, z) - input(x - 1, y + 1, z)
               - input(x + 1, y - 1, z) + input(x - 1, y - 1, z));
          const Real fxz = factx * factz
            * (input(x + 1, y, z + 1) - input(x - 1, y, z + 1)
               - input(x + 1, y, z - 1) + input(x - 1, y, z - 1));
          const Real fyz = facty * factz
            * (input(x, y + 1, z + 1) - input(x, y - 1, z + 1)
               - input(x, y + 1, z - 1) + input(x, y - 1, z - 1));

          result(x, y, z) = (fx * fx * (fyy + fzz) - 2 * fy * fz * fyz +
                             fy * fy * (fxx + fzz) - 2 * fx * fz * fxz +
                             fz * fz * (fxx + fyy) - 2 * fx * fy * fxy)
            / (norm2 * std::sqrt(norm2));
        }
      }

    #pragma omp atomic
    ++slices_done;

    if(verbose) {
      #pragma omp critical(print_stderr)
      clog << "\r  " << slices_done << " / "
           << size_z - 2 << " slices processed. "
           << flush;
    }
  }

  clog << endl;

  return result;
}
} // end of anonymous namespace

int main(const int argc, const char **argv)
{
  typedef float Real;
  // Initialize command-line option parsing
  aims::Reader<VolumeRef<Real> > input_reader;
  aims::Writer<VolumeRef<Real> > output_writer;
  std::string mode = "mean";

  program_name = argv[0];
  aims::AimsApplication app(argc, argv,
                            "Compute the curvature of isosurfaces");
  app.addOption(input_reader, "--input",
                "input image volume (scalar field)");
  app.addOption(mode, "--mode",
                "type of curvature to compute"
                " (mean, sum, geom, pri1, pri2)");
  app.addOption(output_writer, "--output",
                "output image volume containing the curvature of the"
                " isosurfaces of the input field");
  app.alias("-i", "--input");
  app.alias("-m", "--mode");
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
    clog << program_name << ": error processing command-line options: "
         << e.what() << endl;
    return EXIT_USAGE_ERROR;
  }

  if(mode != "sum") {
    clog << "Mode " << mode << " is not implemented (yet)." << endl;
    return EXIT_USAGE_ERROR;
  }

  if(verbose) clog << program_name << ": reading input..." << endl;
  input_reader.setAllocatorContext(
    AllocatorContext(AllocatorStrategy::ReadOnly));
  VolumeRef<Real> input;
  bool success = input_reader.read(input);
  if(!success) {
    clog << program_name << ": error reading file '"
         << input_reader.fileName()
         << "'specified as --input, aborting" << endl;
    return EXIT_FAILURE;
  }

  // TODO implement other curvatures (mean, geom, pri1, pri2)
  VolumeRef<Real> result = sum_curvatures(input);


  if(verbose) clog << program_name << ": writing output..." << endl;
  {
    bool write_success = output_writer.write(result);
    if(!write_success) {
      clog << program_name << ": cannot write output volume" << endl;
    }
    success = success && write_success;
  }

  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
