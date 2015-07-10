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
#include <iostream>
#include <sstream>
#include <cmath>

#include <boost/format.hpp>

#include <cartobase/config/verbose.h>
#include <cartodata/volume/volume.h>
#include <aims/getopt/getopt2.h>
#include <aims/io/reader.h>
#include <aims/io/writer.h>

#include <highres-cortex/laplace_solver.hh>

using std::clog;
using std::endl;
using std::flush;
using carto::VolumeRef;
using carto::Volume;
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
  typedef float Real;
  // Initialize command-line option parsing
  aims::Reader<VolumeRef<int16_t> > classif_reader;
  Real relative_precision = 0.001f;
  float typical_cortical_thickness = 3;
  aims::Writer<VolumeRef<Real> > output_writer;

  program_name = argv[0];
  aims::AimsApplication app(argc, argv,
                            "Solve the Laplacian model in the cortex");
  app.addOption(classif_reader, "--classif",
                "classification image of the cortex (100 inside, 0 in CSF, "
                "200 in white matter)");
  app.addOption(output_writer, "--output",
                "output pseudo-temperature field (from 0 in CSF "
                "to 1 in the white matter)");
  {
    std::ostringstream help_str;
    help_str << "desired maximum relative error in first-order"
      " finite differences (default: " << relative_precision << ")";
    app.addOption(relative_precision, "--precision", help_str.str(), true);
  }
  {
    std::ostringstream help_str;
    help_str << "typical thickness of the cortex (mm), used for accelerating"
      " convergence (default: " << typical_cortical_thickness << "mm)";
    app.addOption(typical_cortical_thickness, "--typical-cortical-thickness",
                  help_str.str(), true);
  }
  app.alias("-i", "--classif");
  app.alias("--input", "--classif");
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

  if(!(relative_precision > 0 && relative_precision < 1)) {
    clog << program_name << ": precision must be strictly between 0 and 1"
         << endl;
    return EXIT_USAGE_ERROR;
  }
  if(!(typical_cortical_thickness > 0)) {
    clog << program_name << ": typical cortical thickness must be"
      " strictly positive" << endl;
    return EXIT_USAGE_ERROR;
  }

  if(verbose) clog << program_name << ": reading classif..." << endl;
  classif_reader.setAllocatorContext(
    AllocatorContext(AllocatorStrategy::ReadOnly));
  VolumeRef<int16_t> classif;
  bool success = classif_reader.read(classif);
  if(!success) {
    clog << program_name << ": error reading file '"
         << classif_reader.fileName()
         << "'specified as --classif, aborting" << endl;
    return EXIT_FAILURE;
  }

  const std::vector<float> voxsize = classif->getVoxelSize();
  const float min_voxsize = *std::min_element(voxsize.begin(),
                                              voxsize.begin() + 3);
  const Real absolute_precision = relative_precision * min_voxsize
    / (2 * typical_cortical_thickness);

  if(absolute_precision < yl::LaplaceSolver<Real>::best_precision()) {
    clog << program_name
         << boost::format(": warning: requested precision (relative %1$.1e, "
                          "absolute %2$.1e) cannot be guaranteed.")
      % relative_precision % absolute_precision << std::endl;
  }

  yl::LaplaceSolver<Real> solver(classif);

  solver.set_verbosity(verbose);
  solver.initialize_solution();
  solver.SOR(absolute_precision, typical_cortical_thickness);
  solver.clamp_to_range(0, 1);
  solver.eliminate_extrema();
  VolumeRef<Real> solution = solver.solution();


  if(verbose) clog << program_name << ": writing output..." << endl;
  {
    bool write_success = output_writer.write(solution);
    if(!write_success) {
      clog << program_name << ": cannot write output volume" << endl;
    }
    success = success && write_success;
  }

  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
