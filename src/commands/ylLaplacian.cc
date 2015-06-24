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

#include <cartobase/config/verbose.h>
#include <cartobase/type/converter.h>
#include <cartodata/volume/volume.h>
#include <aims/utility/converter_volume.h>
#include <aims/getopt/getopt2.h>
#include <aims/io/reader.h>
#include <aims/io/writer.h>

#include <highres-cortex/cortex.hh>

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

VolumeRef<float>
PrepareInitialization(const VolumeRef<int16_t> classif)
{
  // Allocate the solution volume with a border
  static const int border_width = 1;
  VolumeRef<float> initialization(classif->getSizeX() + 2 * border_width,
                                  classif->getSizeY() + 2 * border_width,
                                  classif->getSizeZ() + 2 * border_width);
  initialization.reset( new carto::Volume<float>(
    initialization,
    Volume<float>::Position4Di(border_width, border_width, border_width, 0),
    Volume<float>::Position4Di(
      classif->getSizeX(), classif->getSizeY(), classif->getSizeZ(), 1)));

  // Initialize the output volume to a reasonable value (0 in CSF, 1 in white
  // matter, 0.5 in the cortex)
  {
    carto::RescalerInfo rescaler_info;
    rescaler_info.vmin = CSF_LABEL;  // 0
    rescaler_info.vmax = WHITE_LABEL;  // 200
    rescaler_info.omin = 0;
    rescaler_info.omax = 1;

    carto::Converter<VolumeRef<int16_t>, VolumeRef<float> >
      converter(true, rescaler_info);
    converter.convert(classif, initialization);
  }
  initialization->copyHeaderFrom(classif->header());

  return initialization;
}


// This is a Successive Over-Relaxation method using a Gauss-Seidel iteration
// with the largest time step ensuring stability, and even-odd traversal with
// Chebyshev acceleration. See W.H. Press, W.H Teukolsky, W.T. Vetterling and
// B.P. Flannery; Numerical Recipes in C++: The Art of Scientific Computing;
// (Cambridge University Press, 2002); pages 868--871.
void
SolveLaplace_SOR(VolumeRef<float> solution,
                 const VolumeRef<int16_t> classif,
                 const float precision,
                 const float typical_cortical_thickness)
{
  assert(precision > 0);
  assert(typical_cortical_thickness > 0);

  const int size_x = classif->getSizeX();
  const int size_y = classif->getSizeY();
  const int size_z = classif->getSizeZ();

  const std::vector<float> voxsize = classif->getVoxelSize();
  const float voxsize_x = voxsize[0];
  const float voxsize_y = voxsize[1];
  const float voxsize_z = voxsize[2];

  const float factx = 0.5f / (1 + square(voxsize_x / voxsize_y) +
                              square(voxsize_x / voxsize_z));
  const float facty = 0.5f / (square(voxsize_y / voxsize_x) + 1 +
                              square(voxsize_y / voxsize_z));
  const float factz = 0.5f / (square(voxsize_z / voxsize_x) +
                              square(voxsize_z / voxsize_y) + 1);

  static const float pi = 3.14159265358979323846f;
  const float cortical_voxels =
    typical_cortical_thickness * 3 / (voxsize_x + voxsize_y + voxsize_z);
  const float estimated_spectral_radius =
    (std::cos(pi / cortical_voxels) + 2) / 3;
  const float omega_update_factor = 0.25f * square(estimated_spectral_radius);

  float max_residual;
  float omega = 1.f;
  unsigned int iter = 0;
  do {
    max_residual = 0;
    for(int even_odd = 0; even_odd <= 1; ++even_odd) {
      for(int z = 0 ; z < size_z ; ++z)
        for(int y = 0 ; y < size_y ; ++y)
          for(int x = (y + z + even_odd) % 2 ; x < size_x ; x += 2) {
            if(classif(x, y, z) == yl::CORTEX_LABEL) {
              const float old_value = solution(x, y, z);
              const float gauss_seidel_new_value =
                factx * (solution(x - 1, y, z) + solution(x + 1, y, z)) +
                facty * (solution(x, y - 1, z) + solution(x, y + 1, z)) +
                factz * (solution(x, y, z - 1) + solution(x, y, z + 1));
              const float new_value = omega * gauss_seidel_new_value
                + (1.f - omega) * old_value;
              solution(x, y, z) = new_value;
              const float abs_residual = std::abs(gauss_seidel_new_value
                                                  - old_value);
              if(abs_residual > max_residual) {
                max_residual = abs_residual;
              }
              if(false) // debug
                clog << "change at (" << x << ", " << y << ", " << z
                     << "): " << old_value << " -> " << new_value << endl;
            }
          }

      // Chebyshev acceleration
      if(iter == 0 && even_odd == 0)
        omega = 1 / (1 - 0.5f * square(estimated_spectral_radius));
      else
        omega = 1 / (1 - omega * omega_update_factor);
      assert(omega > 0);
      assert(omega < 2);
    }
    ++iter;
    clog << "\riteration " << iter << ", max residual = " << max_residual
         << "   " << flush;
  } while(max_residual >= precision);
  clog << endl;
}

}

int main(const int argc, const char **argv)
{
  // Initialize command-line option parsing
  aims::Reader<VolumeRef<int16_t> > classif_reader;
  float precision = 0.001f;
  float typical_cortical_thickness = 3;
  aims::Writer<VolumeRef<float> > output_writer;

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
    help_str << "the desired precision (iteration stops when the"
      " residual goes below this value) (default: " << precision << ")";
    app.addOption(precision, "--precision", help_str.str(), true);
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

  if(!(precision > 0 && precision < 1)) {
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

  VolumeRef<float> solution = PrepareInitialization(classif);
  SolveLaplace_SOR(solution, classif, precision, typical_cortical_thickness);


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
