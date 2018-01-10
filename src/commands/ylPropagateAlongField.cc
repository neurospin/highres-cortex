/*
Copyright Forschungszentrum Jülich GmbH (2017).
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

#include <cartobase/config/verbose.h>
#include <cartodata/volume/volume.h>
#include <aims/getopt/getopt2.h>
#include <aims/io/process.h>
#include <aims/io/finder.h>
#include <aims/io/reader.h>
#include <aims/io/writer.h>

#include <highres-cortex/propagate_along_field.hh>
#include <highres-cortex/field.hh>

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

struct PropagateAlongFieldProcess : public aims::Process
{
  PropagateAlongFieldProcess();

  aims::Reader<VolumeRef<float> > fieldx_reader, fieldy_reader, fieldz_reader;
  aims::Reader<VolumeRef<float> > grad_field_reader;
  int32_t target_label;
  float step;
  unsigned int max_iter;
  std::string output_regions_filename;
  aims::Writer<VolumeRef<float> > output_points_writer;
};

template<typename Tlabel>
bool doit(aims::Process& proc_base,
          const std::string &seeds_filename,
          aims::Finder &finder)
{
  PropagateAlongFieldProcess &proc =
    dynamic_cast<PropagateAlongFieldProcess&>(proc_base);

  boost::shared_ptr<yl::VectorField3d> vector_field;
  bool success = true;

  if(!proc.grad_field_reader.fileName().empty()) {
    // --grad-field provided
    if(verbose) clog << program_name << ": reading field..." << endl;
    proc.grad_field_reader.setAllocatorContext(
      AllocatorContext(AllocatorStrategy::ReadOnly));
    VolumeRef<float> grad_field;
    success = proc.grad_field_reader.read(grad_field);
    if(!success) {
      clog << program_name << ": error reading file '"
           << proc.fieldx_reader.fileName()
           << "'specified as --grad-field, aborting" << endl;
      return false;
    }
    vector_field = boost::make_shared<yl::LinearlyInterpolatedScalarFieldGradient>
      (grad_field);
  } else {
    // --fieldx, --fieldy, --fieldz provided
    if(verbose) clog << program_name << ": reading field..." << endl;
    VolumeRef<float> fieldx, fieldy, fieldz;
    proc.fieldx_reader.setAllocatorContext(
      AllocatorContext(AllocatorStrategy::ReadOnly));
    success = proc.fieldx_reader.read(fieldx);
    if(!success) {
      clog << program_name << ": error reading file '"
           << proc.fieldx_reader.fileName() << "'specified as --fieldx, aborting"
           << endl;
      return false;
    }
    proc.fieldy_reader.setAllocatorContext(
      AllocatorContext(AllocatorStrategy::ReadOnly));
    success = proc.fieldy_reader.read(fieldy);
    if(!success) {
      clog << program_name << ": error reading file '"
           << proc.fieldy_reader.fileName() << "'specified as --fieldy, aborting"
           << endl;
      return false;
    }
    proc.fieldz_reader.setAllocatorContext(
      AllocatorContext(AllocatorStrategy::ReadOnly));
    success = proc.fieldz_reader.read(fieldz);
    if(!success) {
      clog << program_name << ": error reading file '"
           << proc.fieldz_reader.fileName() << "'specified as --fieldz, aborting"
           << endl;
      return false;
    }
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
    vector_field = boost::make_shared<yl::LinearlyInterpolatedVectorField3d>
      (fieldx, fieldy, fieldz);
  }

  VolumeRef<Tlabel> seeds;
  {
    if(verbose) clog << program_name << ": reading seeds..." << endl;
    aims::Reader<VolumeRef<Tlabel> > seeds_reader(seeds_filename);
    seeds_reader.setAllocatorContext(
      AllocatorContext(AllocatorStrategy::ReadOnly));
    seeds_reader.read(seeds);
  }

  static const bool test_vector_field = false;
  if(test_vector_field) {
    const int size_x = seeds.getSizeX(),
      size_y = seeds.getSizeY(),
      size_z = seeds.getSizeZ();
    const std::vector<float> voxel_size = seeds->getVoxelSize();
    const float voxel_size_x = voxel_size[0];
    const float voxel_size_y = voxel_size[1];
    const float voxel_size_z = voxel_size[2];

    carto::VolumeRef<float> vector_field_volume(size_x, size_y, size_z, 3);
    vector_field_volume.header() = seeds.header();
    vector_field_volume.fill(0);
    Point3df pos, field_value;
    for(int z = 0; z < size_z; ++z) {
      clog << "test slice " << z << " / " << size_z << endl;
      pos[2] = z * voxel_size_z;
      for(int y = 0; y < size_y; ++y) {
        pos[1] = y * voxel_size_y;
        for(int x = 0; x < size_x; ++x) {
          if(seeds(x, y, z) >= 0) {
            pos[0] = x * voxel_size_x;
            try {
              vector_field->evaluate(pos, field_value);
            } catch(const yl::Field::UndefinedField&) {
              clog << "undefined field at " << pos << endl;
              continue;
            }
            vector_field_volume(x, y, z, 0) = field_value[0];
            vector_field_volume(x, y, z, 1) = field_value[1];
            vector_field_volume(x, y, z, 2) = field_value[2];
          }
        }
      }
    }
    aims::Writer<VolumeRef<float> > writer("vector_field.nii");
    writer.write(vector_field_volume);
  }

  yl::PropagateAlongField propagator(vector_field);
  propagator.setVerbose(verbose);
  propagator.setStep(proc.step);
  propagator.setMaxIter(proc.max_iter);

  Tlabel target_label = static_cast<Tlabel>(proc.target_label);
  if(target_label != proc.target_label) {
    clog << program_name << ": --target-label value is out of the range "
      "permitted by the label data type " << finder.dataType() << endl;
    return false;
  }

  // Initialize with null pointers
  VolumeRef<Tlabel> regions(static_cast<carto::Volume<Tlabel>*>(0));
  VolumeRef<float> points(static_cast<carto::Volume<float>*>(0));

  if(!proc.output_points_writer.fileName().empty()) {
    std::pair<VolumeRef<Tlabel>, VolumeRef<float> > outputs =
      propagator.propagate_regions_keeping_dests(seeds, target_label);
    regions = outputs.first;
    points = outputs.second;
  } else {
    regions = propagator.propagate_regions(seeds, target_label);
  }

  if(!regions.isNull()) {
    if(verbose) clog << program_name << ": writing output regions "
                     << proc.output_regions_filename << " as Volume of "
                     << finder.dataType() << "..." << endl;
    aims::Writer<VolumeRef<Tlabel> > writer(proc.output_regions_filename);
    success = writer.write(regions) && success;
  }

  if(!points.isNull()) {
    if(verbose) clog << program_name
                     << ": writing destination points " << endl;
    success = proc.output_points_writer.write(points) && success;
  }

  return true;
}

PropagateAlongFieldProcess::PropagateAlongFieldProcess()
  : target_label(0),
    step(yl::PropagateAlongField::default_step),
    max_iter(yl::PropagateAlongField::default_max_iter)
{
  registerProcessType("Volume", "S16", doit<int16_t>);
  registerProcessType("Volume", "S32", doit<int32_t>);
};

} // End of anonymous namespace


int main(const int argc, const char **argv)
{
  // Initialize command-line option parsing
  PropagateAlongFieldProcess proc;
  std::string seeds_filename;

  program_name = argv[0];
  aims::AimsApplication app(argc, argv,
"Propagate regions down a vector field\n"
"\n"
"From all voxels of the seed volume that have the target label, the field\n"
"is followed in fixed-size steps until a seed is reached, or the maximum\n"
"number of iterations is exceeded.");
  app.addOption(seeds_filename, "--seeds",
                "volume of labels (either S16 or S32):\n"
                "- positive labels are seeds,\n"
                "- zero is the region of propagation,\n"
                "- negative labels are forbidden regions.");
  app.addOption(proc.grad_field_reader, "--grad-field",
                "scalar field whose gradient is to be advected along", true);
  app.addOption(proc.fieldx_reader, "--fieldx",
                "x component of vector field to advect along", true);
  app.addOption(proc.fieldy_reader, "--fieldy",
                "y component of vector field to advect along", true);
  app.addOption(proc.fieldz_reader, "--fieldz",
                "z component of vector field to advect along", true);
  app.addOption(proc.target_label, "--target-label",
                "voxels having this label are used as starting points "
                "[default: 0]", true);
  {
    std::ostringstream help_str;
    help_str << "size of the advection step (millimetres) [default: "
             << proc.step << "]";
    app.addOption(proc.step, "--step", help_str.str(), true);
  }
  {
    std::ostringstream help_str;
    help_str << "abort after so many iterations [default: "
             << proc.max_iter << "]";
    app.addOption(proc.max_iter, "--max-iter", help_str.str(), true);
  }
  app.addOption(proc.output_regions_filename, "--output",
                "output the propagated regions");
  app.addOption(proc.output_points_writer, "--dest-points",
                "output the destination points for each propagated voxel",
                true);


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

  bool fieldx_provided = !proc.fieldx_reader.fileName().empty();
  bool fieldy_provided = !proc.fieldy_reader.fileName().empty();
  bool fieldz_provided = !proc.fieldz_reader.fileName().empty();
  bool grad_field_provided = !proc.grad_field_reader.fileName().empty();

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

  bool success = proc.execute(seeds_filename);

  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
