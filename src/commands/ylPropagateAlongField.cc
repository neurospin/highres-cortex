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

#include <yleprince/propagate_along_field.hh>
#include <yleprince/field.hh>

using std::clog;
using std::endl;
using carto::VolumeRef;
using carto::verbose;


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

  bool fieldx_provided = !proc.fieldx_reader.fileName().empty();
  bool fieldy_provided = !proc.fieldy_reader.fileName().empty();
  bool fieldz_provided = !proc.fieldz_reader.fileName().empty();
  bool grad_field_provided = !proc.grad_field_reader.fileName().empty();

  boost::shared_ptr<yl::VectorField3d> vector_field;
  bool success = true;

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
    success = proc.grad_field_reader.read(grad_field);
    if(!success) {
      clog << program_name << ": error reading file '"
           << proc.fieldx_reader.fileName()
           << "'specified as --grad-field, aborting" << endl;
      return EXIT_FAILURE;
    }
    vector_field = boost::make_shared<yl::LinearlyInterpolatedScalarFieldGradient>
      (grad_field);
  } else {
    // --fieldx, --fieldy, --fieldz provided
    if(verbose) clog << program_name << ": reading field..." << endl;
    VolumeRef<float> fieldx, fieldy, fieldz;
    success = proc.fieldx_reader.read(fieldx);
    if(!success) {
      clog << program_name << ": error reading file '"
           << proc.fieldx_reader.fileName() << "'specified as --fieldx, aborting"
           << endl;
      return EXIT_FAILURE;
    }
    success = proc.fieldy_reader.read(fieldy);
    if(!success) {
      clog << program_name << ": error reading file '"
           << proc.fieldy_reader.fileName() << "'specified as --fieldy, aborting"
           << endl;
      return EXIT_FAILURE;
    }
    success = proc.fieldz_reader.read(fieldz);
    if(!success) {
      clog << program_name << ": error reading file '"
           << proc.fieldz_reader.fileName() << "'specified as --fieldz, aborting"
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
    vector_field = boost::make_shared<yl::LinearlyInterpolatedVectorField3d>
      (fieldx, fieldy, fieldz);
  }

  VolumeRef<Tlabel> seeds;
  {
    if(verbose) clog << program_name << ": reading seeds..." << endl;
    aims::Reader<VolumeRef<Tlabel> > seeds_reader(seeds_filename);
    seeds_reader.read(seeds);
  }

  yl::PropagateAlongField propagator(vector_field);
  propagator.setVerbose(verbose);
  propagator.setStep(proc.step);
  propagator.setMaxIter(proc.max_iter);

  Tlabel target_label = static_cast<Tlabel>(proc.target_label);
  if(target_label != proc.target_label) {
    clog << program_name << ": --target-label value is out of the range "
      "permitted by the label data type " << finder.dataType() << endl;
    return EXIT_FAILURE;
  }

  VolumeRef<Tlabel> regions;
  VolumeRef<float> points;

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

  return EXIT_SUCCESS;
};

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
                "use the gradient of this scalar field", true);
  app.addOption(proc.fieldx_reader, "--fieldx",
                "x component of vector field", true);
  app.addOption(proc.fieldy_reader, "--fieldy",
                "y component of vector field", true);
  app.addOption(proc.fieldz_reader, "--fieldz",
                "z component of vector field", true);
  app.addOption(proc.target_label, "--target-label",
                "voxels having this label are used as starting points "
                "[default: 0]", true);
  {
    std::ostringstream help_str;
    help_str << "move in steps this big (millimetres) [default: "
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

  return proc.execute(seeds_filename);
}
