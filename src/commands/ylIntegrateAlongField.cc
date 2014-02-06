#include <cstdlib>
#include <iostream>
#include <sstream>

#include <cartobase/config/verbose.h>
#include <cartodata/volume/volume.h>
#include <aims/getopt/getopt2.h>
#include <aims/io/process.h>
#include <aims/io/finder.h>
#include <aims/io/reader.h>
#include <aims/io/writer.h>

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

#include <yleprince/propagate_along_field.hh>

using std::clog;
using std::endl;
using carto::VolumeRef;
using carto::verbose;


// Anonymous namespace for file-local symbols
namespace
{

std::string program_name;

struct IntegrateAlongFieldProcess : public aims::Process
{
  IntegrateAlongFieldProcess();

  aims::Reader<VolumeRef<float> > fieldx_reader, fieldy_reader, fieldz_reader;
  aims::Reader<VolumeRef<float> > scalar_field_reader;
  int32_t target_label;
  float step;
  unsigned int max_iter;
  aims::Writer<VolumeRef<float> > output_volume_writer;
};

template<typename Tlabel>
bool doit(aims::Process& proc_base,
          const std::string &seeds_filename,
          aims::Finder &finder)
{
  IntegrateAlongFieldProcess &proc =
    dynamic_cast<IntegrateAlongFieldProcess&>(proc_base);
  if(verbose) clog << program_name << ": reading field..." << endl;
  VolumeRef<float> fieldx, fieldy, fieldz;
  proc.fieldx_reader.read(fieldx);
  proc.fieldy_reader.read(fieldy);
  proc.fieldz_reader.read(fieldz);
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

  VolumeRef<Tlabel> seeds;
  {
    if(verbose) clog << program_name << ": reading seeds..." << endl;
    aims::Reader<VolumeRef<Tlabel> > seeds_reader(seeds_filename);
    seeds_reader.read(seeds);
  }

  VolumeRef<float> scalar_field_volume;
  if(!proc.scalar_field_reader.read(scalar_field_volume)) {
    return false; // Failure
  }
  boost::shared_ptr<yl::LinearlyInterpolatedScalarField> scalar_field =
    boost::make_shared<yl::LinearlyInterpolatedScalarField>(scalar_field_volume);

  if(seeds.getSizeX() != fieldx.getSizeX() ||
     seeds.getSizeY() != fieldx.getSizeY() ||
     seeds.getSizeZ() != fieldx.getSizeZ()) {
    clog << program_name << ": the size of the seed volume does not match "
      "the size of the field" << endl;
    return false;
  }

  yl::PropagateAlongField propagator(fieldx, fieldy, fieldz);
  propagator.setVerbose(verbose);
  propagator.setStep(proc.step);
  propagator.setMaxIter(proc.max_iter);

  Tlabel target_label = static_cast<Tlabel>(proc.target_label);
  if(target_label != proc.target_label) {
    clog << program_name << ": --target-label value is out of the range "
      "permitted by the label data type " << finder.dataType() << endl;
    return false;
  }

  VolumeRef<float> output_volume =
    propagator.evolve_unit_surface_from_region(seeds, scalar_field, target_label);

  bool success = proc.output_volume_writer.write(output_volume);

  return success;
};

IntegrateAlongFieldProcess::IntegrateAlongFieldProcess()
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
  IntegrateAlongFieldProcess proc;
  std::string seeds_filename;

  program_name = argv[0];
  aims::AimsApplication app(argc, argv,
"Integrate a scalar field along trajectories defined by vector field advection\n"
"\n"
"From all voxels of the seed volume that have the target label, the field\n"
"is followed in fixed-size steps until a seed is reached, or the maximum\n"
"number of iterations is exceeded.");
  app.addOption(seeds_filename, "--seeds",
                "volume of labels (either S16 or S32):\n"
                "- positive labels are seeds,\n"
                "- zero is the region of integration,\n"
                "- negative labels are forbidden regions.");
  app.addOption(proc.fieldx_reader, "--fieldx", "x component of vector field");
  app.addOption(proc.fieldy_reader, "--fieldy", "y component of vector field");
  app.addOption(proc.fieldz_reader, "--fieldz", "z component of vector field");
  app.addOption(proc.scalar_field_reader, "--field",
                "scalar field to be integrated");
  app.addOption(proc.target_label, "--target-label",
                "voxels having this label are used as starting points "
                "[default: 0]", true);
  app.addOption(proc.output_volume_writer, "--output",
                "output volume containing the integral results");
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

  bool success = proc.execute(seeds_filename);
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
