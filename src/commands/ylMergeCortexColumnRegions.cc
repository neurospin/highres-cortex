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


namespace {

template <class VolumeType>
float get_smallest_voxel_spacing(const VolumeType& volume)
{
  const carto::Object& voxel_size = volume.header().getProperty("voxel_size");
  assert(voxel_size->isArray());
  float voxel_sizes[] = {voxel_size->getArrayItem(0)->value<float>(),
                         voxel_size->getArrayItem(1)->value<float>(),
                         voxel_size->getArrayItem(2)->value<float>()};
  std::sort(voxel_sizes, voxel_sizes + 3);
  return voxel_sizes[0];
}

} // end of anonymous namespace


int main(const int argc, const char **argv)
{
  typedef yl::CortexColumnRegionQuality QualityCriterion;

  // Initialize command-line option parsing
  aims::Reader<VolumeRef<int32_t> > input_reader;
  aims::Reader<VolumeRef<float> > CSF_projections_reader;
  aims::Reader<VolumeRef<float> > white_projections_reader;
  float goal_diametre = QualityCriterion::default_goal_diametre();
  float max_thickness = QualityCriterion::default_max_thickness();
  aims::Writer<VolumeRef<int32_t> > output_writer;
  aims::AimsApplication app(argc, argv,
    "Aggregate oversegmented cortex column regions");
  app.addOption(input_reader, "--input", "input label volume");
  app.addOption(CSF_projections_reader, "--proj-csf",
                "projected coordinates of the CSF surfate");
  app.addOption(white_projections_reader, "--proj-white",
                "projected coordinates of the white surfate");
  {
    std::ostringstream help_str;
    help_str << "goal region diametre (millimetres) [default: "
             << goal_diametre << "]";
    app.addOption(goal_diametre, "--goal-diametre",
                  help_str.str(), true);
  }
  {
    std::ostringstream help_str;
    help_str << "maximum cortex thickness (millimetres) [default: "
             << max_thickness << "]";
    app.addOption(max_thickness, "--max-thickness",
                  help_str.str(), true);
  }
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

  VolumeRef<float> CSF_projections;
  CSF_projections_reader.read(CSF_projections);
  VolumeRef<float> white_projections;
  white_projections_reader.read(white_projections);

  QualityCriterion quality_criterion(CSF_projections,
                                     white_projections);
  quality_criterion.setShapeParametres(goal_diametre, max_thickness);

  yl::IterativeRegionMerger<int32_t, QualityCriterion>
    region_merger(input_regions, quality_criterion, verbose);
  const std::size_t max_region_size
    = std::ceil(2 * max_thickness / get_smallest_voxel_spacing(CSF_projections));
  region_merger.set_max_region_size(max_region_size);

  region_merger.merge_worst_regions_iteratively();

  VolumeRef<int32_t> output_volume = region_merger.volume();

  bool success = output_writer.write(output_volume);
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
