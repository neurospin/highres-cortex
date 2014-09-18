#include <cstdlib>
#include <iostream>
#include <sstream>

#include <cartobase/config/verbose.h>
#include <cartodata/volume/volume.h>
#include <aims/getopt/getopt2.h>
#include <aims/io/reader.h>
#include <aims/io/writer.h>

#include <yleprince/label_volume.hh>
#include <yleprince/cortex_column_region_quality.hh>

using std::clog;
using std::endl;
using carto::VolumeRef;
using carto::verbose;


int main(const int argc, const char **argv)
{
  typedef yl::CortexColumnRegionQuality QualityCriterion;

  // Initialize command-line option parsing
  aims::Reader<VolumeRef<int32_t> > input_reader;
  aims::Reader<VolumeRef<float> > CSF_projections_reader;
  aims::Reader<VolumeRef<float> > white_projections_reader;
  float goal_diametre = QualityCriterion::default_goal_diametre();
  float max_thickness = QualityCriterion::default_max_thickness();
  aims::Writer<carto::Volume<float> > output_writer;
  aims::AimsApplication app(argc, argv,
    "TODO");
  app.addOption(input_reader, "--input", "input label volume");
  app.addOption(CSF_projections_reader, "--proj-csf",
                "projected coordinates of the CSF surfate");
  app.addOption(white_projections_reader, "--proj-white",
                "projected coordinates of the white surfate");
  {
    std::ostringstream help_str;
    help_str << "goal region diametre (in millimetres) [default: "
             << goal_diametre << "]";
    app.addOption(goal_diametre, "--goal-diametre",
                  help_str.str(), true);
  }
  {
    std::ostringstream help_str;
    help_str << "maximum cortex thickness [default: "
             << max_thickness << "]";
    app.addOption(max_thickness, "--max-thickness",
                  help_str.str(), true);
  }
  app.addOption(output_writer, "--output", "output quality map");
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
  yl::LabelVolume<int32_t> label_volume(input_regions);

  VolumeRef<float> CSF_projections;
  CSF_projections_reader.read(CSF_projections);
  VolumeRef<float> white_projections;
  white_projections_reader.read(white_projections);

  const int extension_x = input_regions.getSizeX();
  const int extension_y = input_regions.getSizeY();
  const int extension_z = input_regions.getSizeZ();
  carto::Volume<float> output_volume(extension_x, extension_y, extension_z);
  output_volume.copyHeaderFrom(input_regions->header());
  output_volume.fill(0);

  QualityCriterion quality_criterion(CSF_projections,
                                     white_projections);
  quality_criterion.setShapeParametres(goal_diametre, max_thickness);

  for(typename yl::LabelVolume<int32_t>::const_regions_iterator
        labels_it = label_volume.regions_begin(),
        labels_end = label_volume.regions_end();
      labels_it != labels_end;
      ++labels_it)
  {
    const int32_t label = *labels_it;
    const float quality = quality_criterion.evaluate(label_volume, label);
    for(yl::LabelVolume<int32_t>::const_point_iterator
          point_it = label_volume.region_begin(label),
          point_end = label_volume.region_end(label);
        point_it != point_end;
        ++point_it)
    {
      const Point3d& point = *point_it;
      const int x = point[0];
      const int y = point[1];
      const int z = point[2];
      output_volume(x, y, z) = quality;
    }
  }

  bool success = output_writer.write(output_volume);
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
