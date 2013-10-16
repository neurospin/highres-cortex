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
  aims::Reader<VolumeRef<float> > distCSF_reader;
  aims::Reader<VolumeRef<float> > distwhite_reader;
  float border_prox_weight = QualityCriterion::default_border_prox_weight();
  float compacity_weight = QualityCriterion::default_compacity_weight();
  float size_weight = QualityCriterion::default_size_weight();
  aims::Writer<carto::Volume<float> > output_writer;
  aims::AimsApplication app(argc, argv,
    "TODO");
  app.addOption(input_reader, "--input", "input label volume");
  app.addOption(distCSF_reader, "--dist-csf",
                "distance map to CSF");
  app.addOption(distwhite_reader, "--dist-white",
                "distance map to white matter");
  {
    std::ostringstream help_str;
    help_str << "weight of the border proximity criterion [default: "
             << border_prox_weight << "]";
    app.addOption(border_prox_weight, "--border-prox-weight",
                  help_str.str(), true);
  }
  {
    std::ostringstream help_str;
    help_str << "weight of the compacity criterion [default: "
             << compacity_weight << "]";
    app.addOption(compacity_weight, "--compacity-weight",
                  help_str.str(), true);
  }
  {
    std::ostringstream help_str;
    help_str << "weight of the region size criterion [default: "
             << size_weight << "]";
    app.addOption(size_weight, "--size-weight",
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

  VolumeRef<float> distCSF;
  distCSF_reader.read(distCSF);
  VolumeRef<float> distwhite;
  distwhite_reader.read(distwhite);

  const int size_x = input_regions.getSizeX();
  const int size_y = input_regions.getSizeY();
  const int size_z = input_regions.getSizeZ();
  carto::Volume<float> output_volume(size_x, size_y, size_z);
  output_volume.header() = input_regions->header();
  output_volume.fill(0);

  QualityCriterion quality_criterion(distCSF, distwhite);
  quality_criterion.setBorderProxWeight(border_prox_weight);
  quality_criterion.setCompacityWeight(compacity_weight);
  quality_criterion.setSizeWeight(size_weight);

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
