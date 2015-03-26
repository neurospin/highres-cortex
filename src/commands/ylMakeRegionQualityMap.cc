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

#include <highres-cortex/label_volume.hh>
#include <highres-cortex/cortex_column_region_quality.hh>

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
  aims::Reader<VolumeRef<int16_t> > classif_reader;
  float goal_diameter = QualityCriterion::default_goal_diameter();
  aims::Writer<carto::Volume<float> > output_writer;
  aims::AimsApplication app(argc, argv,
    "Map the quality of cortex column regions");
  app.addOption(input_reader, "--input", "input label volume");
  app.addOption(CSF_projections_reader, "--proj-csf",
                "projected coordinates of the CSF surface");
  app.addOption(white_projections_reader, "--proj-white",
                "projected coordinates of the white surface");
  app.addOption(classif_reader, "--classif",
                "grey/white/CSF classification image");
  {
    std::ostringstream help_str;
    help_str << "goal region diameter (in millimetres) [default: "
             << goal_diameter << "]";
    app.addOption(goal_diameter, "--goal-diameter",
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
  VolumeRef<int16_t> classif;
  classif_reader.read(classif);

  const int extension_x = input_regions.getSizeX();
  const int extension_y = input_regions.getSizeY();
  const int extension_z = input_regions.getSizeZ();
  carto::Volume<float> output_volume(extension_x, extension_y, extension_z);
  output_volume.copyHeaderFrom(input_regions->header());
  output_volume.fill(0);

  QualityCriterion quality_criterion(CSF_projections,
                                     white_projections,
                                     classif);
  quality_criterion.setShapeParametres(goal_diameter);

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
