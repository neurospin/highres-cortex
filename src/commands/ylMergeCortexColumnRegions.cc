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

#include <highres-cortex/iterative_region_merger.hh>
#include <highres-cortex/cortex_column_region_quality.hh>

using std::clog;
using std::endl;
using carto::VolumeRef;
using carto::verbose;


namespace {

template <class VolumeType>
float get_smallest_voxel_spacing(const VolumeType& volume)
{
  std::vector<float> voxel_size = volume->getVoxelSize();
  voxel_size.resize(3);
  std::sort(voxel_size.begin(), voxel_size.end());
  return voxel_size[0];
}

} // end of anonymous namespace


int main(const int argc, const char **argv)
{
  typedef yl::CortexColumnRegionQuality QualityCriterion;

  // Initialize command-line option parsing
  aims::Reader<VolumeRef<int32_t> > input_reader;
  aims::Reader<VolumeRef<float> > CSF_projections_reader;
  aims::Reader<VolumeRef<float> > white_projections_reader;
  aims::Reader<VolumeRef<int16_t> > classif_reader;
  float goal_diameter = QualityCriterion::default_goal_diameter();
  aims::Writer<VolumeRef<int32_t> > output_writer;
  aims::AimsApplication app(argc, argv,
    "Aggregate over-segmented cortex column regions");
  app.addOption(input_reader, "--input", "input label volume");
  app.addOption(CSF_projections_reader, "--proj-csf",
                "projected coordinates of the CSF surface");
  app.addOption(white_projections_reader, "--proj-white",
                "projected coordinates of the white surface");
  app.addOption(classif_reader, "--classif",
                "grey/white/CSF classification image");
  {
    std::ostringstream help_str;
    help_str << "goal region diameter (millimetres) [default: "
             << goal_diameter << "]";
    app.addOption(goal_diameter, "--goal-diameter",
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
  VolumeRef<int16_t> classif;
  classif_reader.read(classif);

  QualityCriterion quality_criterion(CSF_projections,
                                     white_projections,
                                     classif);
  quality_criterion.setShapeParametres(goal_diameter);

  yl::IterativeRegionMerger<int32_t, QualityCriterion>
    region_merger(input_regions, quality_criterion, verbose);

  region_merger.merge_worst_regions_iteratively();

  VolumeRef<int32_t> output_volume = region_merger.volume();

  bool success = output_writer.write(output_volume);
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
