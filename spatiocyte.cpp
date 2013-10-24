#include <iostream> 
#include <boost/date_time/posix_time/posix_time.hpp>
#include <Compartment.hpp>
#include <Stepper.hpp>
#include <Model.hpp>

int main()
{
  const double aVoxelRadius(2.5e-9);
  const double aLength(1e-6);
  Compartment aRootComp(aVoxelRadius, aLength, aLength, aLength);
  Stepper aStepper(aRootComp);
  Model aModel(aStepper);
  aModel.initialize();
  /*
  Species A;
  aRootComp.addSpecies(A);
  aModel.addStepper(aStepper);
  VisualizationLogProcess aVisualizer;
  aVisualizer.addSpecies(A);
  */
  boost::posix_time::ptime start(
                 boost::posix_time::microsec_clock::universal_time()); 
  aModel.run(0.1);
  boost::posix_time::ptime end(
                 boost::posix_time::microsec_clock::universal_time());
  std::cout << "duration:" << end-start << std::endl;
}
