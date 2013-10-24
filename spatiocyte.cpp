#include <iostream> 
#include <boost/date_time/posix_time/posix_time.hpp>
#include <Stepper.hpp>
#include <Model.hpp>

int main ()
{
  Stepper aStepper;
  Model aModel(aStepper);
  aModel.initialize();
  /*
  CompartmentProcess aRootComp;
  aModel.addCompartment(aRootComp);
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
