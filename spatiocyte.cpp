#include <iostream> 
#include <boost/date_time/posix_time/posix_time.hpp>
#include <SpatiocyteStepper.hpp>

int main ()
{
  SpatiocyteStepper aStepper;
  aStepper.initialize();
  const unsigned steps(0.1/4.16667e-6);
  boost::posix_time::ptime start(
                 boost::posix_time::microsec_clock::universal_time()); 
  for(unsigned i(0); i !=steps; ++i)
    {
      aStepper.step();
    }
  boost::posix_time::ptime end(
                 boost::posix_time::microsec_clock::universal_time());
  std::cout << "duration:" << end-start << std::endl;
}
