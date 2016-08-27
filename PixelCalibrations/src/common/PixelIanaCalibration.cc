#include "CalibFormats/SiPixelObjects/interface/PixelCalibConfiguration.h"
#include "CalibFormats/SiPixelObjects/interface/PixelDACNames.h"
#include "PixelCalibrations/include/PixelIanaCalibration.h"

//#include <toolbox/convertstring.h>

using namespace pos;
using namespace std;

PixelIanaCalibration::PixelIanaCalibration(const PixelSupervisorConfiguration & tempConfiguration, SOAPCommander* mySOAPCmdr)
  : PixelCalibrationBase(tempConfiguration, *mySOAPCmdr)
{
  std::cout << "Greetings from the PixelIanaCalibration copy constructor." << std::endl;
}

void PixelIanaCalibration::beginCalibration() {
  PixelCalibConfiguration* tempCalibObject = dynamic_cast<PixelCalibConfiguration*>(theCalibObject_);
  assert(tempCalibObject != 0);
  
}

bool PixelIanaCalibration::execute() {
  PixelCalibConfiguration* tempCalibObject = dynamic_cast<PixelCalibConfiguration*>(theCalibObject_);
  assert(tempCalibObject != 0);

  const bool firstOfPattern = event_ % tempCalibObject->nTriggersPerPattern() == 0;
  const unsigned state = event_/(tempCalibObject->nTriggersPerPattern());
  reportProgress(0.05);

  // Configure all TBMs and ROCs according to the PixelCalibConfiguration settings, but only when it's time for a new configuration.
  if (firstOfPattern) {
    //if (ToggleChannels) commandToAllFEDCrates("ToggleChannels");
    commandToAllFECCrates("CalibRunning");
  }

  const vector<PixelROCName>& rocs=tempCalibObject->rocList();
 
  for (unsigned int i=0;i<rocs.size();i++){

    //cout << "Selected ROC: " << rocs[i].rocname() << endl
    // << "Will set Vana = " << vana << endl;
    //setDAC(rocs[i], pos::k_DACAddress_Vana, vana);
    setDAC(rocs[i], pos::k_DACAddress_Readback, 12);

  }

    // Send trigger to all TBMs and ROCs.
    //for (int itrig = 0; itrig < 32; ++itrig) {
     sendTTCCalSync();
     usleep(1000);
    //}

    // Read out data from each FED.
    Attribute_Vector parametersToFED(2);
    parametersToFED[0].name_ = "WhatToDo"; parametersToFED[0].value_ = "RetrieveData";
    parametersToFED[1].name_ = "StateNum"; parametersToFED[1].value_ = itoa(state);
    commandToAllFEDCrates("FEDCalibrations", parametersToFED);
    
  return event_ + 1 < tempCalibObject->nTriggersTotal();

}

void PixelIanaCalibration::endCalibration() {
  PixelCalibConfiguration* tempCalibObject = dynamic_cast<PixelCalibConfiguration*>(theCalibObject_);
  assert(tempCalibObject != 0);
  assert(event_ == tempCalibObject->nTriggersTotal());
	
  Attribute_Vector parametersToFED(2);
  parametersToFED[0].name_ = "WhatToDo"; parametersToFED[0].value_ = "Analyze";
  parametersToFED[1].name_ = "StateNum"; parametersToFED[1].value_ = "0";
  commandToAllFEDCrates("FEDCalibrations", parametersToFED);

}

std::vector<std::string> PixelIanaCalibration::calibrated() {
  std::vector<std::string> tmp;
  return tmp;
}
