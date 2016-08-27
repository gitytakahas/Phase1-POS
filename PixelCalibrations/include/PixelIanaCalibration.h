#ifndef _PixelIanaCalibration_h_
#define _PixelIanaCalibration_h_

#include "PixelCalibrations/include/PixelCalibrationBase.h"

class PixelIanaCalibration : public PixelCalibrationBase {
 public:
  PixelIanaCalibration(const PixelSupervisorConfiguration&, SOAPCommander*);

  void beginCalibration();
  virtual bool execute();
  void endCalibration();
  virtual std::vector<std::string> calibrated();

  /*bool ToggleChannels;
  bool CycleScopeChannels;
  bool DelayBeforeFirstTrigger;
  bool DelayEveryTrigger;*/
};

#endif
