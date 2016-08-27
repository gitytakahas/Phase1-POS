#ifndef _PixelROCDelayCalibration_h_
#define _PixelROCDelayCalibration_h_

#include "PixelCalibrations/include/PixelCalibrationBase.h"

class PixelROCDelayCalibration : public PixelCalibrationBase {
 public:
  PixelROCDelayCalibration(const PixelSupervisorConfiguration&, SOAPCommander*);

  void beginCalibration();
  virtual bool execute();
  void endCalibration();
  virtual std::vector<std::string> calibrated();

  bool ToggleChannels;
  bool CycleScopeChannels;
  bool DelayBeforeFirstTrigger;
  bool DelayEveryTrigger;
};

#endif
