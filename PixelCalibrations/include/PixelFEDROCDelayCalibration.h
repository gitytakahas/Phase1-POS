#ifndef _PixelFEDROCDelayCalibration_h_
#define _PixelFEDROCDelayCalibration_h_

#include "CalibFormats/SiPixelObjects/interface/PixelROCName.h"
#include "PixelCalibrations/include/PixelFEDCalibrationBase.h"
#include "PixelUtilities/PixelFEDDataTools/include/Moments.h"
#include "PixelUtilities/PixelFEDDataTools/include/PixelScanRecord.h"

#include <cstdint>
#include <fstream>

class TFile;
class TH1F;
class TH2F;
class TH3F;

class PixelFEDROCDelayCalibration: public PixelFEDCalibrationBase {
 public:
  PixelFEDROCDelayCalibration(const PixelFEDSupervisorConfiguration&, SOAPCommander*);

  virtual void initializeFED();
  virtual xoap::MessageReference beginCalibration(xoap::MessageReference msg);
  virtual xoap::MessageReference execute(xoap::MessageReference msg);
  virtual xoap::MessageReference endCalibration(xoap::MessageReference msg);

 private:
  void RetrieveData(unsigned int state);
  void Analyze();
  void CloseRootf();
  void BookEm(const TString& path);
  void FillEm(unsigned state,int fedid, int ch, int which, float c, int roc);

  bool DumpFIFOs;
  bool ReadFifo1;
  bool ReadFifo3;
  std::vector<std::string> dacsToScan;
  unsigned lastTBMADelay;
  TFile* rootf;
  bool inject_;

  std::map<int,std::map<int,std::vector<TH2F*> > > scansROCs;
  std::map<int,std::map<std::string,std::vector<TH2F*> > > ROCsHistoSum;

};

#endif
