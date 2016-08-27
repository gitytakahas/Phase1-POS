#ifndef _PixelFEDIanaCalibration_h_
#define _PixelFEDIanaCalibration_h_

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

class PixelFEDIanaCalibration: public PixelFEDCalibrationBase {
 public:
  PixelFEDIanaCalibration(const PixelFEDSupervisorConfiguration&, SOAPCommander*);

  virtual void initializeFED();
  virtual xoap::MessageReference beginCalibration(xoap::MessageReference msg);
  virtual xoap::MessageReference execute(xoap::MessageReference msg);
  virtual xoap::MessageReference endCalibration(xoap::MessageReference msg);

 private:
  void RetrieveData(unsigned int state);
  void Analyze();
  void CloseRootf();
  void BookEm(const TString& path);
  void ReadIana(void);

   struct branch{
    float pass;
    char rocName[38];
  };

  struct branch_sum{
    float deltaVana;
    float newVana;
    float newIana;
    float maxIana;
    float fitChisquare;
    char rocName[38];
  };
   
  std::map<pos::PixelModuleName,pos::PixelDACSettings*> dacsettings_;  
  TFile* rootf;
  bool inject_;
  int CurrentVana_;
  int lastVana;
  std::vector<std::string> dacsToScan;  
  //std::map<int,std::map<int,std::vector<TH1F*> > > scansROCs;
  std::map<int,std::map<int,std::map<int,std::map<int,int> > > > values;
  std::map<int,std::map<int,std::map<int,std::vector<int> > > > last_dacs;

};

#endif
