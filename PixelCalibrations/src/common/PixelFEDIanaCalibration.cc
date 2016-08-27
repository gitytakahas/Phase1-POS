#include "CalibFormats/SiPixelObjects/interface/PixelCalibConfiguration.h"
#include "CalibFormats/SiPixelObjects/interface/PixelDACNames.h"
#include "PixelCalibrations/include/PixelFEDIanaCalibration.h"
#include "PixelConfigDBInterface/include/PixelConfigInterface.h"
#include "PixelUtilities/PixelFEDDataTools/include/PixelFEDDataTypes.h"
#include "PixelUtilities/PixelFEDDataTools/include/ErrorFIFODecoder.h"
#include "PixelUtilities/PixelFEDDataTools/include/ColRowAddrDecoder.h"
#include "PixelUtilities/PixelFEDDataTools/include/DigScopeDecoder.h"
#include "PixelUtilities/PixelFEDDataTools/include/DigTransDecoder.h"
#include "PixelUtilities/PixelFEDDataTools/include/FIFO3Decoder.h"
#include "PixelUtilities/PixelRootUtilities/include/PixelRootDirectoryMaker.h"
#include "PixelUtilities/PixelFEDDataTools/include/DigFIFO1Decoder.h"
#include "PixelCalibrations/include/PixelIanaAnalysis.h"


#include "TH1F.h"
#include "TH2F.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TFile.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TTree.h"

#include <iomanip>
#include <algorithm>

using namespace pos;

///////////////////////////////////////////////////////////////////////////////////////////////
PixelFEDIanaCalibration::PixelFEDIanaCalibration(const PixelFEDSupervisorConfiguration & tempConfiguration, SOAPCommander* mySOAPCmdr)
  : PixelFEDCalibrationBase(tempConfiguration,*mySOAPCmdr), rootf(0)
{
  std::cout << "In PixelFEDIanaCalibration copy ctor()" << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////////
void PixelFEDIanaCalibration::initializeFED() {
  setFEDModeAndControlRegister(0x8, 0x30010);
  //setFEDModeAndControlRegister(0x8, 0x00014);
  printIfSlinkHeaderMessedup_off();
  sendResets();
  //setFEDModeAndControlRegister(0x8, 0x10015);
  
}

///////////////////////////////////////////////////////////////////////////////////////////////
xoap::MessageReference PixelFEDIanaCalibration::beginCalibration(xoap::MessageReference msg) {
  std::cout << "In PixelFEDIanaCalibration::beginCalibration()" << std::endl;

  PixelCalibConfiguration* tempCalibObject = dynamic_cast<PixelCalibConfiguration*>(theCalibObject_);
  assert(tempCalibObject != 0);

  tempCalibObject->writeASCII(outputDir());

  setFIFO1Mode();//jen

  vector<PixelModuleName>::const_iterator module_name = theDetectorConfiguration_->getModuleList().begin();
  for (;module_name!=theDetectorConfiguration_->getModuleList().end();++module_name){
    // First we need to get the DAC settings for the ROCs on this module.
    PixelDACSettings *dacs=0; 
    string modulePath=module_name->modulename();
    PixelConfigInterface::get(dacs, "pixel/dac/"+modulePath, *theGlobalKey_);
    assert(dacs!=0);
    dacsettings_[*module_name]=dacs;
  }
  
  inject_ = false;
  const std::vector<std::vector<unsigned int> > cols = tempCalibObject->columnList();
  const std::vector<std::vector<unsigned int> > rows = tempCalibObject->rowList();
  if( cols[0].size() != 0 && rows[0].size() != 0 ) inject_ = true;

  for (unsigned dacnum = 0; dacnum < tempCalibObject->numberOfScanVariables(); ++dacnum) {
    const std::string& dacname = tempCalibObject->scanName(dacnum);
    std::vector<unsigned int> dacvals = tempCalibObject->scanValues(dacname);
    if (dacvals.size() > 1)
      dacsToScan.push_back(dacname);

    for( unsigned int i = 0; i < dacvals.size(); ++i ) std::cout << " dac value " << i << " is " << dacvals[i] << std::endl;
  }

  BookEm("");

  xoap::MessageReference reply = MakeSOAPMessageReference("BeginCalibrationDone");
  return reply;
}

///////////////////////////////////////////////////////////////////////////////////////////////
xoap::MessageReference PixelFEDIanaCalibration::execute(xoap::MessageReference msg) {
  Attribute_Vector parameters(2);
  parameters[0].name_ = "WhatToDo";
  parameters[1].name_ = "StateNum";
  Receive(msg, parameters);

  const unsigned state = atoi(parameters[1].value_.c_str());
  
  if (parameters[0].value_ == "RetrieveData")
    RetrieveData(state);
  else if (parameters[0].value_ == "Analyze")
    Analyze();
  else {
    cout << "ERROR: PixelFEDIanaCalibration::execute() does not understand the WhatToDo command, "<< parameters[0].value_ <<", sent to it.\n";
    assert(0);
  }

  xoap::MessageReference reply = MakeSOAPMessageReference("FEDCalibrationsDone");
  return reply;
}

///////////////////////////////////////////////////////////////////////////////////////////////
xoap::MessageReference PixelFEDIanaCalibration::endCalibration(xoap::MessageReference msg) {

  std::cout << "In PixelFEDIanaCalibration::endCalibration()" << std::endl;
  xoap::MessageReference reply = MakeSOAPMessageReference("EndCalibrationDone");
  return reply;
}

///////////////////////////////////////////////////////////////////////////////////////////////
void PixelFEDIanaCalibration::RetrieveData(unsigned state) {
  PixelCalibConfiguration* tempCalibObject = dynamic_cast<PixelCalibConfiguration*>(theCalibObject_);
  assert(tempCalibObject != 0);

  std::map<std::string, unsigned int> currentDACValues;
  for (unsigned dacnum = 0; dacnum < tempCalibObject->numberOfScanVariables(); ++dacnum) {
    const std::string& dacname = tempCalibObject->scanName(dacnum);
    const unsigned dacvalue = tempCalibObject->scanValue(tempCalibObject->scanName(dacnum), state);
    currentDACValues[dacname] = dacvalue;
  }
  if(currentDACValues["Vana"] != lastVana){
   event_ = 0;
   lastVana = currentDACValues["Vana"];
   last_dacs.clear();
  }
  
  const std::vector<std::pair<unsigned, std::vector<unsigned> > >& fedsAndChannels = tempCalibObject->fedCardsAndChannels(crate_, theNameTranslation_, theFEDConfiguration_, theDetectorConfiguration_);
  
  for (unsigned ifed = 0; ifed < fedsAndChannels.size(); ++ifed) {

    const unsigned fednumber = fedsAndChannels[ifed].first;
    const unsigned long vmeBaseAddress = theFEDConfiguration_->VMEBaseAddressFromFEDNumber(fednumber);
    PixelFEDInterface* iFED = FEDInterface_[vmeBaseAddress];

    const int MaxChans = 37;
    uint32_t bufferFifo1[MaxChans][1024];
    int statusFifo1[MaxChans] = {0};
    
    for( unsigned int ch = 0; ch < fedsAndChannels[ifed].second.size(); ch++ ){
     statusFifo1[ch] = iFED->drainFifo1(fedsAndChannels[ifed].second[ch], bufferFifo1[ch], 1024);
    }

    for( unsigned int ch = 0; ch < fedsAndChannels[ifed].second.size(); ch++ ){
      
     int channel = (fedsAndChannels[ifed].second)[ch];
           
     if (statusFifo1[ch] > 0){ 

      for( int i = 0; i < statusFifo1[ch]; ++i ){
       uint32_t w = bufferFifo1[ch][i];
       int decCh = (w >> 26) & 0x3f;       
       if( decCh != channel ) continue;       
       uint32_t mk = (w >> 21) & 0x1f;
       uint32_t az = (w >> 8) & 0x1fff;
       uint32_t dc = (w >> 16) & 0x1f;
       uint32_t px = (w >> 8) & 0xff;
       uint32_t f8 = w & 0xff;
       if (az == 0 && mk != 0x1e ){
        if( mk <= 8 ){
         last_dacs[fedsAndChannels[ifed].first][channel][mk].push_back(f8);
        }
       }// close reading dac from roc header
      }// close fifo1 reading  
     }// close if status > 0
     
    }// close loop on channels

  }// close loop on feds

  if( event_ == 31 ){

     ReadIana();
  }

  event_++;
  sendResets();
    
}

///////////////////////////////////////////////////////////////////////////////////////////////
void PixelFEDIanaCalibration::Analyze() {

  ofstream out((outputDir()+"/iana.dat").c_str()); //leave the file method intact for now
  assert(out.good()); //file method
  
  branch theBranch;
  branch_sum theBranch_sum;
  TDirectory* dirSummaries = gDirectory->mkdir("SummaryTrees","SummaryTrees");
  dirSummaries->cd();

  TTree* tree = new TTree("PassState","PassState");
  TTree* tree_sum =new TTree("SummaryInfo","SummaryInfo");
  
  tree->Branch("PassState",&theBranch,"pass/F:rocName/C",4096000);
  tree_sum->Branch("SummaryInfo",&theBranch_sum,"deltaVana/F:newVana/F:newIana/F:maxIana/F:fitChisquare/F:rocName/C",4096000);
  rootf->cd();

  PixelCalibConfiguration* tempCalibObject = dynamic_cast<PixelCalibConfiguration*>(theCalibObject_);
  assert(tempCalibObject != 0);

  PixelRootDirectoryMaker rootDirs(tempCalibObject->rocList(),gDirectory);


  for( std::map<int,std::map<int,std::map<int,std::map<int,int> > > >::iterator it = values.begin(); it!=values.end(); ++it ){
   
   for( std::map<int,std::map<int,std::map<int,int> > >::iterator it2 = (it->second).begin(); it2!=(it->second).end(); ++it2 ){
   
    for( std::map<int,std::map<int,int> >::iterator it3 = (it2->second).begin(); it3!=(it2->second).end(); ++it3 ){
    
       if( theNameTranslation_->ROCNameFromFEDChannelROCExists(it->first,it2->first,it3->first-1) ){
       
        PixelROCName theROC = theNameTranslation_->ROCNameFromFEDChannelROC(it->first,it2->first,it3->first-1);
        PixelModuleName theModule(theROC.rocname());
        
        theBranch.pass = 0;
        strcpy(theBranch.rocName, theROC.rocname().c_str());
        strcpy(theBranch_sum.rocName, theROC.rocname().c_str());
        
        int npoints = it3->second.size();
        std::vector<double> x(npoints+1), y(npoints+1), ey(npoints+1);

        int j = 0;
        for( std::map<int,int>::iterator it4 = (it3->second).begin(); it4!=(it3->second).end(); ++it4 ){
	     y[j] = 0.25*it4->second;
	     x[j] = it4->first;
	     ey[j] = 1.;
	     j++;
        }
        
        const int oldVana = dacsettings_[theModule]->getDACSettings(theROC)->getVana();

        rootDirs.cdDirectory(theROC);
        PixelIanaAnalysis analysis(true);
        analysis.go(theROC.rocname(),oldVana,npoints,x, y, ey,out);

        theBranch_sum.maxIana = analysis.maxIana;
        theBranch_sum.fitChisquare = analysis.fitChisquare;
        theBranch_sum.newVana = analysis.newVana;
        theBranch_sum.newIana = analysis.newIana;
        theBranch_sum.deltaVana = theBranch_sum.newVana - oldVana;
        theBranch.pass = analysis.pass;

        tree->Fill();
        tree_sum->Fill();
                      
       }//close loop 4					      
    }//close for3    
   }//close for2   
  }//close for 1    
  
  CloseRootf();

}

///////////////////////////////////////////////////////////////////////////////////////////////
void PixelFEDIanaCalibration::CloseRootf() {
  if (rootf) {
    rootf->Write();
    rootf->Close();
    delete rootf;
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////
void PixelFEDIanaCalibration::BookEm(const TString& path) {

  values.clear();
  
  TString root_fn;
  root_fn.Form("%s/Iana.root", outputDir().c_str());
  cout << "writing histograms to file " << root_fn << endl;
  CloseRootf();
  rootf = new TFile(root_fn, "create");
  assert(rootf->IsOpen());
  rootf->cd(0);

}
  
///////////////////////////////////////////////////////////////////////////////////////////////
void PixelFEDIanaCalibration::ReadIana( void ) {
  
  std::cout << " **************** VANA = " << lastVana << std::endl;
  for( std::map<int,std::map<int,std::map<int,std::vector<int> > > >::iterator it = last_dacs.begin(); it!=last_dacs.end(); ++it ){
   
   for( std::map<int,std::map<int,std::vector<int> > >::iterator it2 = (it->second).begin(); it2!=(it->second).end(); ++it2 ){
   
    for( std::map<int,std::vector<int> >::iterator it3 = (it2->second).begin(); it3!=(it2->second).end(); ++it3 ){

      int start = -1;
      int stop = -1;
      int iana = 0;
      std::vector<int> bits;
       
      for( unsigned int i=0; i<(it3->second).size(); ++i ){

       //std::cout << " - channel#" << it2->first << " ROC#" << it3->first << " - dac " << i << ":" << (it3->second)[i] << std::endl;
       if( (((it3->second)[i])&(0x2)) == 2 && start == -1 ){start = i;}    
       else if( (((it3->second)[i])&(0x2)) == 2 && start != -1 ){stop = i;}         
      }
      
       std::cout << "    (found start =  " << start << " and stop = " << stop << " ) " << std::endl;     
       if( start != stop && (stop-start)==16){
   
        //now put together the 16 bits
        for( int i = start+1; i <= stop; ++i ){ bits.push_back((it3->second)[i]&0x1);}
    
        //now get Iana
        for( unsigned int b=0; b<bits.size(); ++b ){ std::cout << bits[b];}
        std::cout << std::endl;
        for( unsigned int b=8; b<bits.size(); ++b ){
         //std::cout << bits[b];
         iana+=(bits[b]<<(16-b));
        }
    
       }
   
       //std::cout << std::endl;
       //std::cout << " **************** VANA = " << lastVana << " FOUND IANA: " << iana << std::endl;
       
       values[it->first][it2->first][it3->first][lastVana] = iana;            
      
    }      
   }      
  }//close loop

}
