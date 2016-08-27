#include "CalibFormats/SiPixelObjects/interface/PixelCalibConfiguration.h"
#include "CalibFormats/SiPixelObjects/interface/PixelDACNames.h"
#include "PixelCalibrations/include/PixelFEDTBMDelayCalibration.h"
#include "PixelConfigDBInterface/include/PixelConfigInterface.h"
#include "PixelUtilities/PixelFEDDataTools/include/PixelFEDDataTypes.h"
#include "PixelUtilities/PixelFEDDataTools/include/ErrorFIFODecoder.h"
#include "PixelUtilities/PixelFEDDataTools/include/ColRowAddrDecoder.h"
#include "PixelUtilities/PixelFEDDataTools/include/DigScopeDecoder.h"
#include "PixelUtilities/PixelFEDDataTools/include/DigTransDecoder.h"
#include "PixelUtilities/PixelFEDDataTools/include/FIFO3Decoder.h"
#include "PixelUtilities/PixelRootUtilities/include/PixelRootDirectoryMaker.h"
#include "PixelUtilities/PixelFEDDataTools/include/DigFIFO1Decoder.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include <iomanip>
#include <algorithm>

using namespace pos;

///////////////////////////////////////////////////////////////////////////////////////////////
PixelFEDTBMDelayCalibration::PixelFEDTBMDelayCalibration(const PixelFEDSupervisorConfiguration & tempConfiguration, SOAPCommander* mySOAPCmdr)
  : PixelFEDCalibrationBase(tempConfiguration,*mySOAPCmdr), rootf(0)
{
  std::cout << "In PixelFEDTBMDelayCalibration copy ctor()" << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////////
void PixelFEDTBMDelayCalibration::initializeFED() {
  setFEDModeAndControlRegister(0x8, 0x30010);
  //setFEDModeAndControlRegister(0x8, 0x00014);
  printIfSlinkHeaderMessedup_off();
  sendResets();
  //setFEDModeAndControlRegister(0x8, 0x10015);
  
}

///////////////////////////////////////////////////////////////////////////////////////////////
xoap::MessageReference PixelFEDTBMDelayCalibration::beginCalibration(xoap::MessageReference msg) {
  std::cout << "In PixelFEDTBMDelayCalibration::beginCalibration()" << std::endl;

  PixelCalibConfiguration* tempCalibObject = dynamic_cast<PixelCalibConfiguration*>(theCalibObject_);
  assert(tempCalibObject != 0);

  tempCalibObject->writeASCII(outputDir());

  DumpFIFOs = tempCalibObject->parameterValue("DumpFIFOs") == "yes";
  ReadFifo1 = tempCalibObject->parameterValue("ReadFifo1") == "yes";
  ReadFifo3 = tempCalibObject->parameterValue("ReadFifo3") == "yes";
  //  const std::vector<PixelROCName>& rocs = tempCalibObject->rocList();
  //PixelRootDirectoryMaker rootDirs(rocs, rootf);

  if( ReadFifo1 ) setFIFO1Mode();//jen

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
xoap::MessageReference PixelFEDTBMDelayCalibration::execute(xoap::MessageReference msg) {
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
    cout << "ERROR: PixelFEDTBMDelayCalibration::execute() does not understand the WhatToDo command, "<< parameters[0].value_ <<", sent to it.\n";
    assert(0);
  }

  xoap::MessageReference reply = MakeSOAPMessageReference("FEDCalibrationsDone");
  return reply;
}

///////////////////////////////////////////////////////////////////////////////////////////////
xoap::MessageReference PixelFEDTBMDelayCalibration::endCalibration(xoap::MessageReference msg) {

  std::cout << "In PixelFEDTBMDelayCalibration::endCalibration()" << std::endl;
  xoap::MessageReference reply = MakeSOAPMessageReference("EndCalibrationDone");
  return reply;
}

///////////////////////////////////////////////////////////////////////////////////////////////
void PixelFEDTBMDelayCalibration::RetrieveData(unsigned state) {
  PixelCalibConfiguration* tempCalibObject = dynamic_cast<PixelCalibConfiguration*>(theCalibObject_);
  assert(tempCalibObject != 0);

  /*const std::vector<PixelROCName>& rocs = tempCalibObject->rocList();
  typedef std::set< std::pair<unsigned int, unsigned int> > colrow_t;
  const colrow_t colrows = tempCalibObject->pixelsWithHits(state);
  if (PrintHits) {
    std::cout << "ZZ ";
    for (colrow_t::const_iterator cr = colrows.begin(); cr != colrows.end(); ++cr)
      std::cout << "c " << cr->first << " r " << cr->second << " ";
    std::cout << std::endl;
  }*/

  const std::vector<std::pair<unsigned, std::vector<unsigned> > >& fedsAndChannels = tempCalibObject->fedCardsAndChannels(crate_, theNameTranslation_, theFEDConfiguration_, theDetectorConfiguration_);

  if (DumpFIFOs) std::cout << "NEW FEDTBMDelay TRIGGER " << event_ << " state " << state << " ";
  std::map<std::string, unsigned int> currentDACValues;
  for (unsigned dacnum = 0; dacnum < tempCalibObject->numberOfScanVariables(); ++dacnum) {
    const std::string& dacname = tempCalibObject->scanName(dacnum);
    const unsigned dacvalue = tempCalibObject->scanValue(tempCalibObject->scanName(dacnum), state);
    currentDACValues[dacname] = dacvalue;
    if (DumpFIFOs) std::cout << dacname << " " << dacvalue << " ";
  }
  if (DumpFIFOs) std::cout << std::endl;
  if(dacsToScan.size() < 2 && currentDACValues["TBMPLL"] != lastTBMPLL){
   event_ = 0;
   lastTBMPLL = currentDACValues["TBMPLL"];
  }

  //////

  for (unsigned ifed = 0; ifed < fedsAndChannels.size(); ++ifed) {
    const unsigned fednumber = fedsAndChannels[ifed].first;
    const unsigned long vmeBaseAddress = theFEDConfiguration_->VMEBaseAddressFromFEDNumber(fednumber);
    PixelFEDInterface* iFED = FEDInterface_[vmeBaseAddress];
    iFED->readDigFEDStatus(false, false);

    //const uint32_t fifoStatus = iFED->getFifoStatus();

    const int MaxChans = 37;    
    uint64_t buffer3[2048];
    uint32_t bufferErr[36*1024];
    DigTransDecoder* decodeT[MaxChans] = {0};
    DigScopeDecoder* decodeS[MaxChans] = {0};
    FIFO3Decoder* decode3 = 0;
    //const int status3;
    //const int statusErr;
    uint32_t bufferFifo1[MaxChans][1024];
    int statusFifo1[MaxChans] = {0};

    //iFED->SetFitelFiberSwitchTopDauCard(0); // this should be configurable from outside
    //iFED->SetFitelFiberSwitchBottomDauCard(0);

    /* read fifo1 */
    if( ReadFifo1 ){

      for( unsigned int ch = 0; ch < fedsAndChannels[ifed].second.size(); ch++ ){
       statusFifo1[ch] = iFED->drainFifo1(fedsAndChannels[ifed].second[ch], bufferFifo1[ch], 1024);
      }

      for( unsigned int ch = 0; ch < fedsAndChannels[ifed].second.size(); ch++ ){

       int channel = (fedsAndChannels[ifed].second)[ch];
       bool found_TBMA = false;
       std::vector<int> ch_decodedROCs;
       bool ch_foundHit = false;

       if (statusFifo1[ch] > 0) {

        DigFIFO1Decoder theFIFO1Decoder(bufferFifo1[ch],statusFifo1[ch]);
        if( theFIFO1Decoder.globalChannel() != channel ) continue;
        found_TBMA = theFIFO1Decoder.foundTBM();
        if( !inject_ ) ch_decodedROCs = theFIFO1Decoder.ROCHeaders();
        else{
         for( unsigned int h = 0; h < theFIFO1Decoder.nhits(); ++h ){
          if(std::find(ch_decodedROCs.begin(),ch_decodedROCs.end(),theFIFO1Decoder.rocid(h))==ch_decodedROCs.end()) ch_decodedROCs.push_back(theFIFO1Decoder.rocid(h));
         }
        } 

        if( DumpFIFOs ){
         std::cout << "-----------------------------------" << std::endl;
         std::cout << "Contents of FIFO 1 for channel " << channel << " (status = " << statusFifo1[ch] << ")" << std::endl;
         std::cout << "-----------------------------------" << std::endl;
         theFIFO1Decoder.printToStream(std::cout);
        }

       }

       ch_foundHit = (ch_decodedROCs.size() == 4);
       FillEm(state, fedsAndChannels[ifed].first, channel, 0, (!inject_ && found_TBMA) || (inject_ && found_TBMA && ch_foundHit) );

       if( dacsToScan.size() == 0 ){
        for( int r = 0; r < 8; ++r ){

         bool ch_foundROC = std::find(ch_decodedROCs.begin(),ch_decodedROCs.end(),r+1)!=ch_decodedROCs.end();
         if( ch_foundROC ) FillEm(state, fedsAndChannels[ifed].first, channel, 1, r);
        }
       }
       else if( dacsToScan.size() == 1){
        FillEm(state, fedsAndChannels[ifed].first, channel, 1, ch_decodedROCs.size());
       }

      }// end loop on channels

    }//end readFifo1

    uint32_t bufferS[MaxChans][256];
    int statusS[MaxChans] = {0};

    /* read scope fifo */
    /*if( true ){
     //read
     for( unsigned int ch = 0; ch < fedsAndChannels[ifed].second.size(); ch++ ){     
      iFED->SelectScopeChannel(fedsAndChannels[ifed].second[ch]);
      statusS[ch] = iFED->drainDataChannelFifo2(fedsAndChannels[ifed].second[ch], bufferS[ch]);
      decodeS[ch] = new DigScopeDecoder(bufferS[ch], statusS[ch]);
     }
     //fill histos
     /*for( unsigned int ch = 0; ch < fedsAndChannels[ifed].second.size(); ch++ ){

      DigScopeDecoder* d = decodeS[ch];

      FillEm(state, fedsAndChannels[ifed].first, (fedsAndChannels[ifed].second)[ch], 0, d->tbm_header_found_ && d->tbm_trailer_found_ );
      FillEm(state, fedsAndChannels[ifed].first, (fedsAndChannels[ifed].second)[ch], 1, d->roc_headers_.size());

     }
     //print content
     if (DumpFIFOs) {
      int colS=-1, rowS=-1;
      for( unsigned int ch = 0; ch < fedsAndChannels[ifed].second.size(); ch++ ){
       std::cout << "----------------------------------" << std::endl;
       if (statusS[ch] < 0) std::cout << "Scope FIFO for FED " << fednumber << " channel " << (fedsAndChannels[ifed].second)[ch] << " status = " << statusS[ch] << std::endl;
       else {
        std::cout << "Contents of Scope FIFO for FED " << fednumber << " channel " << (fedsAndChannels[ifed].second)[ch] << " (statusS = " << statusS[ch] << ")" <<std::endl;
       std::cout << "----------------------------------" << std::endl;
        for (int i = 0; i <= statusS[ch]; ++i) {
         uint32_t d = bufferS[ch][i];
         uint32_t dh = d & 0xf0;
         if (dh == 0x70 || dh == 0x10 || dh == 0xc0) std::cout << "\n";
         if (d > 0xFF) std::cout << "\nweird word: " << std::hex << d << "\n";
         else std::cout << std::setw(2) << std::hex << d << std::dec << " ";
        }// end read buffer for fifo2
        std::cout << "\n----------------------------------" << std::endl;
       }//close else
       std::cout << "DigScopeDecoder thinks:\n";
       decodeS[ch]->printToStream(std::cout);
       if (decodeS[ch]->n_hits() > 6) {
        colS = decodeS[ch]->hits()[0].col;
        rowS = decodeS[ch]->hits()[0].row;
       }
      }//end loop on channels
     }// end dumping fifo content
    }//end ReadScopeFifo case*/

    if( ReadFifo3 ){
     //read
     const int status3 = iFED->spySlink64(buffer3);
     const int statusErr = iFED->drainErrorFifo(bufferErr);
     if (status3 > 0) decode3 = new FIFO3Decoder(buffer3);
     ErrorFIFODecoder decodeErr(bufferErr, statusErr);
     //print hits and fill histos
     /*if (PrintHits) std::cout << "F3X ";
     if (status3 <= 0) FillEm(state, F3fifoErr, 1);
     else {
      for (unsigned ihit = 0; ihit < decode3->nhits(); ++ihit) {
	const unsigned channel = decode3->channel(ihit);
	const unsigned rocid = decode3->rocid(ihit);
	assert(rocid > 0);

	const PixelROCName& roc = theNameTranslation_->ROCNameFromFEDChannelROC(fednumber, channel, rocid-1);

	// Skip if this ROC is not on the list of ROCs to calibrate.
	// Also skip if we're in singleROC mode, and this ROC is not being calibrated right now.
	vector<PixelROCName>::const_iterator foundROC = find(rocs.begin(), rocs.end(), roc);
	if (foundROC == rocs.end()) // || !tempCalibObject->scanningROCForState(roc, state))
	  FillEm(state, F3wrongRoc, 1);
	else {
	  const unsigned col = decode3->column(ihit);
	  const unsigned row = decode3->row(ihit);
	  if (PrintHits) std::cout << "c " << col << " r " << row << " ";
	  if (colrows.find(std::make_pair(col, row)) == colrows.end())
	    FillEm(state, F3wrongPix, 1);
	  else
	    FillEm(state, F3rightPix, 1);
	}
      }
     }*/
     //print content
     if (DumpFIFOs) {
      for( unsigned int ch = 0; ch < fedsAndChannels[ifed].second.size(); ch++ ){
       std::cout << "----------------------" << std::endl;
       std::cout << "Contents of Spy FIFO 3" << std::endl;
       std::cout << "----------------------" << std::endl;
       for (int i = 0; i <= status3; ++i) std::cout << "Clock " << std::setw(2) << i << " = 0x " << std::hex << std::setw(8) << (buffer3[i]>>32) << " " << std::setw(8) << (buffer3[i] & 0xFFFFFFFF) << std::dec << std::endl;
       if (status3 > 0) {
        std::cout << "FIFO3Decoder thinks:\n" << "nhits: " << decode3->nhits() << std::endl;
        int hits_by_ch[37] = {0};
        int hits_by_roc[37][8] = {{0}};
        unsigned lastroc = 0;
        for (unsigned i = 0; i < decode3->nhits(); ++i) {
         if( decode3->channel(i) != (fedsAndChannels[ifed].second)[ch] ) continue;
         const PixelROCName& rocname = theNameTranslation_->ROCNameFromFEDChannelROC(fednumber, decode3->channel(i), decode3->rocid(i)-1);
         ++hits_by_ch[decode3->channel(i)];
         ++hits_by_roc[decode3->channel(i)][decode3->rocid(i)-1];
         if (lastroc != 0 && decode3->rocid(i) != lastroc) {
	  std::cout << "\n";
	  lastroc = decode3->rocid(i);
         }
         std::cout << "#" << i << ": ch: " << decode3->channel(i)
	  	 << " rocid: " << decode3->rocid(i)
		 << " (" << rocname << ")"
		 << " dcol: " << decode3->dcol(i)
		 << " pxl: " << decode3->pxl(i) << " pulseheight: " << decode3->pulseheight(i)
		 << " col: " << decode3->column(i) << " row: " << decode3->row(i) << std::endl;
        }//end loop on nhits
        std::cout << "Nhits by channel:\n";
        for( unsigned int i = 0; i < fedsAndChannels[ifed].second.size(); i++ ){
         if (hits_by_ch[i]) std::cout << "Ch " << std::setw(2) << i << ": " << std::setw(3) << hits_by_ch[i] << "\n";
        }
        std::cout << "Nhits by roc:\n";
        for( unsigned int i = 0; i < fedsAndChannels[ifed].second.size(); i++ ){
         for (int j = 0; j < 8; ++j){
	  if (hits_by_roc[i][j]) std::cout << "Ch " << std::setw(2) << i << " roc " << j << ": " << std::setw(3) << hits_by_roc[i][j] << "\n";
         }
        }
        if (decode3->nhits() > 0) std::cout /*<< "(fifo2 col: " << colS << " row: " << rowS*/ << "   fifo3 dcol: " << decode3->dcol(0) << " pxl: " << decode3->pxl(0) << " col: " << decode3->column(0) << " row: " << decode3->row(0) << ")\n";
       }//close if status3 > 0*/       
       std::cout << "Contents of Error FIFO" << std::endl;
       std::cout << "----------------------" << std::endl;
       for (int i = 0; i <= statusErr; ++i) std::cout << "Clock " << i << " = 0x" << std::hex << bufferErr[i] << std::dec << std::endl;
       std::cout << "ErrorFIFODecoder thinks:\n";
       decodeErr.printToStream(std::cout);
      }//close loop on channels    
     }//end dumping fifo content 
    }//end read fifo 3 case
    
    for (int i = 0; i < 36; ++i) {
     delete decodeT[i];
     delete decodeS[i];
    }
    delete decode3;
  }

  event_++;
  sendResets();

}

///////////////////////////////////////////////////////////////////////////////////////////////
void PixelFEDTBMDelayCalibration::Analyze() {

  PixelCalibConfiguration* tempCalibObject = dynamic_cast<PixelCalibConfiguration*>(theCalibObject_);
  assert(tempCalibObject != 0);
  int ntriggers = event_-1;
  std::map<std::string,int> bestTBMPLLSettings;

  if (dacsToScan.size() == 0){

    for( std::map<int,std::map<int,std::vector<TH1F*> > >::iterator it1 = ntrigsTBM.begin(); it1 != ntrigsTBM.end(); ++it1 ){
     for( std::map<int,std::vector<TH1F*> >::iterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2 ){
      for(unsigned int i = 0; i < it2->second.size(); ++i ) it2->second[i]->Scale(1./ntriggers);
     }
    }

  }

  if (dacsToScan.size() == 1){

   //fill histo with sum of channels
   for( std::map<int,std::map<int,std::vector<TH2F*> > >::iterator it1 = scansTBM.begin(); it1 != scansTBM.end(); ++it1 ){
    std::string moduleName = "";
    PixelChannel theChannel;
    for( std::map<int,std::vector<TH2F*> >::iterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2 ){  
     
     if( theNameTranslation_->FEDChannelExist(it1->first, it2->first) ){
      theChannel = theNameTranslation_->ChannelFromFEDChannel(it1->first, it2->first);
      moduleName = theChannel.modulename();
     }

     //std::cout << "*********************** ::analyze() " << it1->first << " " << it2->first << " " << moduleName << std::endl;
     for(unsigned int i = 0; i < it2->second.size(); ++i ) TBMsHistoSum[it1->first][moduleName][i]->Add(it2->second[i]);

    }//close loop on channels
   }//close loop on fed

   //normalize by number of triggers
   for( std::map<int,std::map<int,std::vector<TH2F*> > >::iterator it1 = scansTBM.begin(); it1 != scansTBM.end(); ++it1 ){
    for( std::map<int,std::vector<TH2F*> >::iterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2 ){
     for(unsigned int i = 0; i < it2->second.size(); ++i ) it2->second[i]->Scale(1./ntriggers);
    }
   }

   for( std::map<int,std::map<std::string,std::vector<TH2F*> > >::iterator it1 = TBMsHistoSum.begin(); it1 != TBMsHistoSum.end(); ++it1 ){
    for( std::map<std::string,std::vector<TH2F*> >::iterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2 ){
     for(unsigned int i = 0; i < it2->second.size(); ++i ) it2->second[i]->Scale(1./(4*ntriggers));
    }
   }

   //find best settings for each module
   std::cout << "******************************************************" << std::endl;
   for( std::map<int,std::map<std::string,std::vector<TH2F*> > >::iterator it1 = TBMsHistoSum.begin(); it1 != TBMsHistoSum.end(); ++it1 ){
    for( std::map<std::string,std::vector<TH2F*> >::iterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2 ){
     
     float eff = 0;
     std::map<int,int> bestBins;
     for( int bx = 1; bx < it2->second[0]->GetNbinsX()+1; ++bx ){
      for( int by = 1; by < it2->second[0]->GetNbinsY()+1; ++by ){
       
       //std::cout << bx << " " << by << " " << it2->second[0]->GetBinContent(bx,by) << " " << it2->second[1]->GetBinContent(bx,by) << std::endl;
       if( it2->second[0]->GetBinContent(bx,by) >= eff ){
        eff = it2->second[0]->GetBinContent(bx,by);
        bestBins[bx-1] = by-1;
        //std::cout << "******************* FOUND BEST BIN " << bx-1 << " " << by-1 << " " << eff << " " << it2->second[1]->GetBinContent(bx,by) << std::endl;
       }

      }
     }

     std::cout << "RESULTS OF TBMPLLDELAY SCAN FOR MODULE " << it2->first << std::endl;
     if( bestBins.size() == 0 ) std::cout << " --- NO OPTIMAL TBMPLLDELAY VALUE FOUND!" << std::endl;

     int bestX = 0;
     int bestY = 0;
     for( std::map<int,int>::iterator binsIt = bestBins.begin(); binsIt != bestBins.end(); ++binsIt ){

      //std::cout << "******************* FOUND BEST BIN first " << binsIt->first << " second " << binsIt->second << std::endl;

      if( binsIt->first >= 1 && binsIt->first <= 6 && ( (binsIt->second >= 5 && binsIt->second <= 6) || binsIt->second == 1 ) ){
       bestX = binsIt->first;
       bestY = binsIt->second;
       break;
      }

     }


     if( bestBins.size() != 0 && bestX == 0 && bestY == 0 ){
      bestX = (bestBins.begin())->first;
      bestY = (bestBins.begin())->second;
      std::cout << " --- WARNING: BEST TBMPLLDELAY VALUE IS AT THE EDGES! 160 MHZ PLL DELAY = " << bestY;
      std::cout << " ; 400 MHZ PLL DELAY = " << bestX;
      
      bestX = (bestX<<2);
      bestY = (bestY<<5);
      std::cout << " --> NEW SETTINGS TBMPLLDelay: " << bestX+bestY << std::endl;
     }
     else if( bestX != 0 && bestY != 0 ){
      std::cout << " --- FOUND BEST TBMPLLDELAY VALUE! 160 MHZ PLL DELAY = " << bestY;
      std::cout << " 400 MHZ PLL DELAY = " << bestX;
      
      bestX = (bestX<<2);
      bestY = (bestY<<5);
      std::cout << " --> NEW SETTINGS TBMPLLDelay: " << bestX+bestY << std::endl;
     }

     bestTBMPLLSettings[it2->first] = bestX+bestY;

    }
   }      

  }

  CloseRootf();

  if( dacsToScan.size() == 0 ) return;

  std::vector<PixelModuleName>::const_iterator module_name = theDetectorConfiguration_->getModuleList().begin();
  for (;module_name!=theDetectorConfiguration_->getModuleList().end();++module_name)
  {
    PixelTBMSettings *TBMSettingsForThisModule=0;
    std::string moduleNameString=(module_name->modulename());
    PixelConfigInterface::get(TBMSettingsForThisModule, "pixel/tbm/"+moduleNameString, *theGlobalKey_);
    assert(TBMSettingsForThisModule!=0);

    std::string moduleName = module_name->modulename();
    if( bestTBMPLLSettings[moduleName] != 0 ) TBMSettingsForThisModule->setTBMPLLDelay(bestTBMPLLSettings[moduleName]);

    TBMSettingsForThisModule->writeASCII(outputDir());
    std::cout << "Wrote TBM settings for module:" << moduleName << endl;
			
    delete TBMSettingsForThisModule;
  }

}

///////////////////////////////////////////////////////////////////////////////////////////////
void PixelFEDTBMDelayCalibration::CloseRootf() {
  if (rootf) {
    rootf->Write();
    rootf->Close();
    delete rootf;
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////
void PixelFEDTBMDelayCalibration::BookEm(const TString& path) {

  ntrigsTBM.clear();
  scansTBM.clear();
  TBMsHistoSum.clear();

  TString root_fn;
  if (path == "")
    root_fn.Form("%s/TBMDelay.root", outputDir().c_str());
  else
    root_fn.Form("%s/TBMDelay_%s.root", outputDir().c_str(), path.Data());
  cout << "writing histograms to file " << root_fn << endl;
  CloseRootf();
  rootf = new TFile(root_fn, "create");
  assert(rootf->IsOpen());

  PixelCalibConfiguration* tempCalibObject = dynamic_cast<PixelCalibConfiguration*>(theCalibObject_);
  assert(tempCalibObject != 0);
  const std::vector<std::pair<unsigned, std::vector<unsigned> > >& fedsAndChannels = tempCalibObject->fedCardsAndChannels(crate_, theNameTranslation_, theFEDConfiguration_, theDetectorConfiguration_);

  for (unsigned ifed = 0; ifed < fedsAndChannels.size(); ++ifed) {
   TString FEDdir; FEDdir.Form("FED%i",fedsAndChannels[ifed].first);
   rootf->mkdir(FEDdir);
   rootf->cd(FEDdir);

   TDirectory* dir = rootf->GetDirectory(FEDdir);
   std::map<int,std::vector<TH1F*> > chTBMmap;
   std::map<int,std::vector<TH2F*> > chTBMmap2D;
   for( unsigned int ch = 0; ch < fedsAndChannels[ifed].second.size(); ch++ ){
    TString chDir; chDir.Form("Channel%i",(fedsAndChannels[ifed].second)[ch]);
    dir->mkdir(chDir);
    dir->cd(chDir);

    if (dacsToScan.size() == 0){

     TString hname; hname.Form("Ch%i",(fedsAndChannels[ifed].second)[ch]);
     std::vector<TH1F*> histosTBM;
     TH1F* h_TBM_nDecodes = new TH1F(hname+"_TBM_nDecodes", hname+"_TBM_nDecodes", 2, 0, 2 );
     histosTBM.push_back(h_TBM_nDecodes);
     TH1F* h_nROCHeaders = new TH1F(hname+"_nROCHeaders", hname+"_nROCHeaders", 9, 0, 9 );
     histosTBM.push_back(h_nROCHeaders);  
     chTBMmap[(fedsAndChannels[ifed].second)[ch]] = histosTBM;
  
    }// end book histos for zero dacsToScan case

    if(dacsToScan.size() == 1){

      TString hname; hname.Form("Ch%i",(fedsAndChannels[ifed].second)[ch]);
      std::vector<TH2F*> histosTBM;
      TH2F* h_TBM_nDecodes = new TH2F(hname+"_TBM_nDecodes", hname+"_TBM_nDecodes", 8, 0, 8, 8, 0, 8 );
      histosTBM.push_back(h_TBM_nDecodes);
      TH2F* h_nROCHeaders = new TH2F(hname+"_nROCHeaders", hname+"_nROCHeaders", 8, 0, 8, 8, 0, 8 );
      histosTBM.push_back(h_nROCHeaders);     
      chTBMmap2D[(fedsAndChannels[ifed].second)[ch]] = histosTBM;
  
    }// end book histos for 1 dacsToScan case (TBMPLL scan)

   }//close loop on channels

   if (dacsToScan.size() == 0) ntrigsTBM[fedsAndChannels[ifed].first] = chTBMmap;
   if (dacsToScan.size() == 1) scansTBM[fedsAndChannels[ifed].first] = chTBMmap2D;

  }//close loop on feds

  //book histos with sum of channels
  if (dacsToScan.size() == 1){

   for (unsigned ifed = 0; ifed < fedsAndChannels.size(); ++ifed) {

    TString FEDdir; FEDdir.Form("FED%i",fedsAndChannels[ifed].first);
    rootf->cd(FEDdir);

    std::string moduleName = "";
    for( unsigned int ch = 0; ch < fedsAndChannels[ifed].second.size(); ch++ ){

     PixelChannel theChannel = theNameTranslation_->ChannelFromFEDChannel(fedsAndChannels[ifed].first, (fedsAndChannels[ifed].second)[ch]);
     if( moduleName != theChannel.modulename() ){
      moduleName = theChannel.modulename();
      std::vector<TH2F*> histosTBM;
      TString hname(moduleName);
      TH2F* h_TBM_nDecodes = new TH2F(hname+"_nTBMDecodes", hname+"_nTBMDecodes", 8, 0, 8, 8, 0, 8 );
      histosTBM.push_back(h_TBM_nDecodes);
      TH2F* h_nROCHeaders = new TH2F(hname+"_nROCHeaders", hname+"_nROCHeaders", 8, 0, 8, 8, 0, 8 );
      histosTBM.push_back(h_nROCHeaders);   
      TBMsHistoSum[fedsAndChannels[ifed].first][moduleName] = histosTBM;
     }

    }//end loop on channels
   }//end loop on feds
  }//end booking sum histo

  rootf->cd(0);

}

///////////////////////////////////////////////////////////////////////////////////////////////
void PixelFEDTBMDelayCalibration::FillEm(unsigned state, int fedid, int ch, int which, float c) {
  PixelCalibConfiguration* tempCalibObject = dynamic_cast<PixelCalibConfiguration*>(theCalibObject_);
  assert(tempCalibObject != 0);

  if (event_==0) return;

  if (dacsToScan.size() == 0 ) ntrigsTBM[fedid][ch][which]->Fill(c,1);
  if (dacsToScan.size() == 1){

   const std::string& iname = dacsToScan[0];
   const double ival(tempCalibObject->scanValue(iname, state)); 
   uint32_t tmp = ival; 
   int delay1 = (tmp>>2)&0x7;
   int delay2 = ((tmp>>2)&0x38)>>3;
   scansTBM[fedid][ch][which]->Fill(delay1,delay2,c);

  }

}
