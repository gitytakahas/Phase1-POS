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
  PrintHits = tempCalibObject->parameterValue("PrintHits") == "yes";
  ReadTransFifo = tempCalibObject->parameterValue("ReadTransFifo") == "yes";
  ReadScopeFifo = tempCalibObject->parameterValue("ReadScopeFifo") == "yes";
  ReadTmpFifo = tempCalibObject->parameterValue("ReadTmpFifo") == "yes";
  ReadFifo1 = tempCalibObject->parameterValue("ReadFifo1") == "yes";
  ReadFifo3 = tempCalibObject->parameterValue("ReadFifo3") == "yes";
  //  const std::vector<PixelROCName>& rocs = tempCalibObject->rocList();
  //PixelRootDirectoryMaker rootDirs(rocs, rootf);

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

  if (dacsToScan.empty() && tempCalibObject->parameterValue("NoScanOK") != "yes") {
    cout << "no dacs in scan?" << endl;
    assert(0);
  }

  if (dacsToScan.size() < 3)
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

xoap::MessageReference PixelFEDTBMDelayCalibration::endCalibration(xoap::MessageReference msg) {

  std::cout << "In PixelFEDTBMDelayCalibration::endCalibration()" << std::endl;
  xoap::MessageReference reply = MakeSOAPMessageReference("EndCalibrationDone");
  return reply;
}

///////////////////////////////////////////////////////////////////////////////////////////////
void PixelFEDTBMDelayCalibration::BinDisplay(unsigned long int In_Word, int Bits,int Header,int del)  

{
  int i1;
  unsigned long mask;
  
  mask=1<<(Bits-1);  
  
  if (Header) {
    //printf("\n\233\67\155\n");
    if(del) printf("|");
    for (i1=Bits;i1>0;i1--) {
      printf("%1d",(i1-1)%10);
    if(del) {if ( !((i1-1)%4) ) printf("|");}
    }
    printf("\n");
    //printf("\233\60\155 \n");
  } 
  if (del) printf("|");
  for (i1=0;i1<Bits;i1++) {
    if ((In_Word<<i1) & mask) printf("1");
    else printf("0");
    if(del) {if ( !((i1+1)%4) ) printf("|");  }
  }
  
  return;
  
}

///////////////////////////////////////////////////////////////////////////////////////////////
void PixelFEDTBMDelayCalibration::RetrieveData(unsigned state) {
  PixelCalibConfiguration* tempCalibObject = dynamic_cast<PixelCalibConfiguration*>(theCalibObject_);
  assert(tempCalibObject != 0);

  const std::vector<PixelROCName>& rocs = tempCalibObject->rocList();
  typedef std::set< std::pair<unsigned int, unsigned int> > colrow_t;
  const colrow_t colrows = tempCalibObject->pixelsWithHits(state);
  if (PrintHits) {
    std::cout << "ZZ ";
    for (colrow_t::const_iterator cr = colrows.begin(); cr != colrows.end(); ++cr)
      std::cout << "c " << cr->first << " r " << cr->second << " ";
    std::cout << std::endl;
  }

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
  /*if(dacsToScan.size() == 2 && currentDACValues["TBMADelay"] != lastTBMADelay && currentDACValues["TBMBDelay"] != lastTBMBDelay){
   event_ = 0;
   lastTBMADelay = currentDACValues["TBMADelay"];
   lastTBMBDelay = currentDACValues["TBMBDelay"];
   lastTBMPLL = currentDACValues["TBMPLL"];
  }*/
  /*if (dacsToScan.size() >= 3 && currentDACValues["TBMPLL"] != lastTBMPLL) {
    lastTBMPLL = currentDACValues["TBMPLL"];
    BookEm(TString::Format("TBMPLL%03i", lastTBMPLL));
  }*/

  //////

  for (unsigned ifed = 0; ifed < fedsAndChannels.size(); ++ifed) {
    const unsigned fednumber = fedsAndChannels[ifed].first;
    const unsigned long vmeBaseAddress = theFEDConfiguration_->VMEBaseAddressFromFEDNumber(fednumber);
    PixelFEDInterface* iFED = FEDInterface_[vmeBaseAddress];
    iFED->readDigFEDStatus(false, false);

    //const uint32_t fifoStatus = iFED->getFifoStatus();

    const int MaxChans = 37;    
    uint32_t bufferT[MaxChans][256];
    uint32_t bufferS[MaxChans][256];
    uint64_t buffer3[2048];
    uint32_t bufferErr[36*1024];
    int statusS[MaxChans] = {0};
    DigTransDecoder* decodeT[MaxChans] = {0};
    DigScopeDecoder* decodeS[MaxChans] = {0};
    FIFO3Decoder* decode3 = 0;
    //const int status3;
    //const int statusErr;
    uint32_t bufferFifo1[MaxChans][1024];
    int statusFifo1[MaxChans] = {0};
    uint32_t bufferFifo1Add[MaxChans][1024];
    int statusFifo1Add[MaxChans] = {0};

    iFED->SetFitelFiberSwitchTopDauCard(0); // this should be configurable from outside
    iFED->SetFitelFiberSwitchBottomDauCard(0);

    /* read fifo1 */
    if( ReadFifo1 ){

      for( unsigned int ch = 0; ch < fedsAndChannels[ifed].second.size(); ch++ ){
       statusFifo1[ch] = iFED->drainFifo1(fedsAndChannels[ifed].second[ch], bufferFifo1[ch], 1024);
       statusFifo1Add[ch] = iFED->drainFifo1(fedsAndChannels[ifed].second[ch]+1, bufferFifo1Add[ch], 1024);
      }

      for( unsigned int ch = 0; ch < fedsAndChannels[ifed].second.size(); ch++ ){

       int channel = (fedsAndChannels[ifed].second)[ch];
       bool found_TBMA_H = false;
       bool found_TBMA_T = false;
       int addChannel = (fedsAndChannels[ifed].second)[ch]+1;
       bool found_TBMB_H = false;
       bool found_TBMB_T = false;
       std::vector<int> ch_decodedROCs;
       std::vector<int> addCh_decodedROCs;
       bool ch_foundWrongHit = false;
       bool addCh_foundWrongHit = false;
       bool ch_foundHit = false;
       bool addCh_foundHit = false;

       if( DumpFIFOs ){
       std::cout << "-----------------------------------" << std::endl;
       std::cout << "Contents of FIFO 1 for channel " << channel << " (status = " << statusFifo1[ch] << ")" << std::endl;
       std::cout << "-----------------------------------" << std::endl;
       }

       if (statusFifo1[ch] > 0) {

        for (int i = 0; i < statusFifo1[ch]; ++i) {
	      
         uint32_t w = bufferFifo1[ch][i];
         uint32_t decodeCh = (w >> 26) & 0x3f;
         if( decodeCh != channel ) continue;
         if( found_TBMA_H && found_TBMA_T ) continue;

         if( DumpFIFOs ){
	  std::cout << "Word " << std::setw(4) << std::setfill(' ') << i << " = 0x " << std::hex << std::setw(4) << std::setfill('0') << (bufferFifo1[ch][i]>>16) << " " << std::setw(4) << std::setfill('0') << (bufferFifo1[ch][i] & 0xFFFF) << std::dec << "  ";
	  for (int j = 31; j >= 0; --j){
	   if (w & (1 << j)) std::cout << "1";
	   else std::cout << "0";
	   if (j % 4 == 0) std::cout << " ";
	  }
	  std::cout << std::setfill(' ') << "  " << std::endl;
	  }
 
	 uint32_t mk = (w >> 21) & 0x1f;
	 uint32_t az = (w >> 8) & 0x1fff;
	 uint32_t dc = (w >> 16) & 0x1f;
	 uint32_t px = (w >> 8) & 0xff;
	 uint32_t f8 = w & 0xff;

	 if( DumpFIFOs ) std::cout << "  Decoded channel: " << decodeCh << std::endl;

	 if (!found_TBMA_H && mk == 0x1f) {
          if( DumpFIFOs ) printf("  TBM_H_status:%4x event# %i\n",((w>>1)&0xff00)+(w&0xff),f8);
          found_TBMA_H = true;
	 }
	 else if (!found_TBMA_T && mk == 0x1e) {
           if( DumpFIFOs ){
            printf("  TBM_T_status:%4x\n",((w>>4)&0xff00)+(w&0xff)); 
            if( (w&0x00000080)== 0x00000080 ) std::cout << "  ROCs disabled " << std::endl; 
            else if( (w&0x00000080)== 0 ) std::cout << "  ROCs enabled " << std::endl; 
           }
           found_TBMA_T = true;
	 }
	 else { 
	   if (az == 0 && !inject_ ){
            if( DumpFIFOs && mk <= 8 ) std::cout << "  ROC #" << std::dec << mk << "  lastdac: " << std::hex << f8 << std::dec << std::endl;
            if( mk <= 8 ) ch_decodedROCs.push_back(mk);
           }
           else if( az!=0 ){
            int column = dc*2 + px%2;
            int row = 80 - (px/2);
            if( DumpFIFOs && inject_ && mk <= 8 && column < 53 && row < 81 ){
              std::cout << "  FOUND HIT : ROC # " << std::dec << mk << " dcol " << dc << " pxl " << px;
              std::cout << " col " << column << " row " << row << " pulse height " << f8 << std::endl;
            }
            if( !inject_ ) ch_foundWrongHit = true;
            else{
             if( column == 25 && row == 25 && mk<=8 ){
              ch_decodedROCs.push_back(mk);
              //ch_foundHit = true;
             }
            }
           }
	 }
         if( DumpFIFOs ) std::cout << std::endl;
        }//close loop on status lenght

       }//close if status > 0

       //dump fifo for consecutive channel
       if( channel == 8 ) continue;

       if( DumpFIFOs ){
        std::cout << "-----------------------------------" << std::endl;
        std::cout << "Contents of FIFO 1 for channel " << addChannel << " (status = " << statusFifo1[ch+1] << ")" << std::endl;
        std::cout << "-----------------------------------" << std::endl;
       }
       if (statusFifo1Add[ch] > 0) {

        for (int i = 0; i < statusFifo1Add[ch]; ++i) {
	      
         uint32_t w = bufferFifo1Add[ch][i];
         uint32_t decodeAddCh = (w >> 26) & 0x3f;
         if( decodeAddCh != addChannel ) continue;
         if( found_TBMB_H && found_TBMB_T ) continue;

         if( DumpFIFOs ){
	  std::cout << "Word " << std::setw(4) << std::setfill(' ') << i << " = 0x " << std::hex << std::setw(4) << std::setfill('0') << (bufferFifo1Add[ch][i]>>16) << " " << std::setw(4) << std::setfill('0') << (bufferFifo1Add[ch][i] & 0xFFFF) << std::dec << "  ";
	  for (int j = 31; j >= 0; --j){
	   if (w & (1 << j)) std::cout << "1";
	   else std::cout << "0";
	   if (j % 4 == 0) std::cout << " ";
	  }
	  std::cout << std::setfill(' ') << "  " << std::endl;
	 }

	 uint32_t mk = (w >> 21) & 0x1f;
	 uint32_t az = (w >> 8) & 0x1fff;
	 uint32_t dc = (w >> 16) & 0x1f;
	 uint32_t px = (w >> 8) & 0xff;
	 uint32_t f8 = w & 0xff;

	 if( DumpFIFOs ) std::cout << "  Decoded channel: " << decodeAddCh << std::endl;

	 if (!found_TBMB_H && mk == 0x1f) {
          if( DumpFIFOs ) printf("  TBM_H_status:%4x event# %i\n",((w>>1)&0xff00)+(w&0xff),f8);
          found_TBMB_H = true;
	 }
	 else if (!found_TBMB_T && mk == 0x1e) {
           if( DumpFIFOs ){ 
            printf("  TBM_T_status:%4x\n",((w>>4)&0xff00)+(w&0xff)); 
            if( (w&0x00000080)== 0x00000080 ) std::cout << "  ROCs disabled " << std::endl; 
            else if( (w&0x00000080)== 0 ) std::cout << "  ROCs enabled " << std::endl; 
           }
           found_TBMB_T = true;
	 }
	 else { 
	   if (az == 0 && !inject_){
            if( DumpFIFOs ) std::cout << "  ROC #" << std::dec << mk << "  lastdac: " << std::hex << f8 << std::dec << std::endl;
            if (mk <= 8) addCh_decodedROCs.push_back(mk);
           }
           else if( az!=0 ){
            int column = dc*2 + px%2;
            int row = 80 - (px/2);
            if( DumpFIFOs && inject_ && mk <= 8 && column < 53 && row < 81 ){
              std::cout << "  FOUND HIT : ROC # " << std::dec << mk << " dcol " << dc << " pxl " << px;
              std::cout << " col " << column << " row " << row << " pulse height " << f8 << std::endl;
            }
            if( !inject_ ) addCh_foundWrongHit = true;
            else{
             if( column == 25 && row == 25 && mk<=8 ){
              addCh_decodedROCs.push_back(mk);
              //addCh_foundHit = true;
             }
            }
           }
	 }
         if( DumpFIFOs ) std::cout << std::endl;
        }//close loop on status lenght

       }//close if status > 0

       addCh_foundHit = (addCh_decodedROCs.size() == 4);
       ch_foundHit = (ch_decodedROCs.size() == 4);

       //std::cout << "**************************** " << addCh_decodedROCs.size() << " " << ch_decodedROCs.size() << " ";
       //std::cout << addCh_foundHit << " " << ch_foundHit << std::endl; 

       FillEm(state, fedsAndChannels[ifed].first, channel, 0, (!inject_ && found_TBMA_H && found_TBMA_T) || (inject_ && found_TBMA_H && found_TBMA_T && ch_foundHit) );
       FillEm(state, fedsAndChannels[ifed].first, addChannel, 0, (!inject_ && found_TBMB_H && found_TBMB_T) || (inject_ && found_TBMB_H && found_TBMB_T && addCh_foundHit) );

       if( dacsToScan.size() == 0 ){
        for( int r = 0; r < 8; ++r ){

         bool ch_foundROC = std::find(ch_decodedROCs.begin(),ch_decodedROCs.end(),r+1)!=ch_decodedROCs.end();
         bool addCh_foundROC = std::find(addCh_decodedROCs.begin(),addCh_decodedROCs.end(),r+1)!=addCh_decodedROCs.end();
         if( ch_foundROC ) FillEm(state, fedsAndChannels[ifed].first, channel, 1, r);
         if( addCh_foundROC ) FillEm(state, fedsAndChannels[ifed].first, addChannel, 1, r);
        }
       }
       else if( dacsToScan.size() == 1 && tempCalibObject->scanName(0) == "TBMPLL"){
        FillEm(state, fedsAndChannels[ifed].first, channel, 1, ch_decodedROCs.size());
        FillEm(state, fedsAndChannels[ifed].first, addChannel, 1, addCh_decodedROCs.size());
       }
       else if( dacsToScan.size() == 1 && tempCalibObject->scanName(0) == "TBMADelay" ){

        for( int r = 0; r < 8; ++r ){

         bool ch_foundROC = std::find(ch_decodedROCs.begin(),ch_decodedROCs.end(),r+1)!=ch_decodedROCs.end();
         bool addCh_foundROC = std::find(addCh_decodedROCs.begin(),addCh_decodedROCs.end(),r+1)!=addCh_decodedROCs.end();
         if( ch_foundROC ) FillEm(state, fedsAndChannels[ifed].first, channel, 0, 1, r);
         if( addCh_foundROC ) FillEm(state, fedsAndChannels[ifed].first, addChannel, 1, 1, r);

        }
       }

      }// end loop on channels

    }//end readFifo1

    /* read transparent fifo*/
    if( ReadTransFifo ){
     for( unsigned int ch = 0; ch < fedsAndChannels[ifed].second.size(); ch++ ){
      iFED->SelectTransparnetChannel(fedsAndChannels[ifed].second[ch]);
      iFED->drainDigTransFifoByChannel(fedsAndChannels[ifed].second[ch], bufferT[ch]);
      decodeT[ch] = new DigTransDecoder(bufferT[ch]);
     }
     //fill histos
     for( unsigned int ch = 0; ch < fedsAndChannels[ifed].second.size(); ch++ ){

      DigTransDecoder* d = decodeT[ch];
      int channel = (fedsAndChannels[ifed].second)[ch];

      FillEm(state, fedsAndChannels[ifed].first, channel, 0, bool(d->tbm_header_l[0].size()) && bool(d->tbm_trailer_l[0].size()) );
      FillEm(state, fedsAndChannels[ifed].first, channel+1, 0, bool(d->tbm_header_l[1].size()) && bool(d->tbm_trailer_l[1].size()) );
      FillEm(state, fedsAndChannels[ifed].first, channel, 1, d->roc_header_l[0].size());
      FillEm(state, fedsAndChannels[ifed].first, channel+1, 1, d->roc_header_l[1].size());

     }
     //print content
     if (DumpFIFOs) {
      for( unsigned int ch = 0; ch < fedsAndChannels[ifed].second.size(); ch++ ){

       std::cout << "-----------------------------------------\n";
       std::cout << "Contents of transparent FIFO for FED " << fednumber << " channel " << (fedsAndChannels[ifed].second)[ch] << std::endl;
       std::cout << "-----------------------------------------\n";

       std::cout << "** TBM-A ** " << std::endl;  
       std::cout << " number of tbm headers = " << decodeT[ch]->tbm_header_l[0].size() << std::endl;
       std::cout << " number of tbm trailers = " << decodeT[ch]->tbm_trailer_l[0].size() << std::endl;
       std::cout << " number of roc headers = " << decodeT[ch]->roc_header_l[0].size() << std::endl;
       std::cout << "** TBM-B ** " << std::endl;  
       std::cout << " number of tbm headers = " << decodeT[ch]->tbm_header_l[1].size() << std::endl;
       std::cout << " number of tbm trailers = " << decodeT[ch]->tbm_trailer_l[1].size() << std::endl;
       std::cout << " number of roc headers = " << decodeT[ch]->roc_header_l[1].size() << std::endl;

       std::cout << "---" << std::endl;
       std::cout << "DigTransDecoder thinks:\n";
       std::cout << "---" << std::endl;
       decodeT[ch]->printToStream(std::cout);
       std::cout << "---" << std::endl;
       for( int i = 0; i < 42; ++i ){

        unsigned long int tmp  = bufferT[ch][i];
        printf("%04x %04x              ",int((tmp>>16)&0xffff), int(tmp&0xffff)); //show in hex
        BinDisplay((tmp>>16)&0xffff, 16,0,0); //show in binary
        printf("    ");
        BinDisplay(tmp&0xffff, 16,0,0); //show in binary
        printf("\n");

       }//end read buffer for trans fifo for a given channel
      }// end loop on channel     
     }//end dumping trans fifo
    }// end ReadTransFifo case

    /* read scope fifo */
    if( ReadScopeFifo ){
     //read
     for( unsigned int ch = 0; ch < fedsAndChannels[ifed].second.size(); ch++ ){     
      iFED->SelectScopeChannel(fedsAndChannels[ifed].second[ch]);
      statusS[ch] = iFED->drainDataChannelFifo2(fedsAndChannels[ifed].second[ch], bufferS[ch]);
      decodeS[ch] = new DigScopeDecoder(bufferS[ch], statusS[ch]);
     }
     //fill histos
     for( unsigned int ch = 0; ch < fedsAndChannels[ifed].second.size(); ch++ ){

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
    }//end ReadScopeFifo case

    /*read tmp fifo*/
    if( ReadTmpFifo ){ 
     
     uint32_t bufferTemp[8][256];
     iFED->ReadTemp(bufferTemp[0],1);
     iFED->ReadTemp(bufferTemp[1],5);
     iFED->ReadTemp(bufferTemp[2],10);
     iFED->ReadTemp(bufferTemp[3],14);
     iFED->ReadTemp(bufferTemp[4],19);
     iFED->ReadTemp(bufferTemp[5],23);
     iFED->ReadTemp(bufferTemp[6],28);
     iFED->ReadTemp(bufferTemp[7],32);
     for( unsigned int ch = 0; ch < fedsAndChannels[ifed].second.size(); ch++ ){

      int channel = (fedsAndChannels[ifed].second)[ch];
      if( channel == 8 ) continue;
      int addChannel = (fedsAndChannels[ifed].second)[ch]+1;
      int chip = channel/5; if( channel > 27 && chip < 7 ) chip=chip+1;
      int addChip = (addChannel+1)/5; if( addChannel > 27 && addChip < 7 ) addChip=addChip+1;

      if(DumpFIFOs){
       std::cout << "----------------------------------" << std::endl;
       std::cout << "Contents of Temp FIFO for FED " << fednumber << " channels " << channel << "/" << addChannel;
       std::cout << " chip " << chip << "/" << addChip << std::endl;
       std::cout << "----------------------------------" << std::endl;
      }

      bool foundCh = false;
      bool foundAddCh = false;
      bool found_TBMA_H = false;
      bool found_TBMA_T = false;
      bool found_TBMB_H = false;
      bool found_TBMB_T = false;

      for( int i=0; i < 256; ++i ){
       uint32_t d = bufferTemp[chip][i];
       bool foundChTmp = false;
       bool foundAddChTmp = false;

       if(d){
        int decodeCh = (d>>26)&0x3f;
        if(decodeCh == channel){ foundCh = true; foundChTmp = true; }
        if(foundCh && ((d>>21)&0x1f)== 31 && (d&0x000000ff)==event_+1) found_TBMA_H = true;
        if(foundCh && ((d>>21)&0x1f)== 30) found_TBMA_T = true;
       
        if(DumpFIFOs && foundChTmp){
         printf("check %x\n",d); 
         printf("CH#:%2d\n",((d>>26)&0x3f));
         if(((d>>21)&0x1f)== 31) printf("  TBM_H_status:%4x event# %i\n",((d>>1)&0xff00)+(d&0xff),(d&0x000000ff));
         if(((d>>21)&0x1f)== 30) printf("  TBM_T_status:%4x\n",((d>>4)&0xff00)+(d&0xff)); 
        }//close printing fifo content
       }

       uint32_t d1 = bufferTemp[addChip][i];
       if(d1){
        int decodeAddCh = (d1>>26)&0x3f;
        if(decodeAddCh == addChannel){ foundAddCh = true; foundAddChTmp = true; }
        if(foundAddCh && ((d1>>21)&0x1f)== 31 && (d1&0x000000ff)==event_+1) found_TBMB_H = true;
        if(foundAddCh && ((d1>>21)&0x1f)== 30) found_TBMB_T = true;

        if(DumpFIFOs && foundAddChTmp){
         printf("check %x\n",d1); 
         printf("CH#:%2d\n",((d1>>26)&0x3f));
         if(((d1>>21)&0x1f)== 31) printf("  TBM_H_status:%4x event# %i\n",((d1>>1)&0xff00)+(d1&0xff),(d1&0x000000ff));
         if(((d1>>21)&0x1f)== 30) printf("  TBM_T_status:%4x\n",((d1>>4)&0xff00)+(d1&0xff)); 
        }//close printing fifo content
       }

      }//close loop on buffer for given channel

      if(DumpFIFOs){std::cout << "----------------------------------" << std::endl;}

      FillEm(state, fedsAndChannels[ifed].first, channel, 0, found_TBMA_H && found_TBMA_T && foundCh );
      FillEm(state, fedsAndChannels[ifed].first, channel+1, 0, found_TBMB_H && found_TBMB_T && foundAddCh );

     }//close loop on channels

    }//close readtmpfifo case
    
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
  if (dacsToScan.size() == 0){

    for( std::map<int,std::map<int,std::vector<TH1F*> > >::iterator it1 = ntrigsTBM.begin(); it1 != ntrigsTBM.end(); ++it1 ){
     for( std::map<int,std::vector<TH1F*> >::iterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2 ){
      for(unsigned int i = 0; i < it2->second.size(); ++i ) it2->second[i]->Scale(1./ntriggers);
     }
    }

  }

  if (dacsToScan.size() == 1){

   for( std::map<int,std::map<int,std::vector<TH2F*> > >::iterator it1 = scansTBM.begin(); it1 != scansTBM.end(); ++it1 ){
    for( std::map<int,std::vector<TH2F*> >::iterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2 ){
     for(unsigned int i = 0; i < it2->second.size(); ++i ) it2->second[i]->Scale(1./ntriggers);
    }
   }

   /*for( std::map<int,std::map<int,std::vector<TH2F*> > >::iterator it1 = scansTBM.begin(); it1 != scansROCs.end(); ++it1 ){
    for( std::map<int,std::vector<TH2F*> >::iterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2 ){
     for(unsigned int i = 0; i < it2->second.size(); ++i ) it2->second[i]->Scale(1./ntriggers);
    }
   }*/

  }

  CloseRootf();
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
  scansROCs.clear();

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

     if(tempCalibObject->scanName(0) == "TBMPLL"){
      TString hname; hname.Form("Ch%i",(fedsAndChannels[ifed].second)[ch]);
      std::vector<TH2F*> histosTBM;
      TH2F* h_TBM_nDecodes = new TH2F(hname+"_TBM_nDecodes", hname+"_TBM_nDecodes", 8, 0, 8, 8, 0, 8 );
      histosTBM.push_back(h_TBM_nDecodes);
      TH2F* h_nROCHeaders = new TH2F(hname+"_nROCHeaders", hname+"_nROCHeaders", 8, 0, 8, 8, 0, 8 );
      histosTBM.push_back(h_nROCHeaders);     
      chTBMmap2D[(fedsAndChannels[ifed].second)[ch]] = histosTBM;
     }
  
    }// end book histos for 1 dacsToScan case (TBMPLL scan)

    if(dacsToScan.size() == 1){

     if(tempCalibObject->scanName(0) == "TBMADelay"){

      std::vector<TH2F*> histosROCs;
      for( int r = 0; r < 8; ++r ){
       TString hname; hname.Form("Ch%i_ROC%i",(fedsAndChannels[ifed].second)[ch],r);
       TH2F* h_nROCHeaders = new TH2F(hname+"_nHeaders", hname+"_nHeaders", 8, 0, 8, 8, 0, 8 );
       histosROCs.push_back(h_nROCHeaders);    
      }
      chTBMmap2D[(fedsAndChannels[ifed].second)[ch]] = histosROCs;
     }
  
    }// end book histos for 2 dacsToScan case (TBMA/B delay scan)

    if( ReadScopeFifo ) continue;
    //add the consecutive channel until Danek fixes the name translation
    dir->cd();
    int addCh = (fedsAndChannels[ifed].second)[ch]+1;
    chDir.Form("Channel%i",addCh);
    dir->mkdir(chDir);
    dir->cd(chDir);

    if (dacsToScan.size() == 0){

     TString hname; hname.Form("Ch%i",addCh);
     std::vector<TH1F*> histosTBM;
     TH1F* h_TBM_nDecodes = new TH1F(hname+"_TBM_nDecodes", hname+"_TBM_nDecodes", 2, 0, 2 );
     histosTBM.push_back(h_TBM_nDecodes);
     TH1F* h_nROCHeaders = new TH1F(hname+"_nROCHeaders", hname+"_nROCHeaders", 9, 0, 9 );
     histosTBM.push_back(h_nROCHeaders);     
     chTBMmap[addCh] = histosTBM;
  
    }// end book histos for zero dacsToScan case

    if(dacsToScan.size() == 1 && tempCalibObject->scanName(0) == "TBMPLL" ){

     TString hname; hname.Form("Ch%i",addCh);
     std::vector<TH2F*> histosTBM;
     TH2F* h_TBM_nDecodes = new TH2F(hname+"_TBM_nDecodes", hname+"_TBM_nDecodes", 8, 0, 8, 8, 0, 8 );
     histosTBM.push_back(h_TBM_nDecodes);
     TH2F* h_nROCHeaders = new TH2F(hname+"_nROCHeaders", hname+"_nROCHeaders", 8, 0, 8, 8, 0, 8 );
     histosTBM.push_back(h_nROCHeaders);     
     chTBMmap2D[addCh] = histosTBM;
  
    }// end book histos for 1 dacsToScan case

    if(dacsToScan.size() == 1){

     if(tempCalibObject->scanName(0) == "TBMADelay"){
      std::vector<TH2F*> histosROCs;
      for( int r = 0; r < 8; ++r ){
       TString hname; hname.Form("Ch%i_ROC%i",addCh,r);
       TH2F* h_nROCHeaders = new TH2F(hname+"_nHeaders", hname+"_nHeaders", 8, 0, 8, 8, 0, 8 );
       histosROCs.push_back(h_nROCHeaders);    
      } 
      chTBMmap2D[addCh] = histosROCs;
     }
  
    }// end book histos for 2 dacsToScan case (TBMA/B delay scan)

   }//close loop on channels

   if (dacsToScan.size() == 0) ntrigsTBM[fedsAndChannels[ifed].first] = chTBMmap;
   if (dacsToScan.size() == 1 && tempCalibObject->scanName(0) == "TBMPLL") scansTBM[fedsAndChannels[ifed].first] = chTBMmap2D;
   if (dacsToScan.size() == 1 && tempCalibObject->scanName(0) == "TBMADelay") scansROCs[fedsAndChannels[ifed].first] = chTBMmap2D;

  }//close loop on feds

  rootf->cd(0);

}

///////////////////////////////////////////////////////////////////////////////////////////////
void PixelFEDTBMDelayCalibration::FillEm(unsigned state, int fedid, int ch, int which, float c, int roc) {
  PixelCalibConfiguration* tempCalibObject = dynamic_cast<PixelCalibConfiguration*>(theCalibObject_);
  assert(tempCalibObject != 0);

  if (event_==0) return;

  if (dacsToScan.size() == 1 && tempCalibObject->scanName(0) == "TBMADelay"){

   const std::string& iname = dacsToScan[0];
   const double ival(tempCalibObject->scanValue(iname, state)); 
   uint32_t tmp = ival; 
   int delay1 = (tmp>>3)&0x7;
   int delay2 = tmp&0x7;
   scansROCs[fedid][ch][roc]->Fill(delay1,delay2,c);

  }

}

///////////////////////////////////////////////////////////////////////////////////////////////
void PixelFEDTBMDelayCalibration::FillEm(unsigned state, int fedid, int ch, int which, float c) {
  PixelCalibConfiguration* tempCalibObject = dynamic_cast<PixelCalibConfiguration*>(theCalibObject_);
  assert(tempCalibObject != 0);

  if (event_==0) return;

  if (dacsToScan.size() == 0 ) ntrigsTBM[fedid][ch][which]->Fill(c,1);
  if (dacsToScan.size() == 1 && tempCalibObject->scanName(0) == "TBMPLL"){

   const std::string& iname = dacsToScan[0];
   const double ival(tempCalibObject->scanValue(iname, state)); 
   uint32_t tmp = ival; 
   int delay1 = (tmp>>2)&0x7;
   int delay2 = ((tmp>>2)&0x38)>>3;
   scansTBM[fedid][ch][which]->Fill(delay1,delay2,c);

  }

}
