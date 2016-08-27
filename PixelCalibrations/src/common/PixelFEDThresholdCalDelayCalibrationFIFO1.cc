// Modified by Jennifer Vaughan 2007/06/01
// $Id: PixelThresholdCalDelayCalibrationFIFO1.cc,v 1.1 

/*************************************************************************
 * XDAQ Components for Distributed Data Acquisition                      *
 * Copyright (C) 2000-2004, CERN.			                 *
 * All rights reserved.                                                  *
 * Authors: J. Gutleber and L. Orsini					 *
 *                                                                       *
 * For the licensing terms see LICENSE.		                         *
 * For the list of contributors see CREDITS.   			         *
 *************************************************************************/

#include "PixelCalibrations/include/PixelFEDThresholdCalDelayCalibrationFIFO1.h"
#include "PixelUtilities/PixelFEDDataTools/include/Moments.h"
#include "PixelUtilities/PixelFEDDataTools/include/PixelFEDDataTypes.h"
#include "PixelUtilities/PixelFEDDataTools/include/PixelDecodedFEDRawData.h"
#include "PixelCalibrations/include/PixelEfficiency2DVcThrCalDel.h"
#include "PixelUtilities/PixelFEDDataTools/include/FIFO1Decoder.h"
#include "PixelUtilities/PixelFEDDataTools/include/FIFO2Decoder.h"
#include "PixelUtilities/PixelFEDDataTools/include/FIFO3Decoder.h"
#include "PixelUtilities/PixelFEDDataTools/include/ErrorFIFODecoder.h"
#include "PixelConfigDBInterface/include/PixelConfigInterface.h"
#include "CalibFormats/SiPixelObjects/interface/PixelCalibConfiguration.h"
#include "PixelFEDInterface/include/PixelFEDFifoData.h"

#include <toolbox/convertstring.h>

#include "iomanip"

using namespace pos;


PixelFEDThresholdCalDelayCalibrationFIFO1::PixelFEDThresholdCalDelayCalibrationFIFO1(const PixelFEDSupervisorConfiguration & tempConfiguration, SOAPCommander* mySOAPCmdr) 
  : PixelFEDCalibrationBase(tempConfiguration,*mySOAPCmdr)
{
  std::cout << "Greetings from the PixelFEDThresholdCalDelayCalibrationFIFO1 copy constructor." << std::endl;
}

xoap::MessageReference PixelFEDThresholdCalDelayCalibrationFIFO1::execute(xoap::MessageReference msg)
{

  Attribute_Vector parameters(2);

  parameters[0].name_="Event";
  parameters[1].name_="VMEBaseAddress"; //FIXME should not need the base add.
  Receive(msg, parameters);

  unsigned long vmeBaseAddress = atoi(parameters[1].value_.c_str());
  unsigned int event = atoi(parameters[0].value_.c_str());

  const PixelCalibConfiguration* tempCalibObject=dynamic_cast<PixelCalibConfiguration*>(theCalibObject_);
  assert(tempCalibObject!=0);

  unsigned int state=event/tempCalibObject->nTriggersPerPattern();

  unsigned int ithreshold=tempCalibObject->scanCounter("VcThr",state);
  unsigned int icaldelay=tempCalibObject->scanCounter("CalDel",state);

  static std::map <PixelROCName, PixelEfficiency2DVcThrCalDel> eff;

  if (event==0) {

    const std::vector<PixelROCName> aROC_string=tempCalibObject->rocList();

    assert(aROC_string.size()>0);

    PixelModuleName module(aROC_string[0].rocname());

    if (crate_== theFEDConfiguration_->crateFromFEDNumber(theNameTranslation_->firstHdwAddress(module).fednumber())  ) {


      unsigned int nThr=tempCalibObject->nScanPoints("VcThr");
      unsigned int nCal=tempCalibObject->nScanPoints("CalDel");

      double VcThrMin=tempCalibObject->scanValueMin("VcThr");
      double VcThrMax=tempCalibObject->scanValueMax("VcThr");
      double VcThrStep=tempCalibObject->scanValueStep("VcThr");
      
      double CalDelMin=tempCalibObject->scanValueMin("CalDel");
      double CalDelMax=tempCalibObject->scanValueMax("CalDel");
      double CalDelStep=tempCalibObject->scanValueStep("CalDel");


      for (unsigned int i_aROC=0;i_aROC<aROC_string.size();++i_aROC) {

	//std::cout << "Creating ROC:"<<aROC_string[i_aROC].rocname()<<std::endl;

	PixelEfficiency2DVcThrCalDel tmp(aROC_string[i_aROC].rocname(),"CalDel",nCal,
					 CalDelMin-0.5*CalDelStep,
					 CalDelMax+0.5*CalDelStep,
					 "VcThr",nThr,
					 VcThrMin-0.5*VcThrStep,
					 VcThrMax+0.5*VcThrStep);
	//std::cout << "Made temporary efficiency" << std::endl;
	eff[aROC_string[i_aROC]]=tmp;
	//std::cout << "Copied efficiency" << std::endl;
      }
    }
    //std::cout << "Done with initialization"<<std::endl;
  }

  
  std::vector<std::pair<unsigned int,std::vector<unsigned int> > > fedsAndChannels=tempCalibObject->fedCardsAndChannels(crate_, theNameTranslation_, theFEDConfiguration_, theDetectorConfiguration_ );
  
  for (unsigned int ifed=0; ifed<fedsAndChannels.size(); ++ifed) {
    unsigned int fednumber=fedsAndChannels[ifed].first;
    //FIXME should look at all FEDs at once...
    if (vmeBaseAddress!=theFEDConfiguration_->VMEBaseAddressFromFEDNumber(fednumber)) continue;
    
    for (unsigned int ichannel=0; ichannel<fedsAndChannels[ifed].second.size(); ++ichannel) {
      uint32_t buffer[pos::fifo1TranspDepth];
      unsigned int channel=fedsAndChannels[ifed].second[ichannel];
      
      int status = FEDInterface_[vmeBaseAddress]->drainFifo1(channel, buffer);//drain_transBuffer(channel, buffer);

      /*if (status!=(int)pos::fifo1TranspDataLength) {
        std::cout<<"PixelFEDThresholdCalDelayCalibrationFIFO1::execute status="
                 <<status<<std::endl;
        std::cout<<"PixelFEDThresholdCalDelayCalibrationFIFO1::execute -- Could not drain FIFO 1 of FED Channel "<<channel<<" in transparent mode!"<<std::endl;
        diagService_->reportError("PixelFEDThresholdCalDelayCalibrationFIFO1::execute -- Could not drain FIFO 1 in transparent mode!",DIAGWARN);
      }*/		    

       bool found_TBMA_H = false;
       bool found_TBMA_T = false;
       std::map<int,int> ch_decodedROCs;
       bool ch_foundWrongHit = false;
       bool ch_foundHit = false;
       bool inject_ = true;

       if (status > 0) {

        for (int i = 0; i < status; ++i) {
	      
         uint32_t w = buffer[i];
         uint32_t decodeCh = (w >> 26) & 0x3f;
         //if( decodeCh != channel ) continue;
         //if( found_TBMA_H && found_TBMA_T ) continue;

         /*if( true ){
	  std::cout << "Word " << std::setw(4) << std::setfill(' ') << i << " = 0x " << std::hex << std::setw(4) << std::setfill('0') << (buffer[i]>>16) << " " << std::setw(4) << std::setfill('0') << (buffer[i] & 0xFFFF) << std::dec << "  ";
	  for (int j = 31; j >= 0; --j){
	   if (w & (1 << j)) std::cout << "1";
	   else std::cout << "0";
	   if (j % 4 == 0) std::cout << " ";
	  }
	  std::cout << std::setfill(' ') << "  " << std::endl;
	  }*/

	 uint32_t mk = (w >> 21) & 0x1f;
	 uint32_t az = (w >> 8) & 0x1fff;
	 uint32_t dc = (w >> 16) & 0x1f;
	 uint32_t px = (w >> 8) & 0xff;
	 uint32_t f8 = w & 0xff;

	 //if( DumpFIFOs ) std::cout << "  Decoded channel: " << decodeCh << std::endl;

	 if (!found_TBMA_H && mk == 0x1f) {
          //if( DumpFIFOs ) printf("  TBM_H_status:%4x event# %i\n",((w>>1)&0xff00)+(w&0xff),f8);
          found_TBMA_H = true;
	 }
	 else if (!found_TBMA_T && mk == 0x1e) {
           /*if( DumpFIFOs ){
            printf("  TBM_T_status:%4x\n",((w>>4)&0xff00)+(w&0xff)); 
            if( (w&0x00000080)== 0x00000080 ) std::cout << "  ROCs disabled " << std::endl; 
            else if( (w&0x00000080)== 0 ) std::cout << "  ROCs enabled " << std::endl; 
           }*/
           found_TBMA_T = true;
	 }
	 else { 
	   if (az == 0){
            //if( DumpFIFOs && mk <= 8 ) std::cout << "  ROC #" << std::dec << mk << "  lastdac: " << std::hex << f8 << std::dec << std::endl;
            if( mk <= 8 ) ch_decodedROCs[mk-1] = 0;
           }
           if( az!=0 ){
            int column = dc*2 + px%2;
            int row = 80 - (px/2);
            //if( DumpFIFOs && inject_ && mk <= 8 && column < 53 && row < 81 ){
            //  std::cout << "  FOUND HIT : ROC # " << std::dec << mk << " dcol " << dc << " pxl " << px;
            //  std::cout << " col " << column << " row " << row << " pulse height " << f8 << std::endl;
            //}
            if( !inject_ ) ch_foundWrongHit = true;
            else{
             if( /*column == 25 && row == 25 &&*/ mk<=8 ){
              ch_decodedROCs[mk-1] = 1;
              //ch_foundHit = true;
             }
            }
           }
	 }
         //if( DumpFIFOs ) std::cout << std::endl;
        }//close loop on status lenght
 
        }//close loop on status lenght



      //PixelFEDFifoData::decodeNormalData(buffer);
      //PixelDecodedFEDRawData rawdata(buffer);
      
      unsigned int nROCs=4;//rawdata.numROCs();
      
      for(unsigned int iROC=0;iROC<nROCs;iROC++){
        unsigned int nHits=ch_decodedROCs[iROC];//rawdata.ROCOutput(iROC).numHits();

        PixelROCName roc=
          theNameTranslation_->ROCNameFromFEDChannelROC(fednumber,
                                                        channel,
                                                        iROC);
        
        std::map <PixelROCName, PixelEfficiency2DVcThrCalDel>::iterator it=eff.find(roc);
        if (it!=eff.end()) {
          if (nHits>0) it->second.add(icaldelay,ithreshold);
        }
        else{
          std::cout << "Could not find ROC with fednumber="<<fednumber
                    << " channel="<<channel<<" rocid="<<iROC<<std::endl;
        }
        
      }
    }
  }

  //Are we processing the last event?
  if (event+1==tempCalibObject->nTriggersPerPattern()*
      tempCalibObject->nConfigurations()){
    

    //First we need to get the DAC settings for the ROCs
    
    std::vector<PixelModuleName> modules=theDetectorConfiguration_->getModuleList();
    std::vector<PixelModuleName>::iterator module_name = modules.begin();
    
    std::map<PixelModuleName,PixelDACSettings*> theDACs;
    
    for (;module_name!=modules.end();++module_name){

      if (crate_!= theFEDConfiguration_->crateFromFEDNumber(theNameTranslation_->firstHdwAddress(*module_name).fednumber())  ) {
	continue;
      }

      PixelDACSettings *tempDACs=0;
      
      std::string modulePath=(module_name->modulename());
      
      //std::cout << "Reading DACs for module:"<<modulePath<<std::endl;
      
      PixelConfigInterface::get(tempDACs, "pixel/dac/"+modulePath, *theGlobalKey_);
      assert(tempDACs!=0);
      theDACs[*module_name]=tempDACs;
      
    }


    std::ofstream outsummary;
    
    outsummary.open((outputDir()+"/VcThrCalDelaySummary_"+itoa(crate_)+".txt").c_str());
    
    assert(outsummary.good());
    
    std::map<PixelROCName, PixelEfficiency2DVcThrCalDel>::iterator it=eff.begin();
    
    for(;it!=eff.end();++it){
      std::string filename=outputDir()+"/ThrCalDelScan_";
      filename+=it->first.rocname();
      filename+=".dat";
      //std::cout << "Will determine settings for ROC:"
      //          << it->first.rocname()<<std::endl;
      it->second.findSettings(tempCalibObject->nTriggersPerPattern()*
                              tempCalibObject->nPixelPatterns());

      PixelModuleName moduleName(it->first.rocname());

      //std::cout << "moduleName:"<<moduleName<<std::endl;
      
      std::map<PixelModuleName,PixelDACSettings*>::iterator dacs=theDACs.find(moduleName);
      
      assert(dacs!=theDACs.end());
      
      PixelROCDACSettings *rocDACs=dacs->second->getDACSettings(it->first);
      
      assert(rocDACs!=0);
      
      unsigned int oldVcThr=rocDACs->getVcThr();
      unsigned int oldCalDel=rocDACs->getCalDel();
      
      //std::cout << "Old settings: VcThr="<<oldVcThr<<" CalDel="
      //          << oldCalDel << std::endl;
      
      unsigned int newVcThr=it->second.getThreshold();
      unsigned int newCalDel=it->second.getCalDelay();

      //std::cout << "New settings: VcThr="<<newVcThr<<" CalDel="
      //          << newCalDel << std::endl;

      if (it->second.validSettings()){
	rocDACs->setVcThr(newVcThr);
	rocDACs->setCalDel(newCalDel);
      }
      else{
	cout << "Did not have valid settings for:"
	     << it->first.rocname() << endl;
      }

      it->second.setOldThreshold(oldVcThr);
      it->second.setOldCalDelay(oldCalDel);

      it->second.setTitle(it->first.rocname());

      it->second.saveEff(tempCalibObject->nTriggersPerPattern()*
                         tempCalibObject->nPixelPatterns(),
                         filename);

      outsummary << filename << std::endl;

    }

    outsummary.close();

    std::map<PixelModuleName,PixelDACSettings*>::iterator dacs=theDACs.begin();

    for(;dacs!=theDACs.end();++dacs){
      dacs->second->writeASCII(outputDir());
    }

  }

  xoap::MessageReference reply = MakeSOAPMessageReference("ThresholdCalDelayDone");
  return reply;
}

void PixelFEDThresholdCalDelayCalibrationFIFO1::initializeFED(){

  setFEDModeAndControlRegister(0x0,0x30019);

}

xoap::MessageReference PixelFEDThresholdCalDelayCalibrationFIFO1::beginCalibration(xoap::MessageReference msg){

  xoap::MessageReference reply = MakeSOAPMessageReference("BeginCalibrationDone");
  return reply;
}

xoap::MessageReference PixelFEDThresholdCalDelayCalibrationFIFO1::endCalibration(xoap::MessageReference msg){

  xoap::MessageReference reply = MakeSOAPMessageReference("EndCalibrationDone");
  return reply;
}

