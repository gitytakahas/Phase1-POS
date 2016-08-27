#include "PixelUtilities/PixelFEDDataTools/include/DigFIFO1Decoder.h"
#include <iostream>
#include <assert.h>
#include <iomanip>
#include <algorithm> 
#include "PixelUtilities/PixelFEDDataTools/include/dump_error_files_masks.h"

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
DigFIFO1Decoder::DigFIFO1Decoder(uint32_t *buffer, int nWords)
{

  found_tbm_header_ = false;
  found_tbm_trailer_ = false;
  eventNumber_ = -1;
  rocsEnabled_ = -1;
  errors_ = 0;
  roc_headers_.clear();
  rocsWithPixelHits_.clear();
  hits_.clear();

  for( int i = 0; i < nWords; ++i ){

   uint32_t w = buffer[i];

   channel_ = (w >> 26) & 0x3f;
   uint32_t mk = (w >> 21) & 0x1f;
   uint32_t az = (w >> 8) & 0x1fff;
   uint32_t dc = (w >> 16) & 0x1f;
   uint32_t px = (w >> 8) & 0xff;
   uint32_t f8 = w & 0xff;

         if( false ){
	  std::cout << "Word " << std::setw(4) << std::setfill(' ') << i << " = 0x " << std::hex << std::setw(4) << std::setfill('0') << (buffer[i]>>16) << " " << std::setw(4) << std::setfill('0') << (buffer[i] & 0xFFFF) << std::dec << "  ";
	  for (int j = 31; j >= 0; --j){
	   if (w & (1 << j)) std::cout << "1";
	   else std::cout << "0";
	   if (j % 4 == 0) std::cout << " ";
	  }
	  std::cout << std::setfill(' ') << "  " << std::endl;
	  }

   if (!found_tbm_header_ && mk == 0x1f) {
    eventNumber_ = f8;
    found_tbm_header_ = true;
   }
   else if (!found_tbm_trailer_ && mk == 0x1e) {
     if( (w&0x00000080)== 0x00000080 ) rocsEnabled_ = 1;
     else if( (w&0x00000080)== 0 ) rocsEnabled_ = 0; 
     found_tbm_trailer_ = true;
     DecodeError(w);
   }
   else { 
    if (az == 0 ){
     if( mk <= 8 && std::find( roc_headers_.begin(),roc_headers_.end(),mk)==roc_headers_.end() ){
      roc_headers_.push_back(mk); 
      last_dacs_.push_back(f8);
      }
    }
    else if( az!=0 ){
      int column = dc*2 + px%2;
      int row = 80 - (px/2);
      if( mk <= 8 && column < 53 && row < 81 ){
        if( std::find( rocsWithPixelHits_.begin(),rocsWithPixelHits_.end(),mk)==rocsWithPixelHits_.end() ) rocsWithPixelHits_.push_back(mk);
        hits_.push_back(w);
      }
    }//close hits
   }//close rocs

  }

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void DigFIFO1Decoder::DecodeError(uint32_t w){

 if ( ((w & masks::errorMask) == masks::overflowMask) || 
      ((w & masks::errorMask) == masks::RocErrMask)   || 
      ((w & masks::errorMask) == masks::FsmErrMask)  )
{ 
  
  errors_+=1;

 }

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void DigFIFO1Decoder::printToStream(std::ostream& out){

 out << "  Decoded channel: " << channel_ << "\n";
 out << "  Found TBM header: " << found_tbm_header_ << "\n";
 out << "  Found TBM trailer: " << found_tbm_trailer_ << "\n";
 out << "  Event number: " << eventNumber_ << "\n";
 out << "  ROCs enabled: " << rocsEnabled_ << "\n";
 out << "  TBM trailer errors: " << errors_ << "\n";
 out << "  Found ROCs without hits: " << roc_headers_.size() << "\n";
 for( int r = 0; r < roc_headers_.size(); ++r )
  out << "     - ROC #" << roc_headers_[r] << " - last dac : " << last_dacs_[r] << "\n";
 out << "  Found ROCs with hits: " << rocsWithPixelHits_.size() << "\n";
 out << "  Total number of hits found: " << nhits() << "\n";
 for( int h = 0; h < hits_.size(); ++h ){
  out << "     - Channel #" << channel(h) << " ROC #" << rocid(h) << " dcol " << dcol(h) << " pxl " << pxl(h);
  out << " col " << column(h) << " row " << row(h) << " pulse height " << pulseheight(h) << "\n";
 }

}
