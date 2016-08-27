
#ifndef _DigFIFO1Decoder_h_
#define _DigFIFO1Decoder_h_

#include <vector>
#include <stdint.h>
#include <iostream>

class DigFIFO1Decoder
{
 public:
    
    DigFIFO1Decoder(uint32_t *buffer, int nWords);
    void printToStream(std::ostream& out);
    void DecodeError(uint32_t w);

    unsigned int nhits() const {return hits_.size();}
    unsigned int nerrs() const {return errors_;}

    unsigned int channel(unsigned int ihit) const {return (hits_[ihit]>>26)&0x3f;}
    unsigned int rocid(unsigned int ihit) const {return (hits_[ihit]>>21)&0x1f;}
    
    unsigned int dcol(unsigned int ihit) const {return (hits_[ihit]>>16)&0x1f;}
    unsigned int pxl(unsigned int ihit) const {return (hits_[ihit]>>8)&0xff;}
    unsigned int pulseheight(unsigned int ihit) const {return hits_[ihit]&0xff;}
    
    unsigned int column(unsigned int ihit) const {return dcol(ihit)*2 + pxl(ihit)%2;}
    unsigned int row(unsigned int ihit) const {return 80 - (pxl(ihit)/2);}

    unsigned int globalChannel() const {return channel_;}
    int event() const {return eventNumber_;}
    bool foundTBM() const {return (found_tbm_header_ && found_tbm_trailer_);}
    std::vector<int> ROCHeaders() const {return roc_headers_;}
    std::vector<int> ROCLastDac() const {return last_dacs_;}

 private:
    std::vector<int> roc_headers_;
    std::vector<int> last_dacs_;
    std::vector<int> rocsWithPixelHits_;
    std::vector<uint32_t> hits_;
    int errors_;
    bool found_tbm_header_;
    bool found_tbm_trailer_;
    unsigned int channel_;
    int eventNumber_;
    bool rocsEnabled_;

};

#endif
