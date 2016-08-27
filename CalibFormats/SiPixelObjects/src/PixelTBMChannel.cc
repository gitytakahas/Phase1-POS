#include "CalibFormats/SiPixelObjects/interface/PixelTBMChannel.h"
#include <cassert>

using namespace pos;

PixelTBMChannel::PixelTBMChannel(std::string TBMChannel)
{
   std::cout<<TBMChannel<<std::endl;

  if        ( TBMChannel=="A" ) {
    isChannelB_ = false;
    channel_=0;
  } else if ( TBMChannel=="B" ) {
    isChannelB_ = true;
    channel_=1;
  } else if ( TBMChannel=="A1" ) {
    isChannelB_ = true;
    channel_=10;
  } else if ( TBMChannel=="A2" ) {
    isChannelB_ = true;
    channel_=11;
  } else if ( TBMChannel=="B1" ) {
    isChannelB_ = true;
    channel_=12;
  } else if ( TBMChannel=="B2" ) {
    isChannelB_ = true;
    channel_=13;
  } else if ( TBMChannel=="C1" ) {
    isChannelB_ = true;
    channel_=14;
  } else if ( TBMChannel=="C2" ) {
    isChannelB_ = true;
    channel_=15;
  } else if ( TBMChannel=="D1" ) {
    isChannelB_ = true;
    channel_=16;
  } else if ( TBMChannel=="D2" ) {
    isChannelB_ = true;
    channel_=17;
  } else {
    std::cout << "ERROR in PixelTBMChannel: TBM channel must be A or B, but input value was "<<TBMChannel<<std::endl;
    assert(0);
  }
  std::cout<<"Ch1 "<<TBMChannel<<" "<<channel_<<" "<<isChannelB_<<std::endl;

}

std::string PixelTBMChannel::string() const {
  if(channel_==0) return "A";
  else if(channel_==1) return "B";
  else if(channel_==10) return "A";
  else if(channel_==11) return "A";
  else if(channel_==12) return "B";
  else if(channel_==13) return "B";
  else if(channel_==14) return "A";
  else if(channel_==15) return "A";
  else if(channel_==16) return "B";
  else if(channel_==17) return "B ";
  else std::cout<<" something wrong in PixelTBMCHannel::string() "<<channel_<<std::endl;

  //if ( isChannelB_ ) return "B";
  //else              return "A";
}

std::string PixelTBMChannel::stringFull() const {
  if(channel_==0) return "A";
  else if(channel_==1) return "B";
  else if(channel_==10) return "A1";
  else if(channel_==11) return "A2";
  else if(channel_==12) return "B1";
  else if(channel_==13) return "B2";
  else if(channel_==14) return "C1";
  else if(channel_==15) return "C2";
  else if(channel_==16) return "D1";
  else if(channel_==17) return "D2";
  else std::cout<<" something wrong in PixelTBMCHannel::string() "<<channel_<<std::endl;

  //if ( isChannelB_ ) return "B";
  //else              return "A";
}

std::ostream& pos::operator<<(std::ostream& s, const PixelTBMChannel& TBMChannel)
{
	s << TBMChannel.stringFull();
	return s;
}
