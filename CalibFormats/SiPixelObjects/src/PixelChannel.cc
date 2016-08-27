#include "CalibFormats/SiPixelObjects/interface/PixelChannel.h"

using namespace pos;

PixelChannel::PixelChannel(PixelModuleName module, std::string TBMChannel):
  module_(module), TBMChannel_(TBMChannel)
{
	std::cout<<" pixelchannel1 "<<module_.modulename()<<" "<<TBMChannel_<<std::endl;
}

PixelChannel::PixelChannel(PixelModuleName module, PixelTBMChannel TBMChannel):
  module_(module), TBMChannel_(TBMChannel)
{
	std::cout<<" pixelchannel2 "<<module_.modulename()<<" "<<TBMChannel_<<std::endl;
}

PixelChannel::PixelChannel(std::string name)
{
	module_ = PixelModuleName(name);
	char TBMChannelString[2] = {0,0};
	TBMChannelString[0] = name[name.size()-1]; // take the last character of name
	TBMChannel_ = PixelTBMChannel(TBMChannelString);

}

std::ostream& pos::operator<<(std::ostream& s, const PixelChannel& channel)
{
  s << channel.channelname();
  return s;
}

std::string PixelChannel::channelname() const
{
  //std::cout<<" pixelchannelname "<<TBMChannelString<<std::endl;
  return modulename() + "_ch" + TBMChannelStringFull();
}
