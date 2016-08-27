#ifndef PixelTBMChannel_h
#define PixelTBMChannel_h
/*! \file CalibFormats/SiPixelObjects/interface/PixelTBMChannel.h
*   \brief Simple class to hold either "A" or "B" for the TBM channel
*
*    A longer explanation will be placed here later
*/
 
#include <string>
#include <iostream>

namespace pos{
/*! \class PixelTBMChannel PixelTBMChannel.h "interface/PixelTBMChannel.h"
*   \brief Simple class to hold either "A" or "B" for the TBM channel
*
*   A longer explanation will be placed here later
*/
  class PixelTBMChannel
  {
    public:

    PixelTBMChannel(){;}
    PixelTBMChannel(std::string TBMChannel);

    std::string string() const;
    std::string stringFull() const;

    friend std::ostream& operator<<(std::ostream& s, const PixelTBMChannel& TBMChannel);

    // What does this do? Revisit, dk 6/16
    const bool operator<(const PixelTBMChannel& aTBMChannel) const{
      //return ( isChannelB_ == false && aTBMChannel.isChannelB_ == true ); //old one for bools 
      return ( channel_ < aTBMChannel.channel_ ); // does this work
    }

    const bool operator==(const PixelTBMChannel& aTBMChannel) const{
      return (channel_== aTBMChannel.channel_);
      //return isChannelB_==aTBMChannel.isChannelB_;
    }

    private:
    bool isChannelB_;
    int channel_;
  };
  std::ostream& operator<<(std::ostream& s, const PixelTBMChannel& TBMChannel);
}
#endif
