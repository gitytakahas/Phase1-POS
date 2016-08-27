// #ifndef dump_error_files_masks_h
// #define dump_error_files_masks_h

namespace masks
{
typedef const unsigned long type;

// data words
type plsmsk = 0xff;   // pulse height
type pxlmsk = 0xff00; // pixel index
type dclmsk = 0x1f0000;
type rocmsk = 0x3e00000;
type chnlmsk = 0xfc000000;

// Error words
type errorMask = 0x3e00000;     // (0x1f << 21)
type dummyMask = 0x03600000;    // (0x1b << 21)
type gapMask = 0x03400000;      // (0x1a << 21)
type timeOut = 0x3a00000;       // (0x1d << 21)
type eventNumError = 0x3e00000; // (0x1f << 21)
type trailError = 0x3c00000;    // (0x1e << 21)
type fifoError = 0x3800000;     // (0x1c << 21) nearly full
type channelMask = 0xfc000000;  // channel num mask
type eventNumMask = 0x1fe000;   // event number mask

// for ENE
type tbmEventMask = 0xff;      // tbm event num mask

// for trailer errors
type overflowMask = 0x100;     // data overflow
type FsmErrMask = 0x600;     // pointer to FSM errors
type RocErrMask = 0x800;     // pointer to #Roc errors
type tbmStatusMask = 0xff;   // TBM trailer info

// for timeout words
type BlkNumMask = 0x700;     // pointer to error fifo #
type MSB_counter = 0xfc000000;
type LSB_counter = 0x1800;
type TO_data_chs = 0x1f;

// for FIFO / nearly full words
type Fif2NFMask = 0x40; // mask for fifo2 NF
type TrigNFMask = 0x80; // mask for trigger fifo NF
type ChnFifMask = 0x1f; // channel mask
type ChnFifMask0 = 0x1; // channel mask
type ChnFifMask1 = 0x2; // channel mask
type ChnFifMask2 = 0x4; // channel mask
type ChnFifMask3 = 0x8; // channel mask
type ChnFifMask4 = 0x10;// channel mask

    
// offsets for channels in a block
const int offsets[8] = {0, 4, 9, 13, 18, 22, 27, 31};
}

// #endif