/* Include files */

#include "receiveDBPSK____PLUTO_SDR_cgxe.h"
#include "m_Bul6X3hpblE5kWOkaudxB.h"

unsigned int cgxe_receiveDBPSK____PLUTO_SDR_method_dispatcher(SimStruct* S,
  int_T method, void* data)
{
  if (ssGetChecksum0(S) == 1896788184 &&
      ssGetChecksum1(S) == 1721710460 &&
      ssGetChecksum2(S) == 4176125477 &&
      ssGetChecksum3(S) == 1625364854) {
    method_dispatcher_Bul6X3hpblE5kWOkaudxB(S, method, data);
    return 1;
  }

  return 0;
}
