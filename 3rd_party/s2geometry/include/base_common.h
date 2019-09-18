//
// Copyright 2019 Wayz.ai. All Rights Reserved.
//
#ifndef WAYZ_BASE_COMMON_H_
#define WAYZ_BASE_COMMON_H_

#include <glog/logging.h>

namespace S2 {

using int8 = signed char;
using int16 = short;
using int32 = int;
using int64 = long long;

using uint8 = unsigned char;
using uint16 = unsigned short;
using uint32 = unsigned int;
using uint64 = unsigned long long;
using uword_t = unsigned long;

#define S2_DCHECK(condition) PCHECK(condition)
#define S2_DCHECK_GT(val1, val2) CHECK_GT(val1, val2)
#define S2_DCHECK_LE(val1, val2) CHECK_LE(val1, val2)
#define S2_DLOG_IF(severity, condition) LOG_IF(severity, condition)
#define S2_DCHECK_GE(val1, val2) CHECK_GE(val1, val2)
#define S2_DCHECK_LT(val1, val2) CHECK_LT(val1, val2)
#define S2_DCHECK_EQ(val1, val2) CHECK_EQ(val1, val2)

class Bits 
{
public:
    static int FindLSBSetNonZero64(uint64 n);

};
inline int Bits::FindLSBSetNonZero64(uint64 n) {
  return __builtin_ctzll(n);
}

}

#endif