#include "s2cell_id.h"

#include <algorithm>
#include <mutex>
#include <vector>
#include <cfloat>

using S2::internal::kSwapMask;
using S2::internal::kInvertMask;

using S2::internal::kIJtoPos;
using S2::internal::kPosToIJ;
using S2::internal::kPosToOrientation;
using S2::internal::kFaceUVWAxes;
using S2::internal::kFaceUVWFaces;

using std::max;
using std::min;

const int S2CellId::kFaceBits;
const int S2CellId::kNumFaces;
const int S2CellId::kMaxLevel;
const int S2CellId::kPosBits;
const int S2CellId::kMaxSize;

static const int kLookupBits = 4;
static uint16 lookup_pos[1 << (2 * kLookupBits + 2)];
static uint16 lookup_ij[1 << (2 * kLookupBits + 2)];

static void InitLookupCell(int level, int i, int j, int orig_orientation,
                           int pos, int orientation) {
  if (level == kLookupBits) {
    int ij = (i << kLookupBits) + j;
    lookup_pos[(ij << 2) + orig_orientation] = (pos << 2) + orientation;
    lookup_ij[(pos << 2) + orig_orientation] = (ij << 2) + orientation;
  } else {
    level++;
    i <<= 1;
    j <<= 1;
    pos <<= 2;
    const int* r = kPosToIJ[orientation];
    InitLookupCell(level, i + (r[0] >> 1), j + (r[0] & 1), orig_orientation,
                   pos, orientation ^ kPosToOrientation[0]);
    InitLookupCell(level, i + (r[1] >> 1), j + (r[1] & 1), orig_orientation,
                   pos + 1, orientation ^ kPosToOrientation[1]);
    InitLookupCell(level, i + (r[2] >> 1), j + (r[2] & 1), orig_orientation,
                   pos + 2, orientation ^ kPosToOrientation[2]);
    InitLookupCell(level, i + (r[3] >> 1), j + (r[3] & 1), orig_orientation,
                   pos + 3, orientation ^ kPosToOrientation[3]);
  }
}

static std::once_flag flag;
inline static void MaybeInit() {
  std::call_once(flag, []{
    InitLookupCell(0, 0, 0, 0, 0, 0);
    InitLookupCell(0, 0, 0, kSwapMask, 0, kSwapMask);
    InitLookupCell(0, 0, 0, kInvertMask, 0, kInvertMask);
    InitLookupCell(0, 0, 0, kSwapMask|kInvertMask, 0, kSwapMask|kInvertMask);
  });
}


S2CellId S2CellId::FromFaceIJ(int face, int i, int j) 
{
  MaybeInit();
  uint64 n = static_cast<uint64>(face) << (kPosBits - 1);
  uint64 bits = (face & kSwapMask);

#define GET_BITS(k) do { \
    const int mask = (1 << kLookupBits) - 1; \
    bits += ((i >> (k * kLookupBits)) & mask) << (kLookupBits + 2); \
    bits += ((j >> (k * kLookupBits)) & mask) << 2; \
    bits = lookup_pos[bits]; \
    n |= (bits >> 2) << (k * 2 * kLookupBits); \
    bits &= (kSwapMask | kInvertMask); \
  } while (0)

  GET_BITS(7);
  GET_BITS(6);
  GET_BITS(5);
  GET_BITS(4);
  GET_BITS(3);
  GET_BITS(2);
  GET_BITS(1);
  GET_BITS(0);
#undef GET_BITS

  return S2CellId(n * 2 + 1);
}

void S2CellId::AppendVertexNeighbors(int level,
                                     std::vector<S2CellId>* output) const {
  // "level" must be strictly less than this cell's level so that we can
  // determine which vertex this cell is closest to.
  S2_DCHECK_LT(level, this->level());
  int i, j;
  int face = ToFaceIJOrientation(&i, &j, nullptr);

  // Determine the i- and j-offsets to the closest neighboring cell in each
  // direction.  This involves looking at the next bit of "i" and "j" to
  // determine which quadrant of this->parent(level) this cell lies in.
  int halfsize = GetSizeIJ(level + 1);
  int size = halfsize << 1;
  bool isame, jsame;
  int ioffset, joffset;
  if (i & halfsize) {
    ioffset = size;
    isame = (i + size) < kMaxSize;
  } else {
    ioffset = -size;
    isame = (i - size) >= 0;
  }
  if (j & halfsize) {
    joffset = size;
    jsame = (j + size) < kMaxSize;
  } else {
    joffset = -size;
    jsame = (j - size) >= 0;
  }

  output->push_back(parent(level));
  output->push_back(FromFaceIJSame(face, i + ioffset, j, isame).parent(level));
  output->push_back(FromFaceIJSame(face, i, j + joffset, jsame).parent(level));
  // If i- and j- edge neighbors are *both* on a different face, then this
  // vertex only has three neighbors (it is one of the 8 cube vertices).
  if (isame || jsame) {
    output->push_back(FromFaceIJSame(face, i + ioffset, j + joffset,
                                     isame && jsame).parent(level));
  }
}

void S2CellId::AppendAllNeighbors(int nbr_level,
                                  std::vector<S2CellId>* output) const {
  S2_DCHECK_GE(nbr_level, level());
  int i, j;
  int face = ToFaceIJOrientation(&i, &j, nullptr);

  // Find the coordinates of the lower left-hand leaf cell.  We need to
  // normalize (i,j) to a known position within the cell because nbr_level
  // may be larger than this cell's level.
  int size = GetSizeIJ();
  i &= -size;
  j &= -size;

  int nbr_size = GetSizeIJ(nbr_level);
  S2_DCHECK_LE(nbr_size, size);

  // We compute the top-bottom, left-right, and diagonal neighbors in one
  // pass.  The loop test is at the end of the loop to avoid 32-bit overflow.
  for (int k = -nbr_size; ; k += nbr_size) {
    bool same_face;
    if (k < 0) {
      same_face = (j + k >= 0);
    } else if (k >= size) {
      same_face = (j + k < kMaxSize);
    } else {
      same_face = true;
      // Top and bottom neighbors.
      output->push_back(FromFaceIJSame(face, i + k, j - nbr_size,
                                       j - size >= 0).parent(nbr_level));
      output->push_back(FromFaceIJSame(face, i + k, j + size,
                                       j + size < kMaxSize).parent(nbr_level));
    }
    // Left, right, and diagonal neighbors.
    output->push_back(FromFaceIJSame(face, i - nbr_size, j + k,
                                     same_face && i - size >= 0)
                      .parent(nbr_level));
    output->push_back(FromFaceIJSame(face, i + size, j + k,
                                     same_face && i + size < kMaxSize)
                      .parent(nbr_level));
    if (k >= size) break;
  }
}

S2CellId S2CellId::FromFaceIJWrap(int face, int i, int j) {
  // Convert i and j to the coordinates of a leaf cell just beyond the
  // boundary of this face.  This prevents 32-bit overflow in the case
  // of finding the neighbors of a face cell.
  i = max(-1, min(kMaxSize, i));
  j = max(-1, min(kMaxSize, j));

  // We want to wrap these coordinates onto the appropriate adjacent face.
  // The easiest way to do this is to convert the (i,j) coordinates to (x,y,z)
  // (which yields a point outside the normal face boundary), and then call
  // S2::XYZtoFaceUV() to project back onto the correct face.
  //
  // The code below converts (i,j) to (si,ti), and then (si,ti) to (u,v) using
  // the linear projection (u=2*s-1 and v=2*t-1).  (The code further below
  // converts back using the inverse projection, s=0.5*(u+1) and t=0.5*(v+1).
  // Any projection would work here, so we use the simplest.)  We also clamp
  // the (u,v) coordinates so that the point is barely outside the
  // [-1,1]x[-1,1] face rectangle, since otherwise the reprojection step
  // (which divides by the new z coordinate) might change the other
  // coordinates enough so that we end up in the wrong leaf cell.
  static const double kScale = 1.0 / kMaxSize;
  static const double kLimit = 1.0 + DBL_EPSILON;
  // The arithmetic below is designed to avoid 32-bit integer overflows.
  S2_DCHECK_EQ(0, kMaxSize % 2);
  double u = max(-kLimit, min(kLimit, kScale * (2 * (i - kMaxSize / 2) + 1)));
  double v = max(-kLimit, min(kLimit, kScale * (2 * (j - kMaxSize / 2) + 1)));

  // Find the leaf cell coordinates on the adjacent face, and convert
  // them to a cell id at the appropriate level.
  face = S2::XYZtoFaceUV(S2::FaceUVtoXYZ(face, u, v), &u, &v);
  return FromFaceIJ(face, S2::STtoIJ(0.5*(u+1)), S2::STtoIJ(0.5*(v+1)));
}

inline S2CellId S2CellId::FromFaceIJSame(int face, int i, int j,
                                         bool same_face) {
  if (same_face)
    return S2CellId::FromFaceIJ(face, i, j);
  else
    return S2CellId::FromFaceIJWrap(face, i, j);
}

void S2CellId::GetEdgeNeighbors(S2CellId neighbors[4]) const {
  int i, j;
  int level = this->level();
  int size = GetSizeIJ(level);
  int face = ToFaceIJOrientation(&i, &j, nullptr);
  // Edges 0, 1, 2, 3 are in the down, right, up, left directions.
  neighbors[0] = FromFaceIJSame(face, i, j - size, j - size >= 0)
                 .parent(level);
  neighbors[1] = FromFaceIJSame(face, i + size, j, i + size < kMaxSize)
                 .parent(level);
  neighbors[2] = FromFaceIJSame(face, i, j + size, j + size < kMaxSize)
                 .parent(level);
  neighbors[3] = FromFaceIJSame(face, i - size, j, i - size >= 0)
                 .parent(level);
}

int S2CellId::ToFaceIJOrientation(int* pi, int* pj, int* orientation) const {
  // Initialization if not done yet
  MaybeInit();

  int i = 0, j = 0;
  int face = this->face();
  int bits = (face & kSwapMask);

  // Each iteration maps 8 bits of the Hilbert curve position into
  // 4 bits of "i" and "j".  The lookup table transforms a key of the
  // form "ppppppppoo" to a value of the form "iiiijjjjoo", where the
  // letters [ijpo] represents bits of "i", "j", the Hilbert curve
  // position, and the Hilbert curve orientation respectively.
  //
  // On the first iteration we need to be careful to clear out the bits
  // representing the cube face.
#define GET_BITS(k) do { \
    const int nbits = (k == 7) ? (kMaxLevel - 7 * kLookupBits) : kLookupBits; \
    bits += (static_cast<int>(id_ >> (k * 2 * kLookupBits + 1)) \
             & ((1 << (2 * nbits)) - 1)) << 2; \
    bits = lookup_ij[bits]; \
    i += (bits >> (kLookupBits + 2)) << (k * kLookupBits); \
    j += ((bits >> 2) & ((1 << kLookupBits) - 1)) << (k * kLookupBits); \
    bits &= (kSwapMask | kInvertMask); \
  } while (0)

  GET_BITS(7);
  GET_BITS(6);
  GET_BITS(5);
  GET_BITS(4);
  GET_BITS(3);
  GET_BITS(2);
  GET_BITS(1);
  GET_BITS(0);
#undef GET_BITS

  *pi = i;
  *pj = j;

  if (orientation != nullptr) {
    // The position of a non-leaf cell at level "n" consists of a prefix of
    // 2*n bits that identifies the cell, followed by a suffix of
    // 2*(kMaxLevel-n)+1 bits of the form 10*.  If n==kMaxLevel, the suffix is
    // just "1" and has no effect.  Otherwise, it consists of "10", followed
    // by (kMaxLevel-n-1) repetitions of "00", followed by "0".  The "10" has
    // no effect, while each occurrence of "00" has the effect of reversing
    // the kSwapMask bit.
    S2_DCHECK_EQ(0, kPosToOrientation[2]);
    S2_DCHECK_EQ(kSwapMask, kPosToOrientation[0]);
    if (lsb() & 0x1111111111111110ULL) {
      bits ^= kSwapMask;
    }
    *orientation = bits;
  }
  return face;
}

// lat lon -> cell id
S2CellId::S2CellId(const S2Point& p)
{
  double u, v;
  int face = S2::XYZtoFaceUV(p, &u, &v);
  int i = S2::STtoIJ(S2::UVtoST(u));
  int j = S2::STtoIJ(S2::UVtoST(v));
  id_ = FromFaceIJ(face, i, j).id();
}

S2CellId::S2CellId(const S2LatLng& ll)
  : S2CellId(ll.ToPoint()) {
}

S2Point S2CellId::ToPointRaw() const {
  int si, ti;
  int face = GetCenterSiTi(&si, &ti);
  return S2::FaceSiTitoXYZ(face, si, ti);
}

S2LatLng S2CellId::ToLatLng() const {
  return S2LatLng(ToPointRaw());
}



