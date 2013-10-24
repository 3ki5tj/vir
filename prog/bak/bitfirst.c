/* using _bit_scan_forward */
#ifdef __INTEL_COMPILER /* use intrinsics */


/* index of nonzero bit */
INLINE int bitfirst(code_t x)
{
#if CODEBITS == 32
  return (x) ? _bit_scan_forward(x) : 0;
#elif CODEBITS == 64
  if (x == 0) {
    return 0;
  } else {
    uint32_t low = (uint32_t) (x & 0xffffffff);
    if (low) return _bit_scan_forward(low);
    else return _bit_scan_forward((x >> 32) & 0xffffffff) + 32;
  }
#endif
}


#define BITFIRSTLOW(id, x, b) { (b) = (x) & (-(x)); id = bitfirst(x); }

INLINE int bitfirstlow(code_t x, code_t *b)
{
  *b = x & (-x); /* such that only the lowest 1-bit survives */
  return bitfirst(x);
}



#endif
