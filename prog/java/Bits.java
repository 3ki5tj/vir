/** Bitwise operations.
 *  All functions are `static' for no instance is needed */
class Bits {
  /** Return a long integer with only ith lowest bit being 1 */
  public static long makeBit(int i) {
    return (long) 1 << i;
  }

  /** Return a long integer with the lowest i bits being 1 */
  public static long makeBitsMask(int i) {
    return ((long) 1 << i) - 1;
  }

  private static final int [] bitrem = {
    -1,  0,  1, 39,  2, 15, 40, 23,  3, 12, 16, 59, 41, 19, 24, 54,
     4, -1, 13, 10, 17, 62, 60, 28, 42, 30, 20, 51, 25, 44, 55, 47,
     5, 32, -1, 38, 14, 22, 11, 58, 18, 53, 63,  9, 61, 27, 29, 50,
    43, 46, 31, 37, 21, 57, 52,  8, 26, 49, 45, 36, 56,  7, 48, 35,
     6, 34, 33};

  /** From bit to id */
  public static int bit2id(long bit) {
    return bitrem[(int) (bit % 67)];
  }

  /** Seek the least significant bit */
  public static int bitFirst(long x) {
    long bit = x & (-x);
    return bit2id(bit);
  }

  private static final int [] byteBitCount = {
    0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4,
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
    4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8 };

  /** Count the number of one-bits in x */
  public static int bitCount(long x) {
    return byteBitCount[(int) (x         & 0xffl)] +
           byteBitCount[(int) ((x >>  8) & 0xffl)] +
           byteBitCount[(int) ((x >> 16) & 0xffl)] +
           byteBitCount[(int) ((x >> 24) & 0xffl)] +
           byteBitCount[(int) ((x >> 32) & 0xffl)] +
           byteBitCount[(int) ((x >> 40) & 0xffl)] +
           byteBitCount[(int) ((x >> 48) & 0xffl)] +
           byteBitCount[(int) ((x >> 56) & 0xffl)];
  }
}




