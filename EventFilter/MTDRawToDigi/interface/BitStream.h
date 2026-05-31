#ifndef EventFilter_MTDRawToDigi_interface_BitStream_h
#define EventFilter_MTDRawToDigi_interface_BitStream_h

#include <cstdint>
#include <cstddef>
#include <type_traits>
#include <algorithm>
#include <array>

template <std::size_t DIM>
class BitStream {
  static_assert(DIM > 0, "DIM must be > 0");

public:
  using word_t = uint64_t;
  static constexpr std::size_t WORD_BITS = 64;
  static constexpr std::size_t BITS = DIM * WORD_BITS;

  explicit BitStream();

  // Return stored bits number
  size_t size() const { return BITS; }

  // Clear all bits
  void clear();

  // Write a value (lower 'bitCount' bits of value) starting at bit index 'pos' (0 = LSB of data[0]).
  // Throws std::out_of_range if pos+bitCount > BITS or bitCount > 64.
  void set_bits(std::size_t pos, std::size_t bitCount, uint64_t value);

  // Read bitCount bits starting at pos as a uint64_t (returned in low bits). Throws on OOB or bitCount>64.
  uint64_t get_bits(std::size_t pos, std::size_t bitCount) const;

  // Direct access to underlying array (const and non-const)
  const word_t* raw_data() const { return data; }
  word_t* raw_data() { return data; }

  // For convenience: set a single bit
  void set_bit(std::size_t pos, bool value = true);

  // For convenience: get a single bit
  bool get_bit(std::size_t pos) const;

private:
  std::array<word_t, DIM> data;
};

#include "EventFilter/MTDRawToDigi/interface/BitStream.icc"

#endif
