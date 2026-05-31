#include "EventFilter/MTDRawToDigi/interface/BitStream.h"

#include <iostream>
#include <iomanip>
#include <bitset>
#include <sstream>

int main() {
  const size_t DIM = 2;
  BitStream<DIM> pl;

  auto printBits = [&](const uint64_t& id) {
    std::stringstream ss;
    ss << std::bitset<4>((id >> 60) & 0xF).to_string() << " " << std::bitset<4>((id >> 56) & 0xF).to_string() << " "
       << std::bitset<4>((id >> 52) & 0xF).to_string() << " " << std::bitset<4>((id >> 48) & 0xF).to_string() << " "
       << std::bitset<4>((id >> 44) & 0xF).to_string() << " " << std::bitset<4>((id >> 40) & 0xF).to_string() << " "
       << std::bitset<4>((id >> 36) & 0xF).to_string() << " " << std::bitset<4>((id >> 32) & 0xF).to_string() << " "
       << std::bitset<4>((id >> 28) & 0xF).to_string() << " " << std::bitset<4>((id >> 24) & 0xF).to_string() << " "
       << std::bitset<4>((id >> 20) & 0xF).to_string() << " " << std::bitset<4>((id >> 16) & 0xF).to_string() << " "
       << std::bitset<4>((id >> 12) & 0xF).to_string() << " " << std::bitset<4>((id >> 8) & 0xF).to_string() << " "
       << std::bitset<4>((id >> 4) & 0xF).to_string() << " " << std::bitset<4>(id & 0xF).to_string();
    return ss.str();
  };

  auto printStream = [&](const BitStream<DIM>& bs) {
    std::stringstream ss;
    for (int index = static_cast<int>(bs.size()); index >= 0; --index) {
      ss << bs.get_bit(static_cast<size_t>(index));
    }
    return ss.str();
  };

  uint64_t pippo = 0xFAFAFAAF;

  std::cout << pippo << " --> " << printBits(pippo) << std::endl;

  std::string offset("                       ");
  std::stringstream scale;
  for (int index = static_cast<int>(pl.size()); index >= 0; --index) {
    scale << index % 10;
    std::cout << std::fixed << std::setw(3) << index << " " << scale.str() << std::endl;
  }

  std::cout << "\n" << std::endl;
  std::cout << offset << scale.str() << std::endl;
  std::cout << offset << printStream(pl) << std::endl;

  size_t pos(0), bitCount(32);

  pl.set_bits(pos, bitCount, pippo);
  std::cout << std::fixed << "pos = " << std::setw(2) << pos << " bitCount = " << bitCount << " " << printStream(pl)
            << std::endl;
  pl.clear();
  std::cout << offset << printStream(pl) << std::endl;

  pos = 0;
  bitCount = 16;
  pl.set_bits(pos, bitCount, pippo);
  std::cout << std::fixed << "pos = " << std::setw(2) << pos << " bitCount = " << bitCount << " " << printStream(pl)
            << std::endl;
  pl.clear();

  pos = 32;
  bitCount = 32;
  pl.set_bits(pos, bitCount, pippo);
  std::cout << std::fixed << "pos = " << std::setw(2) << pos << " bitCount = " << bitCount << " " << printStream(pl)
            << std::endl;
  pl.clear();

  pos = 50;
  bitCount = 32;
  pl.set_bits(pos, bitCount, pippo);
  std::cout << std::fixed << "pos = " << std::setw(2) << pos << " bitCount = " << bitCount << " " << printStream(pl)
            << std::endl;
  pl.clear();

  return 0;
}
