#ifndef RSA_BIGINT_
#define RSA_BIGINT_

#include <vector>

namespace bigint {

using BaseType = int;

constexpr int kBaseBit = 30;
constexpr int kBase = (1 << kBaseBit);
constexpr int kBaseMask = kBase - 1;
constexpr int kBaseMSB = (1 << (kBaseBit - 1));

class UnsignedBigInt {
 public:
  UnsignedBigInt();
  UnsignedBigInt(const UnsignedBigInt &);
  UnsignedBigInt(int);

  UnsignedBigInt &operator=(const UnsignedBigInt &);
  BaseType &operator[](unsigned int);
  BaseType operator[](unsigned int) const;

  int Len() const;
  void Clean();

  void ShiftBaseLeft();
  void ShiftBaseRight();
  void ShiftRight();
  void ShiftLeft();

 private:
  std::vector<BaseType> digit_;
};

int Compare(const UnsignedBigInt &, const UnsignedBigInt &);

UnsignedBigInt operator+(const UnsignedBigInt &, const UnsignedBigInt &);
UnsignedBigInt operator-(const UnsignedBigInt &, const UnsignedBigInt &);

}  // namespace bigint
#endif  // RSA_BIGINT_