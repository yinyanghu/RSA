#ifndef RSA_BIGINT_
#define RSA_BIGINT_

#include <vector>

namespace bigint {

class UnsignedBigInt {
  UnsignedBigInt();
  UnsignedBigInt(const UnsignedBigInt &);
  UnsignedBigInt(int);

  UnsignedBigInt &operator=(const UnsignedBigInt &);
  int &operator[](int);
  int operator[](int) const;

 private:
  std::vector<int> digit_;
};

}  // namespace bigint
#endif  // RSA_BIGINT_