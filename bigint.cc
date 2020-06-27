#include "bigint.h"

#include <vector>

namespace bigint {

UnsignedBigInt::UnsignedBigInt() {}

UnsignedBigInt::UnsignedBigInt(const UnsignedBigInt &o) {
  this->digit_ = o.digit_;
}

UnsignedBigInt::UnsignedBigInt(int num) {
  for (; num; num >>= kBaseBit) digit_.push_back(num & kBaseMask);
}

UnsignedBigInt &UnsignedBigInt::operator=(const UnsignedBigInt &o) {
  this->digit_ = o.digit_;
  return *this;
}

BaseType &UnsignedBigInt::operator[](unsigned int index) {
  while (index >= digit_.size()) {
    digit_.push_back(0);
  }
  return digit_[index];
}

BaseType UnsignedBigInt::operator[](unsigned int index) const {
  if (index >= digit_.size()) {
    return 0;
  }
  return digit_[index];
}

int UnsignedBigInt::Len() const { return digit_.size(); }

void UnsignedBigInt::Clean() {
  for (int i = digit_.size() - 1; i >= 0 && digit_[i] == 0; --i) {
    digit_.pop_back();
  }
}

void UnsignedBigInt::ShiftBaseLeft() {
  int l = digit_.size();
  if (l == 0) {
    return;
  }
  digit_.push_back(digit_[l - 1]);
  for (int i = l - 1; i; --i) {
    digit_[i] = digit_[i - 1];
  }
  digit_[0] = 0;
}

void UnsignedBigInt::ShiftBaseRight() {
  int l = digit_.size();
  if (l == 0) {
    return;
  }
  for (int i = 0; i < l - 1; ++i) {
    digit_[i] = digit_[i + 1];
  }
  digit_.pop_back();
}

void UnsignedBigInt::ShiftRight() {
  int l = digit_.size();
  if (l == 0) {
    return;
  }
  for (int i = 0; i < l - 1; ++i) {
    digit_[i] = (digit_[i] >> 1);
    if (digit_[i + 1] & 1) {
      digit_[i] = digit_[i] | kBaseMSB;
    }
  }
  digit_[l - 1] >>= 1;
  if (digit_[l - 1] == 0) {
    digit_.pop_back();
  }
}

void UnsignedBigInt::ShiftLeft() {
  int l = digit_.size();
  if (l == 0) {
    return;
  }
  if (digit_[l - 1] & kBaseMSB) {
    digit_.push_back(1);
  }
  for (int i = l - 1; i >= 0; --i) {
    digit_[i] = (digit_[i] << 1) & kBaseMask;
    if (i && (digit_[i - 1] & kBaseMSB)) {
      ++digit_[i];
    }
  }
}

int Compare(const UnsignedBigInt &x, const UnsignedBigInt &y) {
  if (x.Len() != y.Len()) {
    return x.Len() > y.Len() ? 1 : -1;
  }
  int i;
  for (i = x.Len() - 1; i >= 0 && x[i] == y[i]; --i)
    ;
  if (i < 0) {
    return 0;
  }
  return x[i] > y[i] ? 1 : -1;
}

UnsignedBigInt operator+(const UnsignedBigInt &a, const UnsignedBigInt &b) {
  UnsignedBigInt c;
  BaseType carry = 0;
  for (int i = 0; i < std::max(a.Len(), b.Len()) || carry; ++i) {
    carry += a[i] + b[i];
    c[i] = carry & kBaseMask;
    carry >>= kBaseBit;
  }
  return c;
}

UnsignedBigInt operator-(const UnsignedBigInt &a, const UnsignedBigInt &b) {
  UnsignedBigInt c;
  int carry = 0;
  for (int i = 0; i < a.Len(); ++i) {
    c[i] = a[i] - b[i] - carry;
    if (c[i] < 0) {
      carry = 1, c[i] += kBase;
    } else {
      carry = 0;
    }
  }
  c.Clean();
  return c;
}

}  // namespace bigint

/*

inline int sgn(int key) {
  if (key > 0) return 1;
  if (key < 0) return -1;
  if (key == 0) return 0;
}

int abs(int key) { return (key < 0) ? -key : key; }

UnsignedBigInt operator*(const UnsignedBigInt &A, const int B) {
  int i;
  if (B == 0) return 0;
  UnsignedBigInt R;
  long long Carry = 0;
  for (i = 0; i < A.len || Carry > 0; ++i) {
    if (i < A.len) Carry += (long long)(A[i]) * B;
    R[i] = Carry & base_mod;
    Carry >>= base_bit;
  }
  R.len = i;
  return R;
}

UnsignedBigInt operator*(const UnsignedBigInt &A, const UnsignedBigInt &B) {
  if (B.len == 0) return 0;
  UnsignedBigInt R;
  for (int i = 0; i < A.len; ++i) {
    long long Carry = 0;
    for (int j = 0; j < B.len || Carry > 0; ++j) {
      if (j < B.len) Carry += (long long)(A[i]) * B[j];
      if (i + j < R.len) Carry += R[i + j];
      if (i + j >= R.len)
        R[R.len++] = Carry & base_mod;
      else
        R[i + j] = Carry & base_mod;
      Carry >>= base_bit;
    }
  }
  return R;
}

UnsignedBigInt operator/(const UnsignedBigInt &A, const int B) {
  UnsignedBigInt R;
  long long C = 0;
  for (int i = A.len - 1; i >= 0; --i) {
    C = (C << base_bit) + A[i];
    R[i] = C / B;
    C %= B;
  }
  R.len = A.len;
  while (R.len > 0 && R[R.len - 1] == 0) --R.len;
  return R;
}

UnsignedBigInt operator%(const UnsignedBigInt &A, const int B) {
  long long C = 0;
  for (int i = A.len - 1; i >= 0; --i) {
    C = (C << base_bit) + A[i];
    C %= B;
  }
  return C;
}

UnsignedBigInt operator/(const UnsignedBigInt &A, const UnsignedBigInt &B) {
  if (compare(A, B) < 0) return 0;

  UnsignedBigInt R, Carry = 0;
  int left, right, mid;
  for (int i = A.len - 1; i >= 0; --i) {
    // Carry = Carry * base + A[i];
    shift_left_base(Carry);
    Carry = Carry + A[i];

    left = 0;
    right = base;
    while (left + 1 < right) {
      mid = (left + right) >> 1;
      if (compare(B * mid, Carry) <= 0)
        left = mid;
      else
        right = mid;
    }
    R[i] = left;
    Carry = Carry - B * left;
  }
  R.len = A.len;
  while (R.len > 0 && R[R.len - 1] == 0) --R.len;
  return R;
}

UnsignedBigInt operator%(const UnsignedBigInt &A, const UnsignedBigInt &B) {
  if (compare(A, B) < 0) return A;

  UnsignedBigInt Carry = 0;
  int left, right, mid;
  for (int i = A.len - 1; i >= 0; --i) {
    // Carry = Carry * base + A[i];
    shift_left_base(Carry);
    Carry = Carry + A[i];

    left = 0;
    right = base;
    while (left + 1 < right) {
      mid = (left + right) >> 1;
      if (compare(B * mid, Carry) <= 0)
        left = mid;
      else
        right = mid;
    }
    Carry = Carry - B * left;
  }
  return Carry;
}

const UnsignedBigInt unsigned_Zero = 0;
const UnsignedBigInt unsigned_One = 1;

struct output_BigInt {
  int len;
  int __data[capacity];
  output_BigInt() : len(0) {}
  output_BigInt(const output_BigInt &source) : len(source.len) {
    memcpy(__data, source.__data, len * sizeof *__data);
  }
  output_BigInt(int key) : len(0) {
    for (; key > 0; key /= 1000000000) __data[len++] = key % 1000000000;
  }
  int &operator[](int Index) { return __data[Index]; }
  int operator[](int Index) const { return __data[Index]; }
};

output_BigInt operator+(const output_BigInt &A, const output_BigInt &B) {
  output_BigInt R;
  int i;
  int Carry = 0;
  for (i = 0; i < A.len || i < B.len || Carry > 0; ++i) {
    if (i < A.len) Carry += A[i];
    if (i < B.len) Carry += B[i];
    R[i] = Carry % 1000000000;
    Carry /= 1000000000;
  }
  R.len = i;
  return R;
}

output_BigInt operator*(const output_BigInt &A, const int B) {
  int i;
  if (B == 0) return 0;
  output_BigInt R;
  long long Carry = 0;
  for (i = 0; i < A.len || Carry > 0; ++i) {
    if (i < A.len) Carry += (long long)(A[i]) * B;
    R[i] = Carry % 1000000000;
    Carry /= 1000000000;
  }
  R.len = i;
  return R;
}

istream &operator>>(istream &In, UnsignedBigInt &A) {
  char ch;
  for (A = 0; In >> ch;) {
    A = A * 10 + (ch - '0');
    if (In.peek() <= ' ') break;
  }
  return In;
}

ostream &operator<<(ostream &Out, const UnsignedBigInt &A) {
  if (A.len == 0) {
    Out << 0;
    return Out;
  }

  output_BigInt V = 0;
  for (int i = A.len - 1; i >= 0; --i) V = V * base + A[i];

  Out << V[V.len - 1];
  for (int i = V.len - 2; i >= 0; --i)
    for (int j = 100000000; j > 0; j /= 10) Out << V[i] / j % 10;
  return Out;
}

UnsignedBigInt Modular_Exponentiation(UnsignedBigInt A, UnsignedBigInt B,
                                      const UnsignedBigInt &N) {
  UnsignedBigInt R = 1;
  while (B.len != 0) {
    if ((B[0] & 1) == 1) R = (R * A) % N;
    A = (A * A) % N;
    shift_right(B);
  }
  return R;
}

*/