#include <iostream>
#include <fstream>
#include <cstdio>
#include <algorithm>
#include <cstring>
#include <cstdlib>
//#include <string>
using namespace std;

#define capacity 100
#define base 0x40000000  // 2^30
#define __base 0x20000000  // 2^29
//#define base_dec_bit 9
#define base_bit 30
#define base_mod 0x3FFFFFFF
#define exponentiation_size 1500
#define some_prime 25

using namespace std;

ofstream fout;
ofstream public_out;
ofstream private_out;



const int Prime[some_prime] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97};

inline int sgn(int key)
{
	if (key > 0) return 1;
	if (key < 0) return -1;
	if (key == 0) return 0;
}

inline int abs(int key)
{
	return (key < 0) ? -key : key;
}

struct UnsignedBigInt
{
	int len;
	int __data[capacity];
	UnsignedBigInt() : len(0) {}
	UnsignedBigInt(const UnsignedBigInt &source) : len(source.len)
	{
		memcpy(__data, source.__data, len * sizeof *__data);
	}
	UnsignedBigInt(int key) : len(0)
	{
		for (; key > 0; key >>= base_bit)
			__data[len ++] = key & base_mod;
	}
/*
	UnsignedBigInt(long long key) : len(0)
	{
		for (; key > 0; key /= base)
			__data[len ++] = key % base;
	}
*/
	UnsignedBigInt & operator = (const UnsignedBigInt &key)
	{
		len = key.len;
		memcpy(__data, key.__data, len * sizeof *__data);
		return *this;
	}
	int & operator [] (int Index) { return __data[Index]; }
	int operator [] (int Index) const { return __data[Index]; }
};

inline int compare(const UnsignedBigInt &A, const UnsignedBigInt &B)
{
	if (A.len != B.len) return A.len > B.len ? 1 : -1;
	int i;
	for (i = A.len - 1; i >= 0 && A[i] == B[i]; -- i);
	if (i < 0) return 0;
	return A[i] > B[i] ? 1 : -1;
}


inline void shift_right_base(UnsignedBigInt &A)
{
	if (A.len == 0) return;
	for (int i = 0; i < A.len - 1; ++ i)
		A[i] = A[i + 1];
	-- A.len;
}

inline void shift_left_base(UnsignedBigInt &A)
{
	if (A.len == 0) return;
	for (int i = A.len; i > 0; -- i)
		A[i] = A[i - 1];
	A[0] = 0;
	++ A.len;
}

inline void shift_right(UnsignedBigInt &A)
{
	if (A.len == 0) return;
	for (int i = 0; i < A.len - 1; ++ i)
	{
		A[i] = (A[i] >> 1);
		if ((A[i + 1] & 1) == 1)
			A[i] = A[i] | __base;
	}
	A[A.len - 1] = A[A.len - 1] >> 1;
	if (A[A.len - 1] == 0) -- A.len;
}

inline void shift_left(UnsignedBigInt &A)
{
	if (A.len == 0) return;
	int k = A.len;
	if ((A[A.len - 1] & __base) != 0)
		A[A.len ++] = 1;
	for (int i = k - 1; i > 0; -- i)
	{
		A[i] = (A[i] << 1) & base_mod;
		if ((A[i - 1] & __base) != 0)
			++ A[i];
	}
	A[0] = (A[0] << 1) & base_mod;
}




UnsignedBigInt operator + (const UnsignedBigInt &A, const UnsignedBigInt &B)
{
	UnsignedBigInt R;
	int i;
	int Carry = 0;
	for (i = 0; i < A.len || i < B.len || Carry > 0; ++ i)
	{
		if (i < A.len) Carry += A[i];
		if (i < B.len) Carry += B[i];
		R[i] = Carry & base_mod;
		Carry >>= base_bit;
	}
	R.len = i;
	return R;
	
}

UnsignedBigInt operator - (const UnsignedBigInt &A, const UnsignedBigInt &B)
{
	UnsignedBigInt R;
	int Carry = 0;
	R.len = A.len;
	for (int i = 0; i < R.len; ++ i)
	{
		R[i] = A[i] - Carry;
		if (i < B.len) R[i] -= B[i];
		if (R[i] < 0) Carry = 1, R[i] += base;
		else Carry = 0;
	}
	while (R.len > 0 && R[R.len - 1] == 0) -- R.len;
	return R;
}

UnsignedBigInt operator * (const UnsignedBigInt &A, const int B)
{
	int i;
	if (B == 0) return 0;
	UnsignedBigInt R;
	long long Carry = 0;
	for (i = 0; i < A.len || Carry > 0; ++ i)
	{
		if (i < A.len) Carry += (long long)(A[i]) * B;
		R[i] = Carry & base_mod;
		Carry >>= base_bit;
	}
	R.len = i;
	return R;
}

UnsignedBigInt operator * (const UnsignedBigInt &A, const UnsignedBigInt &B)
{
	if (B.len == 0) return 0;
	UnsignedBigInt R;
	for (int i = 0; i < A.len; ++ i)
	{
		long long Carry = 0;
		for (int j = 0; j < B.len || Carry > 0; ++ j)
		{
			if (j < B.len) Carry += (long long)(A[i]) * B[j];
			if (i + j < R.len) Carry += R[i + j];
			if (i + j >= R.len) R[R.len ++] = Carry & base_mod;
			else R[i + j] = Carry & base_mod;
			Carry >>= base_bit;
		}
	}
	return R;
}



UnsignedBigInt operator / (const UnsignedBigInt &A, const int B)
{

	UnsignedBigInt R;
	long long C = 0;
	for (int i = A.len - 1; i >= 0; -- i)
	{
		C = (C << base_bit) + A[i];
		R[i] = C / B;
		C %= B;
	}
	R.len = A.len;
	while (R.len > 0 && R[R.len - 1] == 0) -- R.len;
	return R;

}

UnsignedBigInt operator % (const UnsignedBigInt &A, const int B)
{
	long long C = 0;
	for (int i = A.len - 1; i >= 0; -- i)
	{
		C = (C << base_bit) + A[i];
		C %= B;
	}
	return C;
}


inline void divide(const UnsignedBigInt &A, const UnsignedBigInt &B, UnsignedBigInt &Q, UnsignedBigInt &R)
{
	static bool flag[exponentiation_size];  //odd ---> true, even ---> false
	static int top;
	static UnsignedBigInt temp_divide;
	
	temp_divide = A;
	for (top = 0; compare(temp_divide, B) >= 0; ++ top)
	{
//		cout << temp_divide.len << endl;
		flag[top] = ((temp_divide[0] & 1) == 1);
		shift_right(temp_divide);
	}
	
	Q = 0; R = temp_divide;
/*
	R = 0;
	int left, right, mid;
	for (int i = temp_divide.len - 1; i >= 0; -- i)
	{
		R = R * base + temp_divide[i];
		left = 0;
		right = base;
		while (left + 1 < right)
		{
			mid = (left + right) >> 1;
			if (compare(B * mid , R) <= 0)
				left = mid;
			else
				right = mid;
		}
		Q[i] = left;
		R = R - B * left;
	}
	Q.len = temp_divide.len;
	while (Q.len > 0 && Q[Q.len - 1] == 0) -- Q.len;
*/	
	for (int i = top - 1; i >= 0; -- i)
	{
		shift_left(Q); shift_left(R);
		if (flag[i])
			R = R + 1;
		if (compare(R, B) >= 0)
		{
			R = R - B;
			Q = Q + 1;
		}
	}
	//cout << "*************" << endl;
}


UnsignedBigInt operator / (const UnsignedBigInt &A, const UnsignedBigInt &B)
{
	if (compare(A, B) < 0) return 0;
/*
	UnsignedBigInt R, Q;
	divide(A, B, Q, R);
	return Q;
*/


	UnsignedBigInt R, Carry = 0;
	int left, right, mid;
	for (int i = A.len - 1; i >= 0; -- i)
	{
		//Carry = Carry * base + A[i];
		shift_left_base(Carry);
		Carry = Carry + A[i];

		left = 0;
		right = base;
		while (left + 1 < right)
		{
			mid = (left + right) >> 1;
			if (compare(B * mid , Carry) <= 0)
				left = mid;
			else
				right = mid;
		}
		R[i] = left;
		Carry = Carry - B * left;
	}
	R.len = A.len;
	while (R.len > 0 && R[R.len - 1] == 0) -- R.len;
	return R;

}

UnsignedBigInt operator % (const UnsignedBigInt &A, const UnsignedBigInt &B)
{
	if (compare(A, B) < 0) return A;
/*
	UnsignedBigInt R, Q;
	divide(A, B, Q, R);
	return R;

*/
	UnsignedBigInt Carry = 0;
	int left, right, mid;
	for (int i = A.len - 1; i >= 0; -- i)
	{
		//Carry = Carry * base + A[i];
		shift_left_base(Carry);
		Carry = Carry + A[i];
		
		left = 0;
		right = base;
		while (left + 1 < right)
		{
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


struct signed_BigInt
{
	UnsignedBigInt data;
	int sign;
	signed_BigInt() : data(0), sign(0) {}
	signed_BigInt(const signed_BigInt &source) : sign(source.sign)
	{
		data = source.data;
	}
	signed_BigInt(int key) : sign(sgn(key)), data(abs(key)) {}
	signed_BigInt & operator = (const signed_BigInt &key)
	{
		sign = key.sign;
		data = key.data;
		return *this;
	}
	int & operator [] (int Index) { return data[Index]; }
	int operator [] (int Index) const { return data[Index]; }
};




const UnsignedBigInt unsigned_Zero = 0;
const UnsignedBigInt unsigned_One = 1;

const signed_BigInt Zero = 0;
const signed_BigInt One = 1;




struct output_BigInt
{
	int len;
	int __data[capacity];
	output_BigInt() : len(0) {}
	output_BigInt(const output_BigInt &source) : len(source.len)
	{
		memcpy(__data, source.__data, len * sizeof *__data);
	}
	output_BigInt(int key) : len(0)
	{
		for (; key > 0; key /= 1000000000)
			__data[len ++] = key % 1000000000;
	}
	int & operator [] (int Index) { return __data[Index]; }
	int operator [] (int Index) const { return __data[Index]; }
};


output_BigInt operator + (const output_BigInt &A, const output_BigInt &B)
{
	output_BigInt R;
	int i;
	int Carry = 0;
	for (i = 0; i < A.len || i < B.len || Carry > 0; ++ i)
	{
		if (i < A.len) Carry += A[i];
		if (i < B.len) Carry += B[i];
		R[i] = Carry % 1000000000;
		Carry /= 1000000000;
	}
	R.len = i;
	return R;	
}

output_BigInt operator * (const output_BigInt &A, const int B)
{
	int i;
	if (B == 0) return 0;
	output_BigInt R;
	long long Carry = 0;
	for (i = 0; i < A.len || Carry > 0; ++ i)
	{
		if (i < A.len) Carry += (long long)(A[i]) * B;
		R[i] = Carry % 1000000000;
		Carry /= 1000000000;
	}
	R.len = i;
	return R;
}



ostream & operator << (ostream &Out, const UnsignedBigInt &A)
{
	if (A.len == 0)
	{
		Out << 0;
		return Out;
	}

	output_BigInt V = 0;
	for (int i = A.len - 1; i >= 0; -- i)
		V = V * base + A[i];


	Out << V[V.len - 1];
	for (int i = V.len - 2; i >= 0; -- i)
		for (int j = 100000000; j > 0; j /= 10)
			Out << V[i] / j % 10;
	Out << endl;
	
	if (Out == fout || Out == cout)
	{
		int digit = 0;
		for (int k = 1; k <= V[V.len - 1]; k *= 10)
			++ digit;
	
		digit += (V.len - 1) * 9;
	
		Out << "Total Digit = " << digit << endl;
	}
	return Out;
}


ostream & operator << (ostream &Out, const signed_BigInt &V)
{
	string c = (V.sign == -1) ? "-" : "";
	Out << c << V.data;
	return Out;
}

bool operator == (const signed_BigInt &A, const signed_BigInt &B)
{
	if (A.sign != B.sign) return false;
	return (compare(A.data, B.data) == 0);
}

signed_BigInt operator + (const signed_BigInt &A, const signed_BigInt &B)
{
	if (A.sign == 0) return B;
	if (B.sign == 0) return A;
	signed_BigInt R;
	if (A.sign * B.sign == 1)
	{
		R.data = A.data + B.data;
		R.sign = A.sign;
	}
	else if (A.sign == 1)
	{
		int k = compare(A.data, B.data);
		if (k == 1)
		{
			R.sign = 1;
			R.data = A.data - B.data;
		}
		else if (k == 0)
		{
			R.sign = 0;
			R.data = 0;
		}
		else if (k == -1)
		{
			R.sign = -1;
			R.data = B.data - A.data;
		}
	}
	else
	{
		int k = compare(B.data, A.data);
		if (k == 1)
		{
			R.sign = 1;
			R.data = B.data - A.data;
		}
		else if (k == 0)
		{
			R.sign = 0;
			R.data = 0;
		}
		else if (k == -1)
		{
			R.sign = -1;
			R.data = A.data - B.data;
		}
	}
	return R;
}

signed_BigInt operator - (const signed_BigInt &A)
{
	signed_BigInt R;
	R.sign = -A.sign;
	R.data = A.data;
	return R;
}

signed_BigInt operator - (const signed_BigInt &A, const signed_BigInt &B)
{
	if (B.sign == 0) return A;
	if (A.sign == 0) return -B;
	signed_BigInt R;
	if (A.sign * B.sign == -1)
	{
		R.sign = A.sign;
		R.data = A.data + B.data;
	}
	else
	{
		int k = compare(A.data, B.data);
		if (k == 0)
		{
			R.sign = 0;
			R.data = 0;
			return R;
		}
		if (A.sign == 1 && B.sign == 1)
		{
			if (k == -1)
			{
				R.data = B.data - A.data;
				R.sign = -1;
			}
			else
			{
				R.data = A.data - B.data;
				R.sign = 1;
			}
		}
		else
		{
			if (k == -1)
			{
				R.data = B.data - A.data;
				R.sign = 1;
			}
			else
			{
				R.data = A.data - B.data;
				R.sign = -1;
			}
		}
	}
	return R;
}

signed_BigInt operator * (const signed_BigInt &A, const signed_BigInt &B)
{
	signed_BigInt R;
	if (A.sign * B.sign == 0)
	{
		R.data = 0;
		R.sign = 0;
		return R;
	}
	R.data = A.data * B.data;
	R.sign = A.sign * B.sign;
	return R;
}

signed_BigInt operator / (const signed_BigInt &A, const signed_BigInt &B)
{
	signed_BigInt R;
	if (A.sign == 0)
	{
		R.data = 0;
		R.sign = 0;
		return R;
	}
	R.data = A.data / B.data;
	if (R.data.len == 0)
		R.sign = 0;
	else
		R.sign = A.sign * B.sign;
	return R;
}

signed_BigInt operator % (const signed_BigInt &A, const signed_BigInt &B)
{
	signed_BigInt R;
	if (A.sign == 0)
	{
		R.data = 0;
		R.sign = 0;
		return R;
	}
	R.data = A.data % B.data;
	if (R.data.len == 0)
	{
		R.sign = 0;
		return R;
	}
	R.sign = 1;
	if (A.sign * B.sign == -1)
		R.data = B.data - A.data % B.data;
	return R;
}




signed_BigInt Euclid_GCD(const signed_BigInt &A, const signed_BigInt &B)
{
	return (B.sign == 0) ? A : Euclid_GCD(B, A % B);
}


signed_BigInt Extended_Euclid_GCD(const signed_BigInt &A, const signed_BigInt &B, signed_BigInt &X, signed_BigInt &Y)
{
	if (B.sign == 0)
	{
		X = 1, Y = 0;
		return A;
	}
	else
	{
		signed_BigInt XX, YY;
		signed_BigInt D = Extended_Euclid_GCD(B, A % B, XX, YY);
		X = YY, Y = XX - (A / B) * YY;
		return D;
	}
	
}



struct BigInt_Exponentiation
{
	UnsignedBigInt N;
	int len, valid;
	int data[exponentiation_size];
	BigInt_Exponentiation() : len(0) { memset(data, 0, sizeof(data)); }
	BigInt_Exponentiation(const BigInt_Exponentiation &source) : len(source.len)
	{
		memcpy(data, source.data, len * sizeof *data);
	}
/*	
	BigInt_Exponentiation(UnsignedBigInt A) : len(0)
	{
		for (; A.len != 0; A = A / 2)
			data[len ++] = (A[0] & 1);
	}
*/
	int & operator [] (int Index) { return data[Index]; }
	int operator [] (int Index) const { return data[Index]; }
};




UnsignedBigInt Bin_To_Int(const BigInt_Exponentiation &W)
{
	UnsignedBigInt R;
	R.len = 0;
	int i, k;
	for (i = 0; i + base_bit < W.len; i += base_bit)
	{
		k = 0;
		for (int j = base_bit - 1; j >= 0; -- j)
			k = (k << 1) + W[i + j];
		R[R.len ++] = k;
	}
	if (i < W.len)
	{
		k = 0;
		for (int j = W.len - 1; j >= i; -- j)
			k = (k << 1) + W[j];
		R[R.len ++] = k;
	}
	return R;
}


/*
UnsignedBigInt __Bin_To_Int(const BigInt_Exponentiation &W)
{
	UnsignedBigInt R = 0;
	for (int i = W.len - 1; i >= 0; -- i)
	{
		R = R * 2;
		if (W[i] == 1)
			R = R + 1;
	}
	return R;
}

*/

UnsignedBigInt Minus_One;

/*
UnsignedBigInt __Modular_Exponentiation(UnsignedBigInt A, UnsignedBigInt B, const UnsignedBigInt &N)
{
	UnsignedBigInt R = 1;
	while (B.len != 0)
	{
		if ((B[0] & 1) == 1)
			R = (R * A) % N;
		A = (A * A) % N;
		B = B / 2;
	}
	return R;
}
*/


UnsignedBigInt Modular_Exponentiation(UnsignedBigInt A, const BigInt_Exponentiation &W)
{
	UnsignedBigInt R = 1;
	for (int i = W.valid; i < W.len; ++ i)
	{
		if (W[i] == 1)
			R = (R * A) % W.N;
		A = (A * A) % W.N;
	}
	return R;
}



/*

int limit;


inline int optimal(const int bits)
{
	if (bits <= 9) return 1;
	if (bits <= 25) return 2;
	if (bits <= 70) return 3;
	if (bits <= 197) return 4;
	if (bits <= 539) return 5;
	if (bits <= 1434) return 6;
	if (bits <= 3715) return 7;
	cout << "********" << endl;
}


UnsignedBigInt Exponentiation_Temp[128];
int Optimal_K;






UnsignedBigInt Modular_Exponentiation_Valid(const UnsignedBigInt &A, const BigInt_Exponentiation &W)
{
//	BigInt_Exponentiation W = B;
//	cout <<"-----------------" << W.len << endl;




	Exponentiation_Temp[1] = A;
	if (limit >= 3)
	{
		UnsignedBigInt temp = (A * A) % W.N;
		for (int i = 3; i <= limit; i += 2)
			Exponentiation_Temp[i] = (Exponentiation_Temp[i - 2] * temp) % W.N;
	}
	UnsignedBigInt R = 1;
	int i = W.len - 1;
	while (i >= W.valid)
		if (W[i] == 0)
		{
			R = (R * R) % W.N;
			-- i;
		}
		else
		{
			int s = i - Optimal_K + 1;
			if (s < W.valid) s = W.valid;
			while (W[s] == 0) ++ s;
			for (int j = 0; j < i - s + 1; ++ j)
				R = (R * R) % W.N;
			int u = 0;
			for (int j = i; j >= s; -- j)
				u = ((u << 1) + W[j]);
			R = (R * Exponentiation_Temp[u]) % W.N;
			i = s - 1;
		}
	return R;
}

*/




BigInt_Exponentiation Get_Random_Binary(const int full, const int empty)   // LSB --> MSB   scheme: (10000)empty(1xxxxx1)full
{
	BigInt_Exponentiation X;
	X.len = full + empty;
	X.valid = empty;  // Start from empty
	X[0] = X[empty] = 1;
	for (int i = empty + 1; i < X.len - 1; ++ i)
		X[i] = (rand() & 1);
	X[X.len - 1] = 1;
	X.N = Bin_To_Int(X);
//
//	if (compare(X.N, __Bin_To_Int(X)) != 0)
//	{
//		cout << X.N << endl;
//		cout << __Bin_To_Int(X) << endl;
//		for (int i = 0; i < X.len; ++ i)
//			cout << X[i];
//		cout << endl;
//		cout << "*******************" << endl;
//	}
	return X;
}




/*

UnsignedBigInt Get_Random(const int bits, const int style)		//style 0 -> all; 1 -> odd; -1 -> even; 2 -> odd and last bit != 5
{
	UnsignedBigInt R = 0;
	if (bits > 1) R = (rand() % 9 + 1);
	int i;
	for (i = 2; i + base_dec_bit - 1 < bits; i += base_dec_bit)
	{
		int k = 0;
		for (int j = 1; j <= base_dec_bit; ++ j)
			k = k * 10 + rand() % 10;
		R = R * base + k;
	}
	for (int j = i; j < bits; ++ j)
		R = R * 10 + (rand() % 10);
	if (style == 2)
	{
		int k = rand() % 5;
		while (k == 2) k = rand() % 5;
		R = R * 10 + (k * 2 + 1);
	}
	else if (style == 0)
		R = R * 10 + (rand() % 10);
	else if (style == 1)
		R = R * 10 + ((rand() % 5) * 2 + 1);
	else if (style == -1)
		R = R * 10 + ((rand() % 5) * 2);
	return R;
}

*/

bool Miller_Rabin_Witness(const UnsignedBigInt &A, const BigInt_Exponentiation &W)
{
	int t = W.valid;
	UnsignedBigInt u = Minus_One;
//	for (t = 0; (u[0] & 1) == 0; ++ t)
//		u = u / 2;
	UnsignedBigInt X[2];
//	X[0] = Modular_Exponentiation_Valid(A, W);
	X[0] = Modular_Exponentiation(A, W);
//	cout << X[0] << endl;
	for (int i = 1; i <= t; ++ i)
	{
		X[i & 1] = (X[(i - 1) & 1] * X[(i - 1) & 1]) % W.N;
		if (compare(X[i & 1], unsigned_One) == 0 && compare(X[(i - 1) & 1], unsigned_One) != 0 && compare(X[(i - 1) & 1], Minus_One) != 0)
//		{
//			cout << "//////" << endl;
			return true;
//		}
	}
	
	if (compare(X[t & 1], unsigned_One) != 0)
//	{
//		cout << "***********" << endl;
		return true;
//	}
	
	return false;
}

bool Miller_Rabin_Primality_Test(const BigInt_Exponentiation &W, const int Trial)
{
/*
	if ((W.N % 5).len == 0) return false;
	if ((W.N % 3).len == 0) return false;
	if ((W.N % 7).len == 0) return false;
*/
	for (int i = 0; i < some_prime; ++ i)
		if ((W.N % Prime[i]).len == 0) return false;
	Minus_One = W.N - 1;
	UnsignedBigInt A;

/*	

	Optimal_K = optimal(W.len - W.valid);
//	int k = 3;

	int cur = W.len - 1;
	limit = 1;
	while (cur >= W.valid)
		if (W[cur] == 0)
			-- cur;
		else
		{
			int s = cur - Optimal_K + 1;
			if (s < W.valid) s = W.valid;
			while (W[s] == 0) ++ s;
			int u = 0;
			for (int j = cur; j >= s; -- j)
				u = ((u << 1) + W[j]);
			if (u >= limit) limit = u;
			cur = s - 1;
		}
*/

	if (Miller_Rabin_Witness(2, W)) return false;
	for (int i = 1; i < Trial; ++ i)
	{
		
		A = Get_Random_Binary(W.len - 2, 1).N;
//		if (compare(A, Minus_One) > 0)
//			cout << "/////////////////////" << endl;
//		A = A % Minus_One + 1;
//		cout << "///////////////" << endl;
		if (Miller_Rabin_Witness(A, W)) return false;
//		k += 2;
	}
	return true;
}


signed_BigInt Get_Prime(const int Prime_Bits, const int Trial)
{
	int delta = rand() % (Prime_Bits / 20) + 1;
	int bits = Prime_Bits + delta;

	int noisy = rand() % 100 + 100;

	int counter = 0;
	BigInt_Exponentiation W;
	while (1)
	{
		++ counter;
//		R.data = Get_Random(bits, 2);
		W = Get_Random_Binary(bits - noisy, noisy);    // empty >= 1
		//cout << W.N << endl;
//		cout << counter << endl;
//		cout << "* " << total_multiply <<  " / " <<  total_divide << " % " << total_mod << endl;
//		if (counter == 15) return 0;
		if (Miller_Rabin_Primality_Test(W, Trial))
		{
			cout << "Total Trial = " << counter << endl;
			signed_BigInt R;
			R.sign = 1;
			R.data = W.N;
			return R;
		}
	}
}


int Big_Prime_Bits, Trial;

void Generating()
{
	
	signed_BigInt P = Get_Prime(Big_Prime_Bits, Trial), Q = Get_Prime(Big_Prime_Bits, Trial);
	signed_BigInt Phi = (P - 1) * (Q - 1);
	signed_BigInt N = P * Q;

	signed_BigInt D, E, temp;
	for (E = 7; !(Extended_Euclid_GCD(E, Phi, D, temp) == One); E = E + 2);
	if (D.sign == -1)
		D = D + Phi;
		

	fout << "P = " << endl;
	fout << P << endl;
	fout << endl;
	fout << "Q = " << endl;
	fout << Q << endl;
	fout << endl;
	
	fout << "N = " << endl;
	fout << N << endl;
	fout << endl;
	fout << "E = " << endl;
	fout << E << endl;
	fout << endl;
	fout << "D = " << endl;
	fout << D << endl;

	public_out << N;
	public_out << E;
	
	
	private_out << N;
	private_out << D;
	
	fout.close();
	public_out.close();
	private_out.close();
}


int main(int argc, char *argv[])
{
	system("clear");
	cout <<	"**************************************************************" << endl;
	cout << "*                   RSA Public Key Generator                 *" << endl;
	cout << "*                                  Copyright © 2011, Jian Li *" << endl;
	cout << "*  for help, please using --help                             *" << endl;
	cout << "**************************************************************" << endl;

	srand((unsigned)time(NULL));


	if (argc > 2)
	{
		if (strcmp(argv[1], "--high-performance") == 0)
		{
			Big_Prime_Bits = 800;
			Trial = 6;			
	
		}
		else if (strcmp(argv[1], "--high-reliability") == 0)
		{
			Big_Prime_Bits = 1200;
			Trial = 10;
		}
		else
		{
			cout << "Format error! Please type --help for help." << endl;
			return 0;
		}
		char c[100];
		strcpy(c, argv[2]);
		strcpy(c + strlen(c), "_RSA_Key");
		fout.open(c);
		strcpy(c, argv[2]);
		strcpy(c + strlen(c), "_RSA_Public_Key");
		public_out.open(c);
		strcpy(c, argv[2]);
		strcpy(c + strlen(c), "_RSA_Private_Key");
		private_out.open(c);
		cout << "Please Wait..." << endl;
		Generating();
		cout << "Successful!" << endl;
		return 0;
	}

	if (argc > 1)
	{
		if (strcmp(argv[1], "--high-performance") == 0)
		{
			Big_Prime_Bits = 800;
			Trial = 6;
			fout.open("RSA_Key");
			public_out.open("RSA_Public_Key");
			private_out.open("RSA_Private_Key");
			cout << "Please Wait..." << endl;
			Generating();
			cout << "Successful!" << endl;
			return 0;			
	
		}
		if (strcmp(argv[1], "--high-reliability") == 0)
		{
			Big_Prime_Bits = 1200;
			Trial = 10;
			fout.open("RSA_Key");
			public_out.open("RSA_Public_Key");
			private_out.open("RSA_Private_Key");
			cout << "Please Wait..." << endl;
			Generating();
			cout << "Successful!" << endl;
			return 0;
		}
		if (strcmp(argv[1], "--help") == 0)
		{
			cout << "(0) By default, it will generate a 250-digits Key" << endl;
			cout << "(1) For high-performance, please using --high-performance, it will generate a 500-digits Key" << endl;
			cout << "(2) For high-reliability, please using --high-reliability, it will generate a 750-digits Key, but it takes a long time..." << endl;
			cout << "(3) For digitial-signature, please using --digital-signature" << endl;
			cout << "(4) For getting a large prime number, please using --prime" << endl;
			cout << "(5) Default output file RSA_Key, RSA_Public_Key, RSA_Private_Key, you can change it by add your name follow your command." << endl;
			cout << endl;
			cout << "Example:" << endl;
			cout << "\"unix> ./rsa Alice\", program will output file Alice_RSA_Key, Alice_RSA_Public_Key, Alice_Private_Key" << endl;
			cout << endl;
			cout << "(?) For help, please using --help" << endl;
			cout << "Good Luck!" << endl;
			return 0;
		}
		if (strcmp(argv[1], "--digital-signature") == 0)
		{
			cout << "Generating Alice RSA Key..." << endl;
			fout.open("Alice_RSA_Key");
			public_out.open("Alice_RSA_Public_Key");
			private_out.open("Alice_RSA_Private_Key");
			Trial = 6;
			Big_Prime_Bits = 300;
			Generating();
			
			cout << "Generating Bob RSA Key..." << endl;
			fout.open("Bob_RSA_Key");
			public_out.open("Bob_RSA_Public_Key");
			private_out.open("Bob_RSA_Private_Key");
			Trial = 6;
			Big_Prime_Bits = 400;
			Generating();
			
			cout << "Successful!" << endl;
			return 0;
		}
		if (strcmp(argv[1], "--prime") == 0)
		{
			cout << "How large prime would you want(in Binary Bits)?" << endl;
			cout << "Please input a number <= 1500" << endl;
			Trial = 10; 
			cin >> Big_Prime_Bits;
			cout << "Please Wait..." << endl;
			cout << Get_Prime(Big_Prime_Bits, Trial) << endl;
			return 0;
		}
		Big_Prime_Bits = 400;		// just for fun
		Trial = 5;
		char c[100];
		strcpy(c, argv[1]);
		strcpy(c + strlen(c), "_RSA_Key");
		fout.open(c);
		strcpy(c, argv[1]);
		strcpy(c + strlen(c), "_RSA_Public_Key");
		public_out.open(c);
		strcpy(c, argv[1]);
		strcpy(c + strlen(c), "_RSA_Private_Key");
		private_out.open(c);
		cout << "Please Wait..." << endl;
		Generating();
		cout << "Successful!" << endl;
		return 0;
	}
	Big_Prime_Bits = 400;		// just for fun
	Trial = 5;
	fout.open("RSA_Key");
	public_out.open("RSA_Public_Key");
	private_out.open("RSA_Private_Key");
	cout << "Please Wait..." << endl;
	Generating();
	cout << "Successful!" << endl;
	return 0;



/*

	UnsignedBigInt A = 30, B = 7;
	cout << A % B << endl;
	cout << A / B << endl;
*/
/*
	UnsignedBigInt A = 1;
	for (int i = 0; i < 40; ++ i)
	{
		A = shift_left(A);
		cout << A << endl;
	}
	for (int i = 0; i < 40; ++ i)
	{
		A = shift_right(A);
		cout << A << endl;
	}
*/
/*
	BigInt_Exponentiation W(1729);

	W.valid = 5;
	W.N = Bin_To_Dec(W);

	cout << W.N << endl;
	cout << Miller_Rabin_Primality_Test(W, Trial) << endl;
*/
//	BigInt_Exponentiation P = Get_Random_Binary(12, 7);
/*
	signed_BigInt P = Get_Prime(Big_Prime_Bits, Trial);
	cout << P << endl;
*/
/*	
	UnsignedBigInt P = Get_Random(200, 0);

	int i = 12343;
		for (int j = 1000000000; ; ++ j)
		{
			cout << i << " " << j << endl;
			if (compare(Modular_Exponentiation(i, j, P),  __Modular_Exponentiation(i, j, P)) != 0)
			{

				return 0;
			}
		}
*/	
//	UnsignedBigInt A = Modular_Exponentiation(3, 6, 1000000);
//	cout << A << endl;
	
	
//     Important  signed_BigInt      A * const   not defined
	
//	BigInt_Exponentiation P =  Get_Random_Binary(10, 5);




	

	//cout << Modular_Exponentiation(3123, E.data * D.data, N.data) << endl;
}

