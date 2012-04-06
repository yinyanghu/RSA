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
#define ascii_size 128


using namespace std;


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

struct unsigned_BigInt
{
	int len;
	int __data[capacity];
	unsigned_BigInt() : len(0) {}
	unsigned_BigInt(const unsigned_BigInt &source) : len(source.len)
	{
		memcpy(__data, source.__data, len * sizeof *__data);
	}
	unsigned_BigInt(int key) : len(0)
	{
		for (; key > 0; key >>= base_bit)
			__data[len ++] = key & base_mod;
	}

	unsigned_BigInt & operator = (const unsigned_BigInt &key)
	{
		len = key.len;
		memcpy(__data, key.__data, len * sizeof *__data);
		return *this;
	}
	int & operator [] (int Index) { return __data[Index]; }
	int operator [] (int Index) const { return __data[Index]; }
};

inline int compare(const unsigned_BigInt &A, const unsigned_BigInt &B)
{
	if (A.len != B.len) return A.len > B.len ? 1 : -1;
	int i;
	for (i = A.len - 1; i >= 0 && A[i] == B[i]; -- i);
	if (i < 0) return 0;
	return A[i] > B[i] ? 1 : -1;
}


inline void shift_right_base(unsigned_BigInt &A)
{
	if (A.len == 0) return;
	for (int i = 0; i < A.len - 1; ++ i)
		A[i] = A[i + 1];
	-- A.len;
}

inline void shift_left_base(unsigned_BigInt &A)
{
	if (A.len == 0) return;
	for (int i = A.len; i > 0; -- i)
		A[i] = A[i - 1];
	A[0] = 0;
	++ A.len;
}

inline void shift_right(unsigned_BigInt &A)
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

inline void shift_left(unsigned_BigInt &A)
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




unsigned_BigInt operator + (const unsigned_BigInt &A, const unsigned_BigInt &B)
{
	unsigned_BigInt R;
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

unsigned_BigInt operator - (const unsigned_BigInt &A, const unsigned_BigInt &B)
{
	unsigned_BigInt R;
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

unsigned_BigInt operator * (const unsigned_BigInt &A, const int B)
{
	int i;
	if (B == 0) return 0;
	unsigned_BigInt R;
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

unsigned_BigInt operator * (const unsigned_BigInt &A, const unsigned_BigInt &B)
{
	if (B.len == 0) return 0;
	unsigned_BigInt R;
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



unsigned_BigInt operator / (const unsigned_BigInt &A, const int B)
{

	unsigned_BigInt R;
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

unsigned_BigInt operator % (const unsigned_BigInt &A, const int B)
{
	long long C = 0;
	for (int i = A.len - 1; i >= 0; -- i)
	{
		C = (C << base_bit) + A[i];
		C %= B;
	}
	return C;
}


unsigned_BigInt operator / (const unsigned_BigInt &A, const unsigned_BigInt &B)
{
	if (compare(A, B) < 0) return 0;

	unsigned_BigInt R, Carry = 0;
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

unsigned_BigInt operator % (const unsigned_BigInt &A, const unsigned_BigInt &B)
{
	if (compare(A, B) < 0) return A;

	unsigned_BigInt Carry = 0;
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




const unsigned_BigInt unsigned_Zero = 0;
const unsigned_BigInt unsigned_One = 1;




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

istream & operator >> (istream &In, unsigned_BigInt &A)
{
	char ch;
	for (A = 0; In >> ch;)
	{
		A = A * 10 + (ch - '0');
		if (In.peek() <= ' ') break;
	}
	return In;
}

ostream & operator << (ostream &Out, const unsigned_BigInt &A)
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
	return Out;
}



unsigned_BigInt Modular_Exponentiation(unsigned_BigInt A, unsigned_BigInt B, const unsigned_BigInt &N)
{
	unsigned_BigInt R = 1;
	while (B.len != 0)
	{
		if ((B[0] & 1) == 1)
			R = (R * A) % N;
		A = (A * A) % N;
		shift_right(B);
	}
	return R;
}


string decoder(unsigned_BigInt A)
{
	unsigned_BigInt temp;
	string S = "";
	for (; A.len != 0; A = A / ascii_size)
	{
		temp = A % ascii_size;
		char ch = char(temp[0]);
		S = ch + S;
	}
	return S;
}

int main(int argc, char *argv[])
{

	system("clear");
	cout <<	"**************************************************************" << endl;
	cout << "*                Digital Signature Decoder                   *" << endl;
	cout << "*                                              Base on RSA   *" << endl;
	cout << "*                                  Copyright © 2011, Jian Li *" << endl;
	cout << "*  for help, please using --help                             *" << endl;
	cout << "**************************************************************" << endl;
		
	unsigned_BigInt R;
	
	char s[100];
	char t[100];	
	char u[100];
	
	if (argc > 1)
	{
		if (strcmp(argv[1], "--name") == 0)
		{
			cout << "Please input digital signature file" << endl;
			cin >> s;


			cout << "Please input your Private Key file" << endl;
			cin >> t;


			cout << "Please input sender's Public Key file" << endl;
			cin >> u;
		}
		else if (strcmp(argv[1], "--help") == 0)
		{
			cout << "(1) For input custom file name, please using --name" << endl;
			cout << "(2) Default input file are \"Signature\", \"Bob_RSA_Private_Key\", \"Alice_RSA_Public_Key\"" << endl;
			cout << "(?) For help, please using --help" << endl;
			cout << "Have Fun!" << endl;
			return 0;
		}
		else
		{
			cout << "Error!" << endl;
			cout << "Please type --help to get more help!" << endl;
			return 0;
		}
	}
	else
	{
		strcpy(s, "Signature");
		strcpy(t, "Bob_RSA_Private_Key");
		strcpy(u, "Alice_RSA_Public_Key");
	}
	
	ifstream DS(s);
	DS >> R;
	DS.close();
//	cout << R << endl;
	unsigned_BigInt N_A, N_B, D, E;
	

	
	ifstream Private_File(t);
	Private_File >> N_A >> D;
	Private_File.close();
	
	R = Modular_Exponentiation(R, D, N_A);
//	cout << R << endl;
	

	ifstream Public_File(u);
	Public_File >> N_B >> E;
	Public_File.close();
	
	R = Modular_Exponentiation(R, E, N_B);
//	cout << R << endl;
	
	cout << "Successful!" << endl;
	cout << "Your digital signature for the sender is" << endl;
	
	cout <<	"**************************************************************" << endl;
	cout << decoder(R) << endl;
	cout <<	"**************************************************************" << endl;	
	cout << endl;
	cout << "Please Check!" << endl;
	return 0;
}

