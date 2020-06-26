#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>

using namespace std;

UnsignedBigInt encoder(const string &X) {
  UnsignedBigInt R = 0;
  for (int i = 0; i < X.size(); ++i) R = (R * ascii_size) + int(X[i]);
  return R;
}

int main(int argc, char *argv[]) {
  system("clear");
  cout << "**************************************************************"
       << endl;
  cout << "*                Digital Signature Encoder                   *"
       << endl;
  cout << "*                                              Base on RSA   *"
       << endl;
  cout << "*                                  Copyright © 2011, Jian Li *"
       << endl;
  cout << "*  for help, please using --help                             *"
       << endl;
  cout << "**************************************************************"
       << endl;

  //	cout << R << endl;
  UnsignedBigInt N_A, N_B, D, E;

  char s[100];
  char t[100];
  if (argc > 1) {
    if (strcmp(argv[1], "--name") == 0) {
      cout << "Please input your Private Key file" << endl;
      cin >> s;

      cout << "Please input receiver's Public Key file" << endl;
      cin >> t;

    } else if (strcmp(argv[1], "--help") == 0) {
      cout << "(1) For input custom file name, please using --name" << endl;
      cout << "(2) Default input file are \"Alice_RSA_Private_Key\", "
              "\"Bob_RSA_Public_Key\""
           << endl;
      cout << "(?) For help, please using --help" << endl;
      cout << "Wow!" << endl;
      return 0;
    } else {
      cout << "Error!" << endl;
      cout << "Please type --help to get more help!" << endl;
      return 0;
    }
  } else {
    strcpy(s, "Alice_RSA_Private_Key");
    strcpy(t, "Bob_RSA_Public_Key");
  }

  string DS;
  cout << "Please sign your name" << endl;

  // cin >> DS;
  getline(cin, DS);

  UnsignedBigInt R = encoder(DS);

  //	char s[] = "Alice_RSA_Private_Key";
  ifstream Private_File(s);

  Private_File >> N_A >> D;
  Private_File.close();

  R = Modular_Exponentiation(R, D, N_A);
  //	cout << R << endl;

  //	char t[] = "Bob_RSA_Public_Key";
  ifstream Public_File(t);

  Public_File >> N_B >> E;
  Public_File.close();

  R = Modular_Exponentiation(R, E, N_B);

  //	cout << R << endl;
  cout << "Successful!" << endl;
  cout << "Your digital signature for the receiver has been output to file "
          "\"Signature\""
       << endl;
  ofstream fout("Signature");

  fout << R << endl;
  fout.close();

  return 0;
}
