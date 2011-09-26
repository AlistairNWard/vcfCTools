#include <string>
#include <cstdlib>
#include <cstdio>
#include <iostream>

using namespace std;

int main(int argc, char* argv[]) {
  string reference = "GACGAGCCTGACTGACGTCAGTGCTCTCTCT";
  string refAllele = "CTCTC";
  string altAllele = "C";

  string deleted   = refAllele.substr(1, refAllele.size() - 1);

  size_t i = reference.find(deleted);

  cout << refAllele << endl << altAllele << endl << deleted << endl << i << endl;
}
