/*#include <iostream>
#include <fstream>


int main(int argc, char* argv[]){
	int size = 1024;
    char buffer[size];
    ifstream myFile ("Elapse_time0.txt", ios::in | ios::binary);
    myFile.read (buffer, size);
    if (!myFile) {
        cout<<"error 1 occured"<<endl;
        int count = myFile.gcount();		// returns the number of bytes read.
        cout<<count<<" elements read successfully"<<endl;
        // calling myFile.clear() will reset the stream state
        // so it is usable again.
    }
    if (!myFile.read (buffer, size)) {
    	cout<<"error 2 occured"<<endl;
    	int count = myFile.gcount();		// returns the number of bytes read.
    	cout<<count<<" elements read successfully"<<endl;
    }
}*/

// read a file into memory
#include <iostream>     // std::cout
#include <fstream>      // std::ifstream
using namespace std;

int main () {

  std::ifstream is ("Elapsed_time0.txt", std::ifstream::binary);
  if (is) {
    // get length of file:
    is.seekg (0, is.end);
    int length = is.tellg();
    is.seekg (0, is.beg);

    char * buffer = new char [length];

    std::cout << "Reading " << length << " characters... ";
    // read data as a block:
    is.read (buffer,length);

    if (is)
      std::cout << "all characters read successfully.";
    else
      std::cout << "error: only " << is.gcount() << " could be read";
    is.close();

    // ...buffer contains the entire file...
    for (int i=0; i<length; i++){
    	printf("%lf\n",buffer[i]);
    }
    delete[] buffer;
  }
  return 0;
}


