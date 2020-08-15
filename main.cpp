#include <iostream>
#include <fstream>

using namespace std;
int main(int argc, char** argv) {
    if (argc != 2) {
        cerr<<"argc must be 2"<<endl;
        return -1;
    }

    string file = argv[1];

    fstream result_file;
    result_file.open(file,ios::in);
    string line;
    size_t hit_num = 0;

    while(getline(result_file,line)){
        size_t pos = line.find("#");
        if(pos != line.npos) hit_num++;
    }

    cout<<"no # line number: "<<hit_num<<". please divide primer_num"<<endl;
    return 0;
}