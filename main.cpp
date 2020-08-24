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
    string last = "#";
    while(getline(result_file,line)){
        size_t pos_last = last.find("#");
        size_t pos_now = line.find("#");
        if(pos_now == line.npos && pos_last != last.npos) hit_num++;
        else if(pos_last == last.npos && pos_now != line.npos) hit_num++;
    }

    cout<<"violated primers: "<<hit_num/2<<endl;
    return 0;
}