#include "object_partition.h"

using namespace std;

int main(int argc, char** argv) {
    if (argc != 2) {
        cerr<<"argc must be 2"<<endl;
        return -1;
    }

    string file = argv[1];
    object_partition partition = object_partition(file);
    partition.data_analysis();
    return 0;
}