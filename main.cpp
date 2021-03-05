#include "include/object_partition.h"
#include "include/global.h"
using namespace std;

int main(int argc, char** argv) {
    if (argc != 2) {
        cerr<<"argc must be 2"<<endl;
        return -1;
    }


    string cfgfile = argv[1];

    if (Parse(cfgfile)) {
        cerr<< "parse config file " << cfgfile << " failed!\n";
        return -1;
    }


    object_partition partition = object_partition(g_blast_Result);
    partition.data_analysis();
    /*if (g_if_baseline)
        partition.baseline();
    else if (g_if_decomposition_on_primer_graph)
        partition.decomposition_primer_graph();*/
    return 0;
}