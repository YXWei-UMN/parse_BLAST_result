//
// Created by wyx on 19-6-20.
//

#include "../include/global.h"



string g_blast_Result;
int g_primer_capacity;
int g_tube_capacity;
long g_chunk_size;
bool g_if_baseline;
bool g_if_decomposition_on_primer_graph;
bool g_if_control_payload_totalsize;
double g_threshold_of_totalsize;
long long int g_total_strand_number;

int Parse(string cfgfile){
    ifstream filestream(cfgfile, ios_base::in);
    if (filestream.fail()) {
        cerr << "open cfgfile:" << cfgfile << " fails!\n";
        return -1;
    }
    string line;

    while(getline(filestream, line)) {
        if (line.size()<=1 || line[0]== '#')
            continue;

        stringstream ss(line);
        string key, value;
        getline(ss, key, ' ');
        getline(ss, value, ' ');

        switch(hash_(key.c_str())){
            case hash_("primer_capacity"):
                g_primer_capacity = stoi(value);
                break;
            case hash_("total_strand_number"):
                g_total_strand_number = stoll(value);
                break;
            case hash_("tube_capacity"):
                g_tube_capacity = stoi(value);
                break;
            case hash_("chunk_size"):
                g_chunk_size = stol(value);
                break;
            case hash_("blast_Result"):
                g_blast_Result = value;
                break;
            case hash_("if_baseline"):
                g_if_baseline = (value=="true");
                break;
            case hash_("if_decomposition_on_primer_graph"):
                g_if_decomposition_on_primer_graph = (value=="true");
                break;
            case hash_("if_control_payload_totalsize"):
                g_if_control_payload_totalsize = (value=="true");
                break;
            case hash_("threshold_of_totalsize"):
                g_threshold_of_totalsize = stod(value);
                break;
            default:
                cout<<"unknown cfg: "<<key<<endl;
                return -1;
        }
    }
    return 0;
}
