//
// Created by eason on 2/18/21.
//

#ifndef PARSE_BLAST_RESULT_OBJECT_PARTITION_H
#define PARSE_BLAST_RESULT_OBJECT_PARTITION_H
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <vector>
#include "global.h"
using namespace std;
struct ObjectRef;

struct PrimerRef {
    //int redundant_collision=0;
    unordered_map<string, ObjectRef*> internal_collided_objects; // <chunk ID, chunkRef>
};

struct ObjectRef {
    //int redundant_collision=0;
    unordered_map<string, PrimerRef*> internal_collided_primers; // <file name , fileRef>
};

class Tube {
public:
    Tube(){
        available_primers= g_tube_capacity;
        free_space_in_current_primer = g_primer_capacity;
    }
    ~Tube(){}
    int available_primers;
    int free_space_in_current_primer;
    unordered_set<string> removed_primer; // number of removed primer + stored strand/primer capacity <= tube capacity

    bool assign_object(ObjectRef* of);
};

class object_partition {
public:
    unordered_map<string , PrimerRef*> primers_;
    unordered_map<string , ObjectRef*> objects_;
    vector<Tube> tubes_;
    // least number of tubes
    int K_;
    // most number of removal
    int S_;

    object_partition(string blastfile);
    ~object_partition(){

    }
    void data_analysis();
    // the baseline is sequentially store strands without distinguishing its collision
    // check how many tubes used
    void baseline();
    void decomposition_primer_graph();
};


#endif //PARSE_BLAST_RESULT_OBJECT_PARTITION_H
