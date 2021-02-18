//
// Created by eason on 2/18/21.
//

#ifndef PARSE_BLAST_RESULT_OBJECT_PARTITION_H
#define PARSE_BLAST_RESULT_OBJECT_PARTITION_H
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <vector>
using namespace std;
struct ObjectRef;

struct PrimerRef {
    unordered_map<string, ObjectRef*> internal_collided_objects; // <chunk ID, chunkRef>
};

struct ObjectRef {
    unordered_map<string, PrimerRef*> internal_collided_primers; // <file name , fileRef>
};

class object_partition {
public:
    unordered_map<string , PrimerRef*> primers_;
    unordered_map<string , ObjectRef*> objects_;

    object_partition(string blastfile);
    ~object_partition(){

    }
    void data_analysis();
};


#endif //PARSE_BLAST_RESULT_OBJECT_PARTITION_H
