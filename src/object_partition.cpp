//
// Created by eason on 2/18/21.
//

#include <sstream>
#include <set>
#include "../include/object_partition.h"

//bool if_count_intra_redundant_collision = false;
object_partition::object_partition(string blastfile){
    fstream result_file(blastfile,ios::in);
    if (result_file.fail()) {
        cerr << "fail to open blastfile:" << blastfile << "!\n";
    }
    string line;
    while(getline(result_file,line)){
        if (line.size()<=1 || line[0]== '#')
            continue;

        string delimiter = "\t";
        string primer, object;
        primer = line.substr(0, line.find(delimiter));
        line.erase(0, line.find(delimiter) + delimiter.length());
        object = line.substr(0, line.find(delimiter));

        //screen out irrelevant strands
        if(g_if_control_payload_totalsize){
            int objectID = stoi(object.substr(7,object.size()));
            if(objectID>g_threshold_of_totalsize*objects_.size()) continue;
        }
        int objectID = stoi(object.substr(7,object.size()));
        object = to_string((objectID/g_chunk_size)*g_chunk_size);


        auto primer_it = primers_.find(primer);
        auto object_it = objects_.find(object);
        if(primer_it==primers_.end()) {
            PrimerRef *p;
            p=new PrimerRef;
            primers_.emplace(primer,p);
            primer_it = primers_.find(primer);
        }

        if(object_it==objects_.end()) {
            ObjectRef *o;
            o=new ObjectRef;
            objects_.emplace(object,o);
            object_it = objects_.find(object);
        }

        /*if(object_it->second->internal_collided_primers.find(primer)!=object_it->second->internal_collided_primers.end()){
                object_it->second->redundant_collision++;
                primer_it->second->redundant_collision++;
            continue;
        }*/

        primer_it->second->internal_collided_objects.emplace(object,object_it->second);
        object_it->second->internal_collided_primers.emplace(primer,primer_it->second);
    }
    // blast result only has collided primer & objects, we complete other free objects here
    for (int i = 1; i < g_total_strand_number; ++i) {
        auto ob_it = objects_.find(to_string(i));
        if (ob_it==objects_.end()){
            ObjectRef *o;
            o=new ObjectRef;
            objects_.emplace(to_string(i),o);
        }
    }
    //initialize K and S
    K_=(g_total_strand_number+primers_.size()*g_primer_capacity)/(g_tube_capacity*g_primer_capacity)+1;
    S_=(K_*g_tube_capacity*g_primer_capacity-g_total_strand_number-primers_.size()*g_primer_capacity)/g_primer_capacity;
    // initialize K tubes
    for (int i = 0; i < K_; ++i) {
        Tube tube;
        tubes_.push_back(tube);
    }
}

// assign single object to the tube, check whether has enough space to store
/*bool Tube::assign_object(ObjectRef *of) {
    for(auto n:of->internal_collided_primers){
        if (tube.removed_primer.find(n.first)==tube.removed_primer.end()){
            if(tube.available_primers>0){
                tube.available_primers--;
                tube.removed_primer.emplace(n.first);
                tube.free_space_in_current_primer--;
            } else{

            }
        }
    }
}*/
/*
void object_partition::baseline(){
    cout<<"collided primer number:"<<primers_.size()<<endl;
    cout<<"collided object number:"<<objects_.size()<<endl;

    for (int i = 1; i < objects_.size(); ++i) {
        auto ob_it = objects_.find(to_string(i));
        if (ob_it==objects_.end()) cout<<"[baseline], doesn't find object "<<i<<endl;
        bool if_assigned = false;
        for (int j = 0; j < tubes_.size(); ++j) {
            if(tubes_[j].assign_object(ob_it->second)){
                if_assigned = true;
                break;
            }
        }
        if (!if_assigned){
            cout<<"[baseline] current tube number is "<<tubes_.size()<<" ,which is not enough. Add new one."<<endl;
            Tube tube;
            if(!tube.assign_object(ob_it->second)) cout<<"[baseline] fail to assign object to the new added tube"<<endl;
            tubes_.push_back(tube);
        }
    }
}

void object_partition::decomposition_primer_graph() {}*/

void object_partition::data_analysis() {
    cout<<"collided primer number:"<<primers_.size()<<endl;
    cout<<"collided object number:"<<objects_.size()<<endl;
    int total_collision = 0;
    int intra_redundant_collision = 0;
    //print primer violation distribution
    set<string> hot_primers_collided_objects;
    int hot_primer = 0;

    vector<int> primer_distribution;
    for(auto i:primers_){
        // primers with different collision number
        int x_axis = i.second->internal_collided_objects.size()-1;
        /*if (if_count_intra_redundant_collision){
            x_axis+=i.second->redundant_collision;
        }*/
        if(primer_distribution.size()<x_axis+1){
            for (int j = primer_distribution.size(); j <= x_axis+1; ++j) {
                primer_distribution.push_back(0);
            }
        }
        primer_distribution[x_axis]++;

        //check hot primer's collided objects
        if(i.second->internal_collided_objects.size()>400){
            hot_primer++;
            for(auto n:i.second->internal_collided_objects){
                hot_primers_collided_objects.emplace(n.first);
            }
        }
    }

    ofstream myfile;
    myfile.open ("primer_distribution.csv",ios::out | ios::trunc);
    for(int i=0; i < primer_distribution.size(); i++){
        myfile<<i+1<<","<<primer_distribution[i]<<endl;
    }
    myfile.close();
    vector<int>().swap(primer_distribution); //release memory

    //print object violation distribution
    vector<int> object_distribution;
    for(auto i:objects_){
        //intra_redundant_collision+=i.second->redundant_collision;
        total_collision+=i.second->internal_collided_primers.size();
        // objects with different collision number
        int x_axis = i.second->internal_collided_primers.size();
        /*if (if_count_intra_redundant_collision){
            x_axis+=i.second->redundant_collision;
        }*/
        if(object_distribution.size()<x_axis+1){
            for (int j = object_distribution.size(); j <= x_axis+1; ++j) {
                object_distribution.push_back(0);
            }
        }
        object_distribution[x_axis]++;
    }


    myfile.open ("object_distribution.csv",ios::out | ios::trunc);
    for(int i=0; i < object_distribution.size(); i++){
        myfile<<i+1<<","<<object_distribution[i]<<endl;
    }
    myfile.close();
    vector<int>().swap(object_distribution); //release memory

    total_collision+=intra_redundant_collision;
    cout<<"redundant collision:"<<intra_redundant_collision<<" total_collision:"<<total_collision<<endl;


    cout<<"number of hot primer:"<<hot_primer<<" number of their collided objects:"<<hot_primers_collided_objects.size()<<endl;
}