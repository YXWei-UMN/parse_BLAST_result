//
// Created by eason on 2/18/21.
//

#include <sstream>
#include "object_partition.h"
bool control_payload_totalsize = false;
int threshold_of_totalsize = 0.01 * 1813257;
bool control_payload_singlesize = true;
int number_of_strands_for_one_object = 2;
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
        if(control_payload_totalsize){
            int objectID = stoi(object.substr(7,object.size()));
            if(objectID>threshold_of_totalsize) continue;
        }
        if(control_payload_singlesize){
            int objectID = stoi(object.substr(7,object.size()));
            cout<<objectID<<" ";
            object = to_string((objectID/number_of_strands_for_one_object)*number_of_strands_for_one_object);
            cout<<object<<endl;
        }

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

        primer_it->second->internal_collided_objects.emplace(object,object_it->second);
        object_it->second->internal_collided_primers.emplace(primer,primer_it->second);
    }

}


void object_partition::data_analysis() {
    cout<<"primer number:"<<primers_.size()<<endl;
    cout<<"object number:"<<objects_.size()<<endl;

    //print primer violation distribution
    vector<int> primer_distribution;
    for(auto i:primers_){
        // primers with different collision number
        if(primer_distribution.size()<i.second->internal_collided_objects.size()){
            for (int j = primer_distribution.size(); j <= i.second->internal_collided_objects.size(); ++j) {
                primer_distribution.push_back(0);
            }
        }
        primer_distribution[i.second->internal_collided_objects.size()-1]++;
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
        // primers with different collision number
        if(object_distribution.size()<i.second->internal_collided_primers.size()){
            for (int j = object_distribution.size(); j <= i.second->internal_collided_primers.size(); ++j) {
                object_distribution.push_back(0);
            }
        }
        object_distribution[i.second->internal_collided_primers.size()-1]++;
    }


    myfile.open ("object_distribution.csv",ios::out | ios::trunc);
    for(int i=0; i < object_distribution.size(); i++){
        myfile<<i+1<<","<<object_distribution[i]<<endl;
    }
    myfile.close();
    vector<int>().swap(object_distribution); //release memory

}