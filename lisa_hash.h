/*************************************************************************************
MIT License

Copyright (c) 2020 Intel Labs

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

Authors: Saurabh Kalikar <saurabh.kalikar@intel.com>; Sanchit Misra <sanchit.misra@intel.com>;
*****************************************************************************************/

#include<map>
#include<vector>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include "RMI.h"
using namespace std;

template<typename rmi_key_t, typename rmi_val_t>
class lisa_hash{
	public:

		rmi_val_t* p;
		rmi_val_t* p_bin;

	private:

		rmi_key_t* keys;

		uint64_t* values_enc;
		uint64_t* values_enc_bin;
		
		uint64_t keys_size;
		uint64_t p_size;

 		RMI<rmi_key_t> *rmi;

		void mem_alloc(uint64_t hash_size, uint64_t p_size){

		// Removed for memory optimization
		#if 0
			keys = (rmi_key_t*)malloc((1+hash_size)*sizeof(rmi_key_t));
			keys = keys + 1; // keys[-1] stores total number of keys
		#endif
			values_enc = (uint64_t*)malloc(hash_size*sizeof(uint64_t));
			
			p = (rmi_val_t*)malloc(p_size*sizeof(rmi_val_t));
			
			//values_enc_bin = (uint64_t*)malloc(hash_size*sizeof(uint64_t));
			
			//p_bin = (rmi_val_t*)malloc(p_size*sizeof(rmi_val_t));

			keys_size = hash_size;
		//	keys[-1] = keys_size;	
			this->p_size = p_size;

			fprintf(stderr, "Memory allocated %lld \n", p_size);

		}
		//This function is no longer used. 	
		void store_sorted_keys(string binFileName){

				
			ofstream wf(binFileName, ios::out | ios::binary);
			wf.write((char*)&keys[-1], (1+keys_size)*sizeof(rmi_key_t));			      

			wf.close();	
			
			rmi_key_t *temp_keys = (rmi_key_t*)malloc((1 + keys_size)*sizeof(rmi_key_t));
			ifstream rf(binFileName, ios::out | ios::binary);
			rf.read((char*)&temp_keys[0], (1 + keys_size)*sizeof(rmi_key_t));

			if(temp_keys[0] == keys[-1]){
				for(int i = 0; i <  keys_size; i++){
					if(temp_keys[i+1] != keys[i]){
						cout<<"Error: "<<temp_keys[i+1]<< " " << keys[i]<<endl;
						break;
					}
				}
			}
			else{
				cout<<"File writing error!!\n"<<temp_keys[0] << " "<<keys[-1] << keys_size;
			}
			free(temp_keys);			      
		}	

		void set_key_value(uint64_t i, rmi_key_t key, uint64_t p_pos, uint64_t p_size){
		
		//		keys[i] = key;//<<1;
				values_enc[i] = p_pos << 32 | p_size;
		}
		
		void get_val(uint64_t i, uint64_t &pos, uint64_t &p_size){
			
			pos = values_enc[i] >> 32;
			p_size = (uint32_t) values_enc[i];
		}

		

		void load_bin(string inputFile){
			fprintf(stderr, "Loading from bin\n");

			string f1_name = (string) inputFile + "_pos_bin";
			string f2_name = (string) inputFile + "_val_bin";
			ifstream instream_f1(f1_name, ifstream::binary);
			ifstream instream_f2(f2_name, ifstream::binary);
			instream_f1.seekg(0);
			instream_f1.read((char*)&values_enc[0], keys_size*sizeof(uint64_t));
			instream_f1.close();
			
			instream_f2.seekg(0);
			instream_f2.read((char*)&p[0], p_size*sizeof(uint64_t));
			instream_f2.close();
			/*
			for(uint64_t i = 0; i < p_size; i++){
				if(p[i] != p_bin[i]){
					fprintf(stderr, "Error!! %lld %lld\n",p[i], p_bin[i] );
				}
			}	
			for(uint64_t i = 0; i < keys_size; i++){
				if(values_enc[i] != values_enc_bin[i]){
					fprintf(stderr, "Error!! %lld %lld\n", values_enc[i], values_enc_bin[i]);
				}
			}
			*/	
		}	
	
		void load(string inputFile){

			ifstream f(inputFile);
			
			rmi_key_t key;
			rmi_val_t val;
			uint64_t n;
			int64_t offset = 0;
			int64_t i = 0;
			while(f>>key){
				f>>n;
				//values_enc[++key_size] = val_size << 32 | n;
				for(int j = 0; j < n; j++){
					f>>val;
					p[offset + j] = val;	
				}
				set_key_value(i, key, offset, n);
				offset +=n;
				i++;
			}
		}	
		


	public:
		lisa_hash(string inputFile, char* rmi_prefix, long leaf = 0){
			uint64_t start_time, load_time, rmi_building_time, rmi_object;
			start_time = __rdtsc();
			
			uint64_t val_count = 0;
			
			ifstream f_size(inputFile+"_size");
			f_size>>keys_size;
			f_size>>p_size;
			mem_alloc(keys_size, p_size);
	
			fprintf(stderr, "Num_keys: %lld, num_values = %lld", keys_size, p_size);
			//load(inputFile);
			load_bin(inputFile);
		
			string prefix = inputFile + "_keys";
#ifdef UINT64	
			string keys_bin_file_name = prefix + ".uint64";		
#else
			string keys_bin_file_name = prefix + ".f64";		
#endif	
	
			ifstream rf(keys_bin_file_name, ios::out | ios::binary);
			if(!rf.good()){
				cout<<"Error: Binary file with keys not found!!\n";
				//call store keys
				//store_sorted_keys(keys_bin_file_name);
				exit(0);
			}


			string keys_rmi_file = prefix + ".rmi_PARAMETERS"; 
			ifstream rmi_f(keys_rmi_file, ios::out | ios::binary);
			if(!rmi_f.good() || leaf != 0){
				cout<<"rmi file not found: "<< keys_rmi_file <<endl;
				if(leaf == 0)
				{
					cout<<"Number of rmi leaf nodes are not provided for "<<keys_size<<" keys\n";
					cout<<"Using default number of leaf nodes: "<<keys_size/32<<"\n";
					leaf = keys_size/32;
				}							
				load_time = __rdtsc() - start_time;
				start_time = __rdtsc();
				fprintf(stderr,"TIMER LOG: kay-val files loading time- %lld\n", load_time);

				//string script = "./scripts/build_rmi.sh";
				string script = "./scripts/build-rmi.linear_spline.linear.sh";
				string keys_path = " " + keys_bin_file_name;
				string rmi_path = " " + keys_rmi_file;
				string num_leaf = " " + std::to_string(leaf); 
#ifdef UINT64	
				string key_type = " UINT64";
#else
				string key_type = " F64";
#endif	
				string cmd = script + keys_path + rmi_path + num_leaf + key_type;				
				system(cmd.c_str());
			        rmi_building_time = __rdtsc() - start_time;
				fprintf(stderr,"TIMER LOG: rmi building time- %lld\n", rmi_building_time);
			}
			rmi = new RMI<rmi_key_t>(&prefix[0]);
			

		}		
		

		rmi_val_t* get_hash_value(rmi_key_t key, int*n){
			uint64_t index = rmi->lookup(key);
			if(index == -1){
				cout<<"Key not found\n";
				return NULL;
			}
			uint64_t pos, p_size;

			get_val(index, pos, p_size);	
			*n = p_size;
			return &p[pos];
		}

		rmi_val_t* get_hash_values_batched(rmi_key_t *keys, uint64_t num_keys, int* &num_values){

			int64_t *pos = (int64_t*) malloc(num_keys*sizeof(int64_t));	
			
			rmi->lookup_batched(keys, num_keys, &pos[0]);


			num_values = (int*) malloc(num_keys*sizeof(int));
			rmi_val_t **p_ptrs = (rmi_val_t**) malloc(num_keys*sizeof(rmi_val_t*));	
			uint64_t total_num_values = 0;
			for(int i = 0; i < num_keys; i++){
				if(pos[i] == -1){
					num_values[i] = 0;
					continue;
				}				
				num_values[i] = (uint32_t) values_enc[pos[i]];
				total_num_values+= num_values[i];
				p_ptrs[i] = &p[values_enc[pos[i]] >> 32];
			}	
			free(pos);
			rmi_val_t* ret_values = (rmi_val_t*) malloc(total_num_values*sizeof(rmi_val_t));	
			
			uint64_t cnt= 0;
			for(int i = 0; i < num_keys; i++){

				rmi_val_t *value_start_ptr = p_ptrs[i];
				uint64_t numhit = num_values[i];

 
				for(int j = 0; j < numhit; j++){
					ret_values[cnt++] = value_start_ptr[j];
				}
			}
			free(p_ptrs);
			return ret_values;
		}


		void mm_idx_get_batched(uint64_t* &minimizers, uint64_t num_minimizers, int64_t* &pos, uint64_t** &p_ptrs, int* &num_hits){
			
			rmi->lookup_batched(minimizers, num_minimizers, &pos[0]);
			
			for(int i = 0; i < num_minimizers; i++){
				int64_t p_i = pos[i];

				if(p_i < 0 || p_i > keys_size )
				num_hits[i] = 0;
				else
				num_hits[i] = (uint32_t) values_enc[p_i];
				p_ptrs[i] = p + (values_enc[p_i] >> 32);
			}	
		}
		~lisa_hash(){
			delete rmi;
			free(values_enc);
			free(p);
		}

};
