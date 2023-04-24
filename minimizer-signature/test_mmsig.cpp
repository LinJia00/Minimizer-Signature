#include <iostream>
#include <fstream>
#include <bitset>
#include <map>
#include <cmath>
#include <chrono>
using namespace std;
using namespace std::chrono;
//kmer size = 5
const uint64_t SIG_SIZE = 1025;
//kmer size = 6
// const uint64_t SIG_SIZE = 4097;
//kmer size = 7
// const uint64_t SIG_SIZE = 16385;
//kmer size = 8
// const uint64_t SIG_SIZE = 65537;
double elapsed(high_resolution_clock::time_point t1, high_resolution_clock::time_point t2) {
	return (duration_cast<duration<double>>(t2 - t1)).count();
}

float calc_jaccard(std::bitset<SIG_SIZE> sig1, std::bitset<SIG_SIZE> sig2){
    return (sig1 & sig2).count() * 1.0 / (sig1 | sig2).count();
}

void compare_sig_file(string test_file, string ref_file1, string ref_file2,string ref_file3, ofstream &outputFile){
    ifstream test(test_file);
    if(!test.good()){
        cerr << "Error opening "<< test_file << std::endl;
        
    }

    
    uint counter = 0;
    std::string line0;
    while(getline(test, line0)){
        if(line0[0] == '>'){
            outputFile << line0 <<endl;
        }else{
            std::bitset<SIG_SIZE> test_signature(line0);
            // comparing 1st ref
            std::string line1;
            double sum_similarity_1 = 0;
            uint64_t count_1 = 0;
            
            ifstream ref1(ref_file1);
            if(!ref1.good()){
                cerr << "Error opening " << ref_file1 << std::endl;
                
            }

            while (getline(ref1, line1))
            {
                std::bitset<SIG_SIZE> ref_signature1(line1);
                float temp = calc_jaccard(test_signature, ref_signature1);
                sum_similarity_1 += temp;
                count_1 += 1;
            
            }
            outputFile<< "P.1: " << sum_similarity_1/count_1 << endl;
            // comparing 2nd ref
            std::string line2;
            double sum_similarity_2 = 0;
            uint64_t count_2 = 0;
            
            ifstream ref2(ref_file2);
            if(!ref2.good()){
                cerr << "Error opening " << ref_file2 << std::endl;
                
            }

            while (getline(ref2, line2))
            {
                std::bitset<SIG_SIZE> ref_signature2(line2);
                float temp = calc_jaccard(test_signature, ref_signature2);
                sum_similarity_2 += temp;
                count_2 += 1;
            }
            outputFile<< "P.1.427: " << sum_similarity_2/count_2 << endl;

             // comparing 3rd ref
            std::string line3;
            double sum_similarity_3 = 0;
            uint64_t count_3 = 0;
            
            ifstream ref3(ref_file3);
            if(!ref3.good()){
                cerr << "Error opening " << ref_file3 << std::endl;
                
            }

            while (getline(ref3, line3))
            {
                std::bitset<SIG_SIZE> ref_signature3(line3);
                float temp = calc_jaccard(test_signature, ref_signature3);
                sum_similarity_3 += temp;
                count_3 += 1;
            }
            outputFile<< "P.1.526: " << sum_similarity_3/count_3 << endl;
            outputFile<< " " << endl;
        }
            
        
    }    
    // return sum_similarity/count; 
}


int main(int argc, char **argv)
{

     if (argc <= 3)
    {
        cerr << "Usage: " << argv[0] << " [in_test_mmsigs] [in_ref_mmsigs] [outfile]" << std::endl;

        return -1;
    }
    ifstream input0(argv[1]);
    if (!input0.good())
    {
        cerr << "Error opening " << argv[1] << endl;
        return -1;
    }
    ifstream input1(argv[2]);
    if (!input1.good())
    {
        cerr << "Error opening " << argv[2] << endl;
        return -1;
    }
    ifstream input2(argv[3]);
    if (!input1.good())
    {
        cerr << "Error opening " << argv[3] << endl;
        return -1;
    }
    ifstream input3(argv[4]);
    if (!input1.good())
    {
        cerr << "Error opening " << argv[4] << endl;
        return -1;
    }
    ofstream output(argv[5]);
    if (!output.good())
    {
        cerr << "Error opening " << argv[3] << endl;
        return -1;
    }
    high_resolution_clock::time_point t1, t2;
    t1 = high_resolution_clock::now();
    compare_sig_file(argv[1], argv[2], argv[3], argv[4], output);
    t2 = high_resolution_clock::now();
    std::cout << "Time to Query: " + std::to_string(elapsed(t1, t2)) + " secs\n";
}