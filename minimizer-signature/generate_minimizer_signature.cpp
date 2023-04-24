/*
    References:
    https://stackoverflow.com/a/37952488

*/

#include <iostream>
#include <fstream>
#include <vector>
#include <deque>
#include <cstring>
#include <tuple>
#include <iterator>
#include <map>
#include <cmath>
#include <bitset>
#include <chrono>


using namespace std;
using namespace std::chrono;

// double elapsed(high_resolution_clock::time_point t1, high_resolution_clock::time_point t2) {
// 	return (duration_cast<duration<double> >(t2 - t1)).count();
// }

// int KMER_SIZE = 6;
// const uint64_t SIG_SIZE = pow(4,KMER_SIZE)+1;
//kmer size = 5
const uint64_t SIG_SIZE = 1025;
//kmer size = 6
// const uint64_t SIG_SIZE = 4097;
//kmer size = 7
// const uint64_t SIG_SIZE = 16385;
//kmer size = 8
// const uint64_t SIG_SIZE = 65537;
std::map<std::string, int> lex_map;

double elapsed(high_resolution_clock::time_point t1, high_resolution_clock::time_point t2) {
	return (duration_cast<duration<double>>(t2 - t1)).count();
}

void populate_map(uint kmer_size) {
    // By storing it in a file we save time required to generate all of the values
    // TODO: Or better to recursively generate it?
    std::string kmer_words_path = "../metadata/kmer_words_" + std::to_string(kmer_size) + ".txt";
    std::ifstream kmer_words_file(kmer_words_path);
    if(!kmer_words_file.good()){
        std::cerr << "Error opening "<< kmer_words_path << std::endl;
        exit(1);
    }

    std::string line;
    uint64_t counter = 0;
    while( std::getline( kmer_words_file, line ) ){
        lex_map[line] = counter; counter++;
    }
}


// std::bitset<SIG_SIZE> return_signature(std::string content){
//     std::bitset<SIG_SIZE> signature;
//     for (uint i=0; i<=content.size()-KMER_SIZE; i++){
//         std::string temp_string = content.substr(i, KMER_SIZE);
//         signature.set(lex_map[temp_string]);
//     }
//     // std::cout << signature << "\n";
//     return signature;
// }

void generate_minimizer_signature(string content, int32_t kmer_size, int32_t window_size, ofstream &outputFile)
{
    populate_map(kmer_size);
    
    int32_t last_start_idx = content.length() - kmer_size;
    deque<tuple<int32_t, string>> temp_mms;
    vector<int32_t> all_mms_idx;
    bitset<SIG_SIZE> signature;


    

    for (int start_idx = 0; start_idx < last_start_idx; start_idx++)
    {
        
       
        while (!temp_mms.empty() and (strcmp(get<1>(temp_mms.back()).c_str(), content.substr(start_idx, kmer_size).c_str()) > 0))
        {
            
           
          
            temp_mms.pop_back();
            
        }
        


        tuple<int32_t, string> inserted_item = make_tuple(start_idx, content.substr(start_idx, kmer_size));
        temp_mms.push_back(inserted_item);
   

        if ((start_idx - window_size) >= 0 && get<0>(temp_mms.front()) <= start_idx - window_size)
        {
            
            while (get<0>(temp_mms.front()) <= start_idx - window_size)
            {
        
                temp_mms.pop_front();
            }
        }


       
        // push back the mm idx
        if (!temp_mms.empty()){
            if( all_mms_idx.size() == 0 || all_mms_idx.back() != get<0>(temp_mms.front())){
              
                
                all_mms_idx.push_back(get<0>(temp_mms.front()));
                string current_minimizer = content.substr(get<0>(temp_mms.front()), kmer_size);
                signature.set(lex_map[current_minimizer]);
            }
        }
    }
    // copy(all_mms_idx.begin(), all_mms_idx.end(), ostream_iterator<int>(outputFile, " "));
    outputFile << signature <<"\n";
    
    

}

int main(int argc, char **argv)
{

    // high_resolution_clock::time_point t1, t2, t3, t4;
    if (argc <= 4)
    {
        std::cerr << "Usage: " << argv[0] << " [infile] [kmer_size] [window_size] [outfile]" << std::endl;

        return -1;
    }

    std::ifstream input(argv[1]);
    if (!input.good())
    {
        std::cerr << "Error opening " << argv[1] << std::endl;
        return -1;
    }
    const char *k(argv[2]);
    if (!input.good())
    {
        std::cerr << "Error opening " << argv[2] << std::endl;
        return -1;
    }
    const char *w(argv[3]);
    if (!input.good())
    {
        std::cerr << "Error opening " << argv[3] << std::endl;
        return -1;
    }
    std::ofstream output(argv[4]);
    if (!output.good())
    {
        std::cerr << "Error opening " << argv[4] << std::endl;
        return -1;
    }

    high_resolution_clock::time_point t1, t2;
    double time_gensig = 0;

    std::string line, name, content;
    while (std::getline(input, line))
    {

        if (line.empty() || line[0] == '>')
        { // Identifier marker
            // cout << name << "\n";
            if (!name.empty())
            { // Print out what we read from the last entry
                // Write signature
                output << '>' << name << "\n";
                t1 = high_resolution_clock::now();
                generate_minimizer_signature(content, atoi(k), atoi(w), output);
                t2 = high_resolution_clock::now();
                time_gensig += elapsed(t1, t2);
                name.clear();
            }
            if (!line.empty())
            {
                name = line.substr(1);
               
            }
            content.clear();
        }
        else if (!name.empty())
        {
            if (line.find(' ') != std::string::npos)
            { // Invalid sequence--no spaces allowed
                name.clear();
                content.clear();
            }
            else
            {
                content += line;
            }
        }
    }
    if (!name.empty())
    { // Print out what we read from the last entry
        output << '>' << name << "\n";
        t1 = high_resolution_clock::now();
        generate_minimizer_signature(content, atoi(k), atoi(w), output);
        t2 = high_resolution_clock::now();
        time_gensig += elapsed(t1, t2);
    }
    // double total_time_consumption = elapsed(t1, t2) + elapsed(t3, t4);
    std::cout << "Time to Generate All Minimizer Signatures: " + std::to_string(time_gensig) + " secs\n";
    return 0;
}