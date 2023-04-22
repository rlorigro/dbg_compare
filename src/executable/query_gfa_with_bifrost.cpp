#include "bifrost/include/bifrost/CompactedDBG.hpp"
#include "Filesystem.hpp"
#include "CLI11.hpp"
#include <iostream>
#include <fstream>

#include <functional>

using std::ifstream;
using std::function;
using std::getline;
using ghc::filesystem::path;


void for_each_sequence_in_fasta(path fasta_path, const function<void(const string& name, const string& sequence)>& f){
    ifstream file(fasta_path);

    if (not (file.is_open() and file.good())){
        throw runtime_error("ERROR: could not read file: " + fasta_path.string());
    }

    string name;
    string sequence;
    string line;
    size_t l = 0;

    while (getline(file, line)){
        cout << line << '\n';
        if (line[0] == '>'){
            if (l > 0){
                f(name, sequence);
            }

            // Trim any trailing tokens from the fasta header, keep only the name
            name = line.substr(1, line.find_first_of(" \t\n") - 1);
            cout << name << '\n';

            // Reset sequence
            sequence.clear();
        }
        else {
            sequence += line;
            cout << sequence << '\n';
        }

        l++;
    }

    f(name,sequence);
}


int query_gfa_with_bifrost(path gfa_path, path fasta_path, size_t n_threads, path output_dir){
    if (exists(output_dir)){
        throw runtime_error("ERROR: output directory exists already");
    }
    else {
        create_directories(output_dir);
    }

    const size_t k = 31;

    CompactedDBG<> cdbg(k);

    const bool verbose = true;

    cdbg.read(gfa_path.string(), n_threads, verbose);

    cout << "K-mer size is " << cdbg.getK() << endl;

    for_each_sequence_in_fasta(fasta_path, [&](const string& name, const string& sequence){
        cout << "---- " << name << " ----" << '\n';

        for (KmerIterator it_km(sequence.c_str()), it_km_end; it_km != it_km_end; ++it_km) { //non-ACGT char. are discarded
            UnitigMap<> um = cdbg.find(it_km->first);

            if (um.isEmpty) {
                cout << "Kmer " << it_km->first.toString() << " was not found" << endl;
            }
            else {
                const string unitig = um.referenceUnitigToString();
                const string strandness = um.strand ? "forward" : "reverse-complement";
                const size_t position = um.dist;

                cout << "Kmer " << it_km->first.toString() << " was found in the " << strandness << " direction of unitig " << "" << "at position " << position << '\n';
            }
        }
    });

    return 0;
}


int main(int argc, char* argv[]){
    path gfa_path;
    path fasta_path;
    path output_dir;
    size_t n_threads = 1;

    CLI::App app{"App description"};

    app.add_option(
            "-g,--gfa",
            gfa_path,
            "Path to GFA")
            ->required();

    app.add_option(
            "-q,--query_fasta",
            fasta_path,
            "Path to Fasta containing sequences to query")
            ->required();

    app.add_option(
            "-o,--output_dir",
            output_dir,
            "Path to (nonexistent) directory where output will be stored")
            ->required();

    app.add_option(
            "-t,--threads",
            n_threads,
            "(Default = " + to_string(n_threads) + ")\tMaximum number of threads to use.");

    CLI11_PARSE(app, argc, argv);

    query_gfa_with_bifrost(gfa_path, fasta_path, n_threads, output_dir);

    return 0;
}