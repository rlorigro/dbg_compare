#include "bifrost/include/bifrost/CompactedDBG.hpp"

// Taken from Bifrost documentation:
// https://github.com/pmelsted/bifrost/blob/master/doc/tutorial/Intro.md#creating-a-graph

int main(){
    const size_t k = 31;

    CompactedDBG<> cdbg(k);

    cout << "K-mer size is " << cdbg.getK() << endl;

    const string seq = "ACGTCGTACGTCCCGTAAACGTTAAACGTAAACGTGTGTGCAAAATGTCTAGTTTTTTTACGCTGATATAGTC";

    cdbg.add(seq);

    const string kmer_sequence = "ACGTCGTACGTCCCGTAAACGTTAAACGTAA";
    const Kmer km = Kmer(kmer_sequence.c_str());

    UnitigMap<> um = cdbg.find(km);

    if (um.isEmpty) cout << "Kmer " << kmer_sequence << " was not found" << endl;
    else {

        const string unitig = um.referenceUnitigToString();
        const string strandness = um.strand ? "forward" : "reverse-complement";
        const size_t position = um.dist;

        cout << "Kmer " << kmer_sequence << " was found in the " << strandness << " direction of unitig " << unitig << "at position " << position << '\n';
    }

    return 0;
}

