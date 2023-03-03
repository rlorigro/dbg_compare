#include "bifrost/include/bifrost/CompactedDBG.hpp"

// Adapted from Bifrost documentation:
// https://github.com/pmelsted/bifrost/blob/master/doc/tutorial/Intro.md#creating-a-graph

int main(){
    const size_t k = 31;

    CompactedDBG<> cdbg(k);

    cout << "K-mer size is " << cdbg.getK() << endl;

    const string reference_seq = "ACGTCGTACGTCCCGTAAACGTTAAACGTAAACGTGTGTGCAAAATGTCTAGTTTTTTTACGCTGATATAGTC";
    const string query_sequence = "ACGTCGTACGTCCCGTAAACGTTAAACGTAAACGTGTGTGCAAAATGTCTAGTTTTTTTACGCTGATATAGTC";

    cdbg.add(reference_seq);

    for (KmerIterator it_km(query_sequence.c_str()), it_km_end; it_km != it_km_end; ++it_km) { //non-ACGT char. are discarded
        UnitigMap<> um = cdbg.find(it_km->first);

        if (um.isEmpty) cout << "Kmer " << it_km->first.toString() << " was not found" << endl;
        else {
            const string unitig = um.referenceUnitigToString();
            const string strandness = um.strand ? "forward" : "reverse-complement";
            const size_t position = um.dist;

            cout << "Kmer " << it_km->first.toString() << " was found in the " << strandness << " direction of unitig " << unitig << "at position " << position << '\n';
        }
    }

    return 0;
}

