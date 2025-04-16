#include <iostream>
#include <string>
#include <vector>
#include <cstring>
#include <fstream>
#include <unordered_map>

std::string filename;
std::string out_filename;
std::string seq_name;

// Structure to hold ORF details (protein added later for longest ORF)
struct ORF {
    size_t start_pos;
    size_t end_pos;
    size_t length;
    std::string sequence;
};

bool isStopCodon(const std::string &codon) {
    return codon == "TAA" || codon == "TAG" || codon == "TGA";
}

// Standard codon table
std::unordered_map<std::string, char> codonTable = {
    {"TTT", 'F'}, {"TTC", 'F'}, {"TTA", 'L'}, {"TTG", 'L'},
    {"TCT", 'S'}, {"TCC", 'S'}, {"TCA", 'S'}, {"TCG", 'S'},
    {"TAT", 'Y'}, {"TAC", 'Y'}, {"TAA", '*'}, {"TAG", '*'},
    {"TGT", 'C'}, {"TGC", 'C'}, {"TGA", '*'}, {"TGG", 'W'},
    {"CTT", 'L'}, {"CTC", 'L'}, {"CTA", 'L'}, {"CTG", 'L'},
    {"CCT", 'P'}, {"CCC", 'P'}, {"CCA", 'P'}, {"CCG", 'P'},
    {"CAT", 'H'}, {"CAC", 'H'}, {"CAA", 'Q'}, {"CAG", 'Q'},
    {"CGT", 'R'}, {"CGC", 'R'}, {"CGA", 'R'}, {"CGG", 'R'},
    {"ATT", 'I'}, {"ATC", 'I'}, {"ATA", 'I'}, {"ATG", 'M'},
    {"ACT", 'T'}, {"ACC", 'T'}, {"ACA", 'T'}, {"ACG", 'T'},
    {"AAT", 'N'}, {"AAC", 'N'}, {"AAA", 'K'}, {"AAG", 'K'},
    {"AGT", 'S'}, {"AGC", 'S'}, {"AGA", 'R'}, {"AGG", 'R'},
    {"GTT", 'V'}, {"GTC", 'V'}, {"GTA", 'V'}, {"GTG", 'V'},
    {"GCT", 'A'}, {"GCC", 'A'}, {"GCA", 'A'}, {"GCG", 'A'},
    {"GAT", 'D'}, {"GAC", 'D'}, {"GAA", 'E'}, {"GAG", 'E'},
    {"GGT", 'G'}, {"GGC", 'G'}, {"GGA", 'G'}, {"GGG", 'G'}
};

// Translate an ORF sequence to its protein sequence
std::string translateORF(const std::string &orf_sequence) {
    std::string protein;
    for (size_t i = 0; i + 3 <= orf_sequence.length(); i += 3) {
        std::string codon = orf_sequence.substr(i, 3);
        if (isStopCodon(codon)) {
            break; // Stop at the stop codon, don’t include it
        }
        protein += codonTable[codon];
    }
    return protein;
}

// Find all ORFs in a given frame and store them in a vector (no translation yet)
void findORFs(const std::string &sequence, int frame, std::vector<ORF> &orfs) {
    size_t start_pos = std::string::npos;

    for (size_t i = frame; i + 3 <= sequence.length(); i += 3) {
        std::string codon = sequence.substr(i, 3);

        if (start_pos == std::string::npos && codon == "ATG") {
            start_pos = i;
        } 
        else if (start_pos != std::string::npos && isStopCodon(codon)) {
            size_t end_pos = i + 3;
            std::string orf_seq = sequence.substr(start_pos, end_pos - start_pos);
            orfs.push_back({start_pos + 1, end_pos, orf_seq.length(), orf_seq});
            start_pos = std::string::npos; // Reset for next ORF
        }
    }
}

// Search all frames, select the longest ORF, and translate it
void searchORFs(const std::string &sequence, std::ofstream &outfile, std::string seq_name) {
    std::vector<ORF> all_orfs;

    // Check all three reading frames
    for (int frame = 0; frame < 3; frame++) {
        findORFs(sequence, frame, all_orfs);
    }

    // Find the longest ORF and translate it
    if (!all_orfs.empty()) {
        ORF longest_orf = all_orfs[0];
        for (const auto &orf : all_orfs) {
            if (orf.length > longest_orf.length) {
                longest_orf = orf;
            }
        }

        // Translate only the longest ORF
        std::string protein_seq = translateORF(longest_orf.sequence);

        // Write the longest ORF and its protein sequence to the output file
        //outfile << seq_name << "\t" << longest_orf.start_pos << "\t" << longest_orf.end_pos 
          //      << "\t" << longest_orf.length << "\t" << longest_orf.sequence 
            //    << "\t" << protein_seq << "\n";
                
        outfile << ">" << seq_name << "\n" << protein_seq << "\n";
    }
}

void print_info() {
    std::cout << "\n";
    std::cout << "***** Longest ORF finder with Protein Translation *****" << "\n";
    std::cout << "Input requires the input file argument:" << "\n";
    std::cout << "  -i    input file" << "\n";
    std::cout << "\n";
}

int main(int argc, char **argv) {
    if (argc == 1) {
        print_info();
        return 1;
    }

    for (int i = 1; i < argc; i += 2) {
        if (strcmp(argv[i], "-i") == 0 && i + 1 < argc) {
            filename = argv[i + 1];
            std::cout << "filename: " << filename << "\n";
        } else {
            print_info();
            return 1;
        }
    }
    
    std::ifstream infile(filename);
    out_filename = "Longest_ORFs_translated.fasta";
    std::ofstream outfile(out_filename);
  
    if (!infile.is_open() || !outfile.is_open()) {
        std::cerr << "Error: Unable to open file(s)!" << std::endl;
        return 1;
    }
    
    std::string sequence;
    std::string line;

    while (std::getline(infile, line)) {
        if (!line.empty() && line[0] == '>') {
            // Process the previous sequence if it exists
            if (!sequence.empty()) {
                searchORFs(sequence, outfile, seq_name);
                sequence.clear(); // Reset for the next sequence
            }
            seq_name = line.substr(1); // Store the new sequence name
        } else if (!line.empty()) {
            sequence += line; // Append to the current sequence
        }
    }

    // Process the last sequence
    if (!sequence.empty()) {
        searchORFs(sequence, outfile, seq_name);
    }
    
    infile.close();
    outfile.close();
    return 0;
}
