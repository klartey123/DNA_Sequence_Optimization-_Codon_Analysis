// Author: Kelvin Lartey
// This is my first c++ personal project on bioinformatics and working on biological data
// This project, is to help soildy my skills in parsing and introducing me in coding with biological data

// So the purpose of this program is to parse a DNA sequence from a fasta file and perform multiple analyze on it
// The ana;yze is to help detemine the numer of codon frequency and the nucleotides bases. Then base on all that,
// the results of the codon will be returned base on the a scoring system. 

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <algorithm>
#include <map>
#include <unordered_map>
#include <queue>

using namespace std;

// Nucletiode counter. This file will count the Nucleotide. A,T,G,C

// GC content, we know that GC bonds, to account for the strength of the sequence. The more GC bonds the stronger the sturcture. 
// So this function will analyze the percent of GC content in the sequence. 

void analysisGC (const unordered_map <string, string>& sequence){
    for (auto& pair: sequence){
        string id = pair.first; // Id
        string dna = pair.second; // Dna

        int countA = 0, countT = 0, countG = 0, countC = 0;

        for (char base : dna) {
            if (base == 'A') countA++;
            else if (base == 'T') countT++;
            else if (base == 'G') countG++;
            else if (base == 'C') countC++;
        }

        int total = countA + countT + countG + countC;
        double gcContent = 0;

        // get in case to avaoid dividing by 0
        if (total > 0){
            gcContent = ((double)(countG + countC) / total) * 100.0;
        }
        
        cout << "Sequence: " << id << endl;
        cout << "A: " << countA << "  T: " << countT
             << "  G: " << countG << "  C: " << countC << endl;
        cout << "GC Content: " << gcContent << "%" << endl;
        cout << "---------------------------" << endl;

    }
}

// Now to count how many condons appear in a sequence. so a condon is Three base pairs such as AGT, CAA

unordered_map<string, int> codonFrq(const unordered_map<string, string>& sequence){
    unordered_map<string, int> count_codon;
    for (auto& pair: sequence){
        string id = pair.first; // Id
        string dna = pair.second; // Dna
        
        for(int i = 0; i + 2 < dna.length(); i += 3){
            string codon = dna.substr(i, 3);
            count_codon[codon]++;
            }
        }
        return count_codon;

    }

    // Two specific organism for rank score base on how those organisim codon appear in the file so e_coli and yeast are the exmple
    // This is a foundation the potentail of this project, instead of having two organisim, there can be 
    // a machine learning algorithm that take a huge data set of multiple organism and there seqenences and that can be
    // use in the scoring. 

    vector<string> egOrganisim(string organism) {
        if (organism == "E_coli") {
            return {"ATG", "GCT", "TGG"};
        } else if (organism == "Yeast") {
            return {"ATG", "GGA", "TGC"};
        }
        return {}; 
    }

    // this is dependent on the function above as this functon is only used if the individual decides to have an organism to 
    // compare there sequence too. 

    unordered_map<string, int> bonus_points(const unordered_map<string, string>& sequences, const vector<string>& preferred_codons) {
        unordered_map<string, int> codon_bonus;
    
        for (auto& pair : sequences) {
            string dna = pair.second;
            for (int i = 0; i + 2 < dna.length(); i += 3) {
                string codon = dna.substr(i, 3);
                bool is_preferred = false;
                for (const string& preferred : preferred_codons) {
                    if (codon == preferred) {
                        is_preferred = true;
                        break;
                    }
                }
                if (is_preferred) {
                    codon_bonus[codon] += 2;
                }
            }
        }
    
        return codon_bonus;
    }

    // Now we are analyzing the FASTA file to see how good this dna sequence is. This is the base scoring algorithm 
    // before an organisim is used as a refrence of scoring. 

    void basescore(unordered_map<string, string>& sequence, const vector<string>& preferred_codons, bool use_bonus){
        unordered_map<string, int> new_codoncount = codonFrq(sequence);
        unordered_map<string, int> codon_scores;
        unordered_map<string, int> codon_bonus_map;
    
        if (use_bonus) {
            codon_bonus_map = bonus_points(sequence, preferred_codons);
        }
    
        int total_codons = 0;
        for (auto& pair : new_codoncount) {
            total_codons += pair.second;
        }
    
        for (auto& pair: new_codoncount){
            string codon = pair.first;
            int count = pair.second;
    
            // checking if the codon is valid. 
            bool is_valid = true;
            for (char c : codon) {
                if (c != 'A' && c != 'T' && c != 'G' && c != 'C') {
                    is_valid = false;
                    break;
                }
            }
    
            if (!is_valid) {
                cout << "Codon " << codon << " is invalid. Skipping...\n";
                continue; 
            }
    
            // GC content score
            int gc = 0;
            int gc_score = 0;
            for (char base : codon) {
                if (base == 'G' || base == 'C') 
                    gc++;
            }
            double gc_content_weight = (gc / 3.0) * 100 * 0.02; 
            if (gc_content_weight > 0) {
                gc_score += 3;  
            }
    
            // Frequency percentage
            double freq_percent = (count / (double)total_codons) * 100;
            int freq_score = 0;
            if (freq_percent >= 11 && freq_percent <= 30) freq_score = 1;
            else if (freq_percent >= 31 && freq_percent <= 60) freq_score = 2;
            else if (freq_percent >= 61 && freq_percent <= 90) freq_score = 1;
    
            int total_score = freq_score + gc_score;
            if (use_bonus && codon_bonus_map.count(codon)) {
                total_score += codon_bonus_map[codon];
            }
    
            codon_scores[codon] = total_score;
        }
    
        //  using s priority queue as a source to get the biggest score.

        // Change this to fit more of what was done in class. For better optimization. 
        priority_queue<pair<int, string>> ranked_codons;
        for (const auto& pair : codon_scores) {
            ranked_codons.push({pair.second, pair.first});
        }
    
        cout << "Ranked Codons (Highest to Lowest): "<< endl ;
        while (!ranked_codons.empty()) {
            auto top = ranked_codons.top();
            cout << "Codon: " << top.second << " | Score: " << top.first << endl;
            ranked_codons.pop();
        }
    }
    


int main(){
    string file_name;
    cout << "Enter file name: ";
    cin >> file_name;
    
    ifstream inFS(file_name);
       if (!inFS.is_open()) {
           cout << "Could not open file: " << file_name << endl;
           return 1;
       }
    string line;
    string current_line;
    //getline(inFS, line);
    unordered_map<string , string> sequence;
    while (getline(inFS, line)) {
        if (line[0] == '>'){
            current_line = line.erase(0,1);
            sequence[current_line]= "";

        }
        else{
            sequence[current_line] += line;
        }     
    }

    analysisGC(sequence);
    string organism;
    cout << "Enter organism name for bonus scoring (or leave blank): ";
    // need sp that the organism name can be implemented or not. 
    cin.ignore();
    getline(cin, organism);

    vector<string> preferred_codons;
    bool use_bonus = false;
    if (!organism.empty()) {
        preferred_codons = egOrganisim(organism);
        use_bonus = true;
    }

    basescore(sequence, preferred_codons, use_bonus);
    if (inFS.is_open()) {
        inFS.close();
     }
    return 0;
}


  
