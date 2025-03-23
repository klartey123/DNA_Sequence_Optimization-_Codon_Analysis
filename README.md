# DNA Codon Analysis and Optimization
**Author**: Kelvin Lartey

##  About the Project
This is my first C++ bioinformatics project. It parses DNA sequences from a FASTA file and performs codon frequency analysis, GC content calculation, and scoring based on biological properties.

The goal is to simulate real-world codon optimization. A crucial process in synthetic biology and gene expression studies.

## What It Does
- Parses DNA sequences from a FASTA-formatted file.
- Calculates nucleotide counts (A, T, G, C).
- Computes GC content percentage for each sequence( Determines the stability)
- Splits sequences into codons and calculates their frequency.
- Scores codons based on:
  - Frequency of appearance
  - GC content strength
  - Bonus points for organism-specific preferred codons (e.g., E. coli or Yeast)
- Ranks codons using a priority queue.

##  Example FASTA Input
```fasta
>sequence1
ATGGCTTGGATG
```

##  Sample Output
```
Sequence: sequence1
A: 2  T: 4  G: 5  C: 1
GC Content: 50%
---------------------
Ranked Codons (Highest to Lowest):
Codon: ATG | Score: 9
Codon: GCT | Score: 6
Codon: TGG | Score: 6
```

## How to Use It
1. Clone or download this project.
2. Compile the C++ file using g++ or any compiler:
   ```
   g++ -std=c++11 bioinformatics_parsing#1.cpp -o bioinfo
   ```
3. Run it:
   ```
   ./bioinfo
   ```
4. Enter the name of your FASTA file.
5. Optionally, enter `E_coli` or `Yeast` to include organism-specific scoring.

##  Future Ideas
- Export results to CSV
- Support multiple organisms from real codon usage tables
- Add a GUI or web interface
- Machine learning-based codon scoring

## Skills Demonstrated
- File parsing and FASTA handling
- Use of unordered maps, priority queues, and modular functions
- Custom scoring algorithms
- Organism-aware optimization logic

##  Acknowledgements
Special thanks to GeeksforGeeks, brocode on youtube and my own late-night debugging for helping me understand how to build and think through this project.
