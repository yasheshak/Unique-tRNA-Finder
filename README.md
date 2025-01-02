## **Overview**
This Python program identifies unique subsequences in mitochondrial tRNA (mt.tRNA) sequences. These unique tags help assess tRNA molecule abundance. The program:
- Reads FASTA-formatted sequences from STDIN.
- Finds unique and minimal subsequences for each tRNA.
- Outputs results sorted by tRNA header with aligned formatting.

---

## **Deliverables**
- `findUnique.py`: Script for identifying unique subsequences.
- Notebook: Includes **Inspection Materials** and **Results**.
- Output: Written to STDOUT, sorted by tRNA header.

---

## **Features**
1. **Input**:
   - Reads FASTA sequences via STDIN.
   - Removes alignment characters (`-`, `_`, `.`).
   - Example: `python findUnique.py < tRNA_sequences.fa`.

2. **Processing**:
   - Uses Python sets to compute unique subsequences.
   - Ensures minimal subsequences (no larger redundant substrings).

3. **Output**:
   - Sorted by tRNA header.
   - Aligned subsequences with dots (`.`) for visual alignment.
