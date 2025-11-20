# **ARP:A Configurable Heuristic for the MLCS Problem**

**Authors:** Farhana Akter Tumpa, Rebin Silva Valan Arasu, Rajiv Gupta

This repository provides the official implementation of **ARP**, an efficient heuristic algorithm for the **Multiple Longest Common Subsequence (MLCS)** problem.  
ARP introduces selection-based search that constructs high-quality MLCS solutions through three key components:

- **Addition-Δs** — incrementally adding subsequences from a chosen primary string  
- **Replacement** — replacing subsequences across candidate solutions to increase diversity  
- **Prioritization** — ranking character groups to guide the Δ-addition order  
ARP provides configurable trade-offs between runtime and solution quality and consistently produces strong MLCS results across multiple benchmark datasets.

---
## Environment
- **Supported Python versions:** 3.9 – 3.13  
- **Tested on:** Python 3.13.5  
---
## **Datasets**

This repository includes three benchmark datasets commonly used in MLCS research:

- **ACO-Random**  
- **ACO-Virus**  
- **ACO-Rat**  

Each dataset contains **N sequences of equal length**, with one sequence per line.  

---
## **How to Run ARP**

ARP supports multiple execution modes. 

---

### **Addition Algorithm A(0)**  

Runs the core Δ-addition heuristic without Replacement or Prioritization.

```bash
python3 Code/Addition.py --input path/to/dataset/file.txt

```
### **Replacement Algorithm ARP(0)**  
Runs the complete ARP(0) as described in the paper.
```bash
python3 Code/ARP.py --input path/to/dataset/file.txt --maxdepth <VALUE>

```
### **ARP Algorithm ARP(2000)**  
This mode runs the full ARP heuristic using a beam-guided search.
This is the strongest configuration used in the paper and typically yields the best MLCS quality.
```bash
python3 Code/reorder.py --input path/to/dataset/file.txt --beam <BEAM_WIDTH>
```
Runs the complete ARP heuristic as described in the paper. The dataset is structured so that Sequence 1 is the best primary (highest MLCS), so the script outputs results only for Primary 1.


### **Primary Sequence Ranking**

Predicts which sequences are most likely to be best primaries and runs ARP only on the top-k.
```bash
python3 Code/primary_sequnced_rank_ARP.py --input path/to/dataset/file.txt
```
This significantly reduces runtime while still producing best solution quality.

