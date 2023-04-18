# Codon Usage Analysis

## Context

- There are many web servers allowing for calculating CAI score or relative synonymous codon usage (RSCU). However, most of them are designed for model organisms such as E. coli. 
- Web servers usually limit the size of your uploading files
- Biopython is hard to use for calculating CAI, and also limited to specific organisms. 

## Main function of this simple script

### Sequence feature extraction

- Calculate `RSCU` and relative adaptiveness (`w`) of codons based on a set of coding sequences.
- Calculate `CAI` of any given coding sequences using `RSCU`
- Calculate `GC%` and `GC3` of any given sequences

### Current usage [2023-04-18]

**What needs to be prepared?**

- A list of genes as reference. This list will be used to calculate RSCU 
- A list of genes as target. This list is the genes you want to calculate CAI
- A path of the fasta file of the coding sequences of reference genes (`cds.fna`)
- (option) a path of the fasta of the coding sequences of target genes, by default, it is the same with reference one.
- (option) genetic code where the key is the amino acid and the value is a list of its corresponding codons. The standard genetic code is implemented by default

**One line operation**

```python
CAI_dict = CodonUsageAnalysis.all(refer_list, ref_path, target_list, target_path=target_path)
```

This operation gives you a dict variable where the key is the gene in your target list and the value is the CAI. The class `all` is made up of following methods: 

- `RSCU_dict()` : calculate the RSCU values for each codon based on the reference genes `{codon : RSCU}`
- `RSCU_max_dict()` : get the maximum RSCU among synomyous codons `{AA : RSCU_max}`
- `w_dict()` : get the relative adaptiveness for each codon `{codon : w}`
- `CAI_calculate()` : get the CAI for each gene 
- `get_CAI_dict()` : get the CAI dictionary `{gene : CAI}
- `get_error_ref_list()` and `get_error_target_list()` : these two methods return the genes that can not be divided by 3

```python
analysis = CodonUsageAnalysis(refer_list, ref_path, target_list, target_path=target_path)
RSCU_dict = analysis.RSCU_dict()
RSCU_max_dict = analysis.RSCU_max_dict()
w_dict = analysis.w_dict()
analysis.CAI_calculate()
CAI_dict = analysis.get_CAI_dict()
error_ref_list = analysis.get_error_ref_list()
error_target_list = analysis.get_error_target_list()
```

