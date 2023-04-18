from Bio import SeqIO
from collections import defaultdict
import math

class CodonUsageAnalysis:
    """
    Modified : 2023-04-18
    Author: Ruixuan

    Input:
    - gene_list : is the list of genes in selection, such as ribosomal proteins, etc.
    - cds_path : is the path to the cds file of the corresponding organisms
    - codon_table : is the codon table where the key is amino acid and the value is a list of correponding codon. A standard genetic code has been implemented

    Dependencies:

    - from Bio import SeqIO: This function needs Biopython for reading fasta file
    - from Bio.Data import CodonTable: This is the codon table
    - from collections import defaultdict: for dictionaries
    """

    standard_genetic_code = {
        'Ala': ['GCT', 'GCC', 'GCA', 'GCG'],
        'Arg': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
        'Asn': ['AAT', 'AAC'],
        'Asp': ['GAT', 'GAC'],
        'Cys': ['TGT', 'TGC'],
        'Gln': ['CAA', 'CAG'],
        'Glu': ['GAA', 'GAG'],
        'Gly': ['GGT', 'GGC', 'GGA', 'GGG'],
        'His': ['CAT', 'CAC'],
        'Ile': ['ATT', 'ATC', 'ATA'],
        'Leu': ['CTT', 'CTC', 'CTA', 'CTG', 'TTA', 'TTG'],
        'Lys': ['AAA', 'AAG'],
        'Met': ['ATG'],
        'Phe': ['TTT', 'TTC'],
        'Pro': ['CCT', 'CCC', 'CCA', 'CCG'],
        'Ser': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
        'Thr': ['ACT', 'ACC', 'ACA', 'ACG'],
        'Trp': ['TGG'],
        'Tyr': ['TAT', 'TAC'],
        'Val': ['GTT', 'GTC', 'GTA', 'GTG'],
        'Stop': ['TAA', 'TAG', 'TGA']
    }

    def __init__(self, reference_list, ref_path, target_list, codon_table=None, target_path=None):
        self.gene_list = reference_list
        self.target_list = target_list
        self.cds_path = ref_path
        self.target_path = target_path if target_path is not None else ref_path
        self.codon_table = codon_table if codon_table is not None else self.standard_genetic_code
        self.record_dict = SeqIO.to_dict(SeqIO.parse(ref_path, "fasta"))
    
    
    def RSCU_dict(self):
        """
        This function aims to calculate the relative synonymous codon usage of the reference gene set
        The logic of this part is 
        1. Load sequences from the reference gene list and the corresponding cds_path
        2. Condition 1: Check if the gene can be exactly divided by 3. If not, save in error_seq_list. If true, continue.
        3. Define dictionary codon_counter, in which the key is codon and the value is its frequency in all reference genes
        4. Define dicionary AA_counter, in which the key is the AA and the value is the total frequency of all synonymous codons
        5. Define RSCU_dict, in which the key is the codon and the value is the RSCU
        6. Return RSCU_dict and error sequence list
        """

        codon_counter = defaultdict(int) # the key is the codon and the value is the frequency
        AA_counter = defaultdict(int)
        self.RSCU_dict = defaultdict(float)
        self.error_ref_list = []

        # Load genes and their sequences and then count the codon if the sequences can be exactly divided by 3 

        for gene in self.gene_list:
            sequences = str(self.record_dict[gene].seq)
            if len(sequences) % 3 != 0:
                print(f"The length of {gene} cannot be exactly divided by 3")
                self.error_ref_list.append(gene)
                continue
            else:
                for i in range(0, len(sequences), 3):
                    codon_code = sequences[i: i+3]
                    codon_counter[codon_code] +=1

        for AA in self.codon_table.keys():
            total_AA = 0
            for codon in self.codon_table[AA]:
                freq = codon_counter[codon]
                total_AA += freq
            AA_counter[AA] = total_AA

        for AA in self.codon_table.keys():
            if AA_counter[AA] != 0:
                for codon in self.codon_table[AA]:
                    freq = codon_counter[codon]
                    ratio = freq / AA_counter[AA]
                    self.RSCU_dict[codon] = ratio
            else:
                num_codon = len(self.codon_table[AA])
                for codon in self.codon_table[AA]:
                    ratio = 1 / num_codon
                    self.RSCU_dict[codon] = ratio

        print(f"There are {len(self.error_ref_list)} has wrong length in the reference gene set")
        print("The error gene list is:", self.error_ref_list)
        return self.RSCU_dict
    
    def get_error_ref_list(self):
        return self.error_ref_list

    
    def RSCU_max_dict(self):
        self.RSCU_max_dict = defaultdict(float)
        for AA in self.codon_table.keys():
            AA_RSCU_list = []
            for codon in self.codon_table[AA]:
                AA_RSCU_list.append(self.RSCU_dict[codon])
            self.RSCU_max_dict[AA] = max(AA_RSCU_list)

        return self.RSCU_max_dict
    
    def w_dict(self):
        self.w_dict = defaultdict(float)

        for AA in self.codon_table.keys():
            for codon in self.codon_table[AA]:
                self.w_dict[codon] = self.RSCU_dict[codon] / self.RSCU_max_dict[AA]

        return self.w_dict

    # This part aims to calculate the CAI value of given gene sets
    # 1. This function will load genes in the target_list and read a separated file of cds fasta if provided
    # 2. This function check if the cds region can be divided by 3,
    #   If false, error gene id will be saved in error_gene_list
    #   If true, genes will be separated every 3 nucleotide for check the corresponding codon, count the frequency and calculate CAI based on inherited RSCU_dict
    # 3. Finally, when all genes have been treated, output two following things
    #   3.1 Print number of genes which had a wrong length and the error list
    #   3.2 Output the CAI dictionary where the key is the gene and the value is the CAI value

    def CAI_calculate(self):
        self.CAI_dict = defaultdict(float)
        self.error_target_list = []

        target_record_dict = SeqIO.to_dict(SeqIO.parse(self.target_path, "fasta"))
        
        for gene in self.target_list:
            sequences = str(target_record_dict[gene].seq)
            product_w = 1
            n_codons = 0
            
            if len(sequences) % 3 != 0:
                print(f"The length of {gene} cannot be exactly divided by 3")
                self.error_target_list.append(gene)
                continue
            else:
                codon_counter = defaultdict(int)
                for i in range(0, len(sequences), 3):
                    codon_code = sequences[i: i+3]
                    codon_counter[codon_code] +=1
            
                for codon in self.w_dict.keys():
                    if codon_counter[codon] > 0:
                        product_w *= self.w_dict[codon] ** codon_counter[codon]
                        n_codons += codon_counter[codon]

                CAI = product_w ** (1 / n_codons)
                self.CAI_dict[gene] = CAI
        
        print(f"There are {len(self.error_target_list)} has wrong length in the target gene set")
        print("The error gene list is:", self.error_target_list)
        
    
    def get_error_seq_list(self):
        return self.error_target_list
    
    def get_CAI_dict(self):
        return self.CAI_dict
    
# What is the meaning of class method?

    @classmethod
    def all(cls, ref_list, ref_path, target_list, codon_table=None, target_path=None):
        analysis = cls(ref_list, ref_path, target_list, codon_table=codon_table, target_path=target_path)
        analysis.RSCU_dict()
        analysis.RSCU_max_dict()
        analysis.w_dict()
        analysis.CAI_calculate()
        CAI_dict = analysis.get_CAI_dict()
        return CAI_dict