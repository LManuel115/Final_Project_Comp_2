# Final_Project_Comp_2
This is the portfolio of my code for the Final Project of Computational Biology 2 in Winter 2024.


# Sequence Objects
## Part 1
```python
from Bio.Seq import Seq
```


```python
my_seq = Seq("GATCG")
```


```python
for index, letter in enumerate(my_seq):
    print("%i %s" % (index, letter))
```

    0 G
    1 A
    2 T
    3 C
    4 G



```python
# We can also print the length of each sequence
print(len(my_seq))
```

    5



```python
#We can access the different positions
print(my_seq[0])
```

    G



```python
print(my_seq[4])
```

    G



```python
print(my_seq[2])
```

    T



```python
# We can also make a count function to count the occurances of something
Seq("AAAA").count("AA")
```




    2




```python
# We can put in a new sequence
my_seq = Seq("GTCCATTCAGTCCCCAAAAACTTT")
```


```python
# We can find the length of the new sequence
len(my_seq)
```




    24




```python
#We can count the specific letters in the sequence
my_seq.count("G")
```




    2




```python
# We can use this to measure the GC content
100 * (my_seq.count("G") + my_seq.count("C"))/len(my_seq)
```




    41.666666666666664




```python
# We can use Biopython to access common calculations
from Bio.SeqUtils import gc_fraction
```


```python
my_seq = Seq("GTCCATTCAGTCCCCAAAAACTTT")
```


```python
gc_fraction(my_seq)
```




    0.4166666666666667




```python
# We can slice sequences into different parts
my_seq[4:12]
```




    Seq('ATTCAGTC')




```python
# We can do slices with a start and a stop to pick different points in sequence
my_seq[0::3]
```




    Seq('GCTGCAAT')




```python
my_seq[1::3]
```




    Seq('TACTCAAT')




```python
my_seq[2::3]
```




    Seq('CTACCACT')




```python
# We can use negative values to start at the other end of sequence
my_seq[::-1]
```




    Seq('TTTCAAAAACCCCTGACTTACCTG')




```python
# We can turn a sequence back into a string
str(my_seq)
```




    'GTCCATTCAGTCCCCAAAAACTTT'




```python
# We canuse a verb placeholder to show what our string of letters is
fasta_format_string = ">Name/n%s/n" % my_seq
```


```python
print(fasta_format_string)
```

    >Name/nGTCCATTCAGTCCCCAAAAACTTT/n



```python
# We can add two scripts together
seq1 = Seq("ACGT")
seq2 = Seq("AACCGG")
```


```python
seq1 + seq2
```




    Seq('ACGTAACCGG')




```python
# We can add in either order
seq2 + seq1
```




    Seq('AACCGGACGT')




```python
# We can manipulate our string even more
contigs = [Seq("ATG"), Seq("ATCCCG"), Seq("TTTCGA")]
```


```python
# We can make a spacer
spacer = Seq("N" *10)
```


```python
# We can use a join function to take the spacer and join it with the contigs
spacer.join(contigs)
```




    Seq('ATGNNNNNNNNNNATCCCGNNNNNNNNNNTTTCGA')




```python
# We can use this function if we have case sensitive data
dna_seq = Seq("acgtACGT")
```


```python
dna_seq
```




    Seq('acgtACGT')




```python
# We can make the whole sequence either upper or lower case
dna_seq.upper()
```




    Seq('ACGTACGT')




```python
dna_seq.lower()
```




    Seq('acgtacgt')




```python
dna_seq.upper()
```




    Seq('ACGTACGT')




```python
# We can search for certain codons in a sequence
"gtac" in dna_seq
```




    False




```python
#We must save the capitalized sequence as a variable
dna_seq = dna_seq.upper()
```


```python
"GTAC" in dna_seq
```




    True




```python
my_seq = Seq("GATCGAAAATTTCGCGCGAATT")
```


```python
# We can create the compliment of the sequence
my_seq.complement()
```




    Seq('CTAGCTTTTAAAGCGCGCTTAA')




```python
# We can also create the reverse compliment of the sequence
my_seq.reverse_complement()
```




    Seq('AATTCGCGCGAAATTTTCGATC')




```python
# We can also work with protein seqeunces
protein_seq = Seq("EVRNAK")
protein_seq.complement()
```




    Seq('EBYNTM')




```python
# We can create a code of DNA
coding_dna = Seq("ATGGCCATTGTAATCCCCCGCTGAAAGGGTGCCCGATAG")
```


```python
coding_dna
```




    Seq('ATGGCCATTGTAATCCCCCGCTGAAAGGGTGCCCGATAG')




```python
# We're creating a new variable
template_dna = coding_dna.reverse_complement()
```


```python
template_dna
```




    Seq('CTATCGGGCACCCTTTCAGCGGGGGATTACAATGGCCAT')




```python
# we can convert DNA into RNA
messenger_rna = coding_dna.transcribe()
```


```python
messenger_rna
```




    Seq('AUGGCCAUUGUAAUCCCCCGCUGAAAGGGUGCCCGAUAG')




```python
# This is another way to transcribe DNA into RNA
template_dna.reverse_complement().transcribe()
```




    Seq('AUGGCCAUUGUAAUCCCCCGCUGAAAGGGUGCCCGAUAG')




```python
# We transcribe RNA back into DNA
messenger_rna.back_transcribe()
```




    Seq('ATGGCCATTGTAATCCCCCGCTGAAAGGGTGCCCGATAG')




```python
# We can translate RNA into a protein sequence
messenger_rna
```




    Seq('AUGGCCAUUGUAAUCCCCCGCUGAAAGGGUGCCCGAUAG')




```python
messenger_rna.translate()
```




    Seq('MAIVIPR*KGAR*')




```python
# We can specify the codon tables between a nuclear and mitochondrial codon
coding_dna.translate(table="Vertebrate Mitochondrial")
```




    Seq('MAIVIPRWKGAR*')




```python
# we can do this using a table number instead also
coding_dna.translate(table = 2)
```




    Seq('MAIVIPRWKGAR*')




```python
# We can stop at the certain stop codons in a seqeunce
coding_dna.translate(to_stop = True)
```




    Seq('MAIVIPR')




```python
coding_dna.translate(table = 2, to_stop=True)
```




    Seq('MAIVIPRWKGAR')




```python
# We can specify our stop codon via a symbol
coding_dna.translate(table = 2, stop_symbol = "!")
```




    Seq('MAIVIPRWKGAR!')




```python
# We can also handle nonstandard stop codons
gene = Seq("GTGAAAAAGATGCAATCTATCGTACTCGCACTTTCCCTGGTTCTGGTCGCTCCCATGGCAGCACAGGCTGCGGAAATTACGTTAGTCCCGTCAGTAAAATTACAGATAGGCGATCGTGATAATCGTGGCTATTACTGGGATGGAGGTCACTGGCGCGACCACGGCTGGTGGAAACAACATTATGAATGGCGAGGCAATCGCTGGCACCTACACGGACCGCCGCCACCGCCGCGCCACCATAAGAAAGCTCCTCATGATCATCACGGCGGTCATGGTCCAGGCCAAACATCACCGCTATAG")
```


```python
#We can translate this sequence into a bacterial sequence
gene.translate(table = "Bacterial")
```




    Seq('VKKMQSIVLALSLVLVAPMAAQAAEITLVPSVKLQIGDRDNRGYYWDGGHWRDH...PL*')




```python
# we can tell the bacterial sequence where to stop
gene.translate(table = "Bacterial", to_stop = True)
```




    Seq('VKKMQSIVLALSLVLVAPMAAQAAEITLVPSVKLQIGDRDNRGYYWDGGHWRDH...SPL')




```python
# We can tell biopython that a sequence is a complete gene sequence
gene.translate(table = "Bacterial", cds = True)
```




    Seq('MKKMQSIVLALSLVLVAPMAAQAAEITLVPSVKLQIGDRDNRGYYWDGGHWRDH...SPL')




```python
# We can also import a codon table
from Bio.Data import CodonTable
```


```python
standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
```


```python
# We can import table 2 and label it as "Vertebrate Mitochondrial"
mito_table = CodonTable.unambiguous_dna_by_name["Vertebrate Mitochondrial"]
```


```python
# We can print the codon tables
print(standard_table)
```

    Table 1 Standard, SGC0
    
      |  T      |  C      |  A      |  G      |
    --+---------+---------+---------+---------+--
    T | TTT F   | TCT S   | TAT Y   | TGT C   | T
    T | TTC F   | TCC S   | TAC Y   | TGC C   | C
    T | TTA L   | TCA S   | TAA Stop| TGA Stop| A
    T | TTG L(s)| TCG S   | TAG Stop| TGG W   | G
    --+---------+---------+---------+---------+--
    C | CTT L   | CCT P   | CAT H   | CGT R   | T
    C | CTC L   | CCC P   | CAC H   | CGC R   | C
    C | CTA L   | CCA P   | CAA Q   | CGA R   | A
    C | CTG L(s)| CCG P   | CAG Q   | CGG R   | G
    --+---------+---------+---------+---------+--
    A | ATT I   | ACT T   | AAT N   | AGT S   | T
    A | ATC I   | ACC T   | AAC N   | AGC S   | C
    A | ATA I   | ACA T   | AAA K   | AGA R   | A
    A | ATG M(s)| ACG T   | AAG K   | AGG R   | G
    --+---------+---------+---------+---------+--
    G | GTT V   | GCT A   | GAT D   | GGT G   | T
    G | GTC V   | GCC A   | GAC D   | GGC G   | C
    G | GTA V   | GCA A   | GAA E   | GGA G   | A
    G | GTG V   | GCG A   | GAG E   | GGG G   | G
    --+---------+---------+---------+---------+--



```python
print(mito_table)
```

    Table 2 Vertebrate Mitochondrial, SGC1
    
      |  T      |  C      |  A      |  G      |
    --+---------+---------+---------+---------+--
    T | TTT F   | TCT S   | TAT Y   | TGT C   | T
    T | TTC F   | TCC S   | TAC Y   | TGC C   | C
    T | TTA L   | TCA S   | TAA Stop| TGA W   | A
    T | TTG L   | TCG S   | TAG Stop| TGG W   | G
    --+---------+---------+---------+---------+--
    C | CTT L   | CCT P   | CAT H   | CGT R   | T
    C | CTC L   | CCC P   | CAC H   | CGC R   | C
    C | CTA L   | CCA P   | CAA Q   | CGA R   | A
    C | CTG L   | CCG P   | CAG Q   | CGG R   | G
    --+---------+---------+---------+---------+--
    A | ATT I(s)| ACT T   | AAT N   | AGT S   | T
    A | ATC I(s)| ACC T   | AAC N   | AGC S   | C
    A | ATA M(s)| ACA T   | AAA K   | AGA Stop| A
    A | ATG M(s)| ACG T   | AAG K   | AGG Stop| G
    --+---------+---------+---------+---------+--
    G | GTT V   | GCT A   | GAT D   | GGT G   | T
    G | GTC V   | GCC A   | GAC D   | GGC G   | C
    G | GTA V   | GCA A   | GAA E   | GGA G   | A
    G | GTG V(s)| GCG A   | GAG E   | GGG G   | G
    --+---------+---------+---------+---------+--



```python
# We can ask to show all of the stop codons for a table
mito_table.stop_codons
```




    ['TAA', 'TAG', 'AGA', 'AGG']




```python
# We can also get it to show all of the start codons
mito_table.start_codons
```




    ['ATT', 'ATC', 'ATA', 'ATG', 'GTG']




```python
#We can compare sequences to see what sequences equal
seq1 = Seq("ACGT")
```


```python
"ACGT" == seq1
```




    True




```python
seq1 == "ACGT"
```




    True




```python
# We can create a sequence with only its length known by using the word "None"
unknown_seq = Seq(None, 10)
```


```python
unknown_seq
```




    Seq(None, length=10)




```python
len(unknown_seq)
```




    10




```python
# We can pull certain sequences from a data table and work with a partially defined sequence
seq = Seq({117512683: "TTGAAAACCTGAATGTGAGAGTCAGTCAAGGATAGT"}, length = 159345973)
```


```python
seq [1000:1020]
```




    Seq(None, length=20)




```python
# We can print out a certain section of a sequence if we have it defined
seq[117512690:117512700]
```




    Seq('CCTGAATGTG')




```python
# We can show the sequence from a start until the end
seq[117512670:]
```




    Seq({13: 'TTGAAAACCTGAATGTGAGAGTCAGTCAAGGATAGT'}, length=41833303)




```python
seq= Seq("ACGT")
```


```python
undefined_seq = Seq(None, length = 10)
```


```python
# We can add our undefined sequence into the middle of our defined sequence
seq + undefined_seq + seq
```




    Seq({0: 'ACGT', 14: 'ACGT'}, length=18)




```python
my_seq = Seq("GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA")
```


```python
# We can create a mutable sequence by importing one from Biopython
from Bio.Seq import MutableSeq
```


```python
mutable_seq = MutableSeq(my_seq)
```


```python
mutable_seq
```




    MutableSeq('GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA')




```python
# We can change a specific nucleotide in a mutable sequence
mutable_seq[5] = "C"
```


```python
mutable_seq
```




    MutableSeq('GCCATCGTAATGGGCCGCTGAAAGGGTGCCCGA')




```python
#We can remove a certain nucleotide from a sequence
mutable_seq.remove("T")
```


```python
mutable_seq
```




    MutableSeq('GCCACGTAATGGGCCGCTGAAAGGGTGCCCGA')




```python
# We can reverse our sequence
mutable_seq.reverse()
```


```python
mutable_seq
```




    MutableSeq('AGCCCGTGGGAAAGTCGCCGGGTAATGCACCG')




```python
#We can reprotect our code from mutations
new_seq = Seq(mutable_seq)
```


```python
new_seq
```




    Seq('AGCCCGTGGGAAAGTCGCCGGGTAATGCACCG')




```python
# We can import commmands to help not use the sequence objects from Biopython
from Bio.Seq import reverse_complement, transcribe, back_transcribe, translate
```


```python
my_string = "GTCGATTTGCGAGTTGCTGA"
```


```python
reverse_complement(my_string)
```




    'TCAGCAACTCGCAAATCGAC'




```python
transcribe(my_string)
```




    'GUCGAUUUGCGAGUUGCUGA'




```python
back_transcribe(my_string)
```




    'GTCGATTTGCGAGTTGCTGA'




```python
translate(my_string)
```




    'VDLRVA'




