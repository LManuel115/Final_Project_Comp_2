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
# We can use a verb placeholder to show what our string of letters is
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



## Part 2
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



## Part 3
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



## Part 4
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



# Sequence Annotations
## Part 1

```python
from Bio.SeqRecord import SeqRecord
```


```python
# We can create a seq record by creating a seq object
from Bio.Seq import Seq
```


```python
simple_seq = Seq("GATC")
```


```python
# We are creating the variable we need to create a seq record
simple_seq_r = SeqRecord(simple_seq)
```


```python
# We have created a record with no information in it
simple_seq_r
```




    SeqRecord(seq=Seq('GATC'), id='<unknown id>', name='<unknown name>', description='<unknown description>', dbxrefs=[])




```python
# We can pass an ID to our record
simple_seq_r.id = "AC12345"
```


```python
# We can also pass a description to the record
simple_seq_r.description = "Made up sequence for the VDB Computational Biology Class"
```


```python
# We can print just the description of the record
print(simple_seq_r.description)
```

    Made up sequence for the VDB Computational Biology Class



```python
# We can print out just our recorded sequence
simple_seq_r.seq
```




    Seq('GATC')




```python
# We can print out everything attached to the record
simple_seq_r
```




    SeqRecord(seq=Seq('GATC'), id='AC12345', name='<unknown name>', description='Made up sequence for the VDB Computational Biology Class', dbxrefs=[])




```python
# We can write our seq record to a certain file
simple_seq_r.annotations["evidence"] = "None. This is just an example"
```


```python
print(simple_seq_r.annotations["evidence"])
```

    None. This is just an example



```python
simple_seq_r
```




    SeqRecord(seq=Seq('GATC'), id='AC12345', name='<unknown name>', description='Made up sequence for the VDB Computational Biology Class', dbxrefs=[])




```python
# We can do per letter annotations
simple_seq_r.letter_annotations["phred_quality"] = [40, 40, 38, 30]
```


```python
print(simple_seq_r.letter_annotations)
```

    {'phred_quality': [40, 40, 38, 30]}



```python
# We can import a fasta file
from Bio import SeqIO
```


```python
# We can save the fasta file as a record
record = SeqIO.read("NC_005816.fna.txt","fasta")
```


```python
record
```




    SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='gi|45478711|ref|NC_005816.1|', name='gi|45478711|ref|NC_005816.1|', description='gi|45478711|ref|NC_005816.1| Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=[])




```python
# We can look at the sequence individually
record.seq
```




    Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG')




```python
# We can look at the ID individually
record.id
```




    'gi|45478711|ref|NC_005816.1|'




```python
# We can look at the description individually
record.description
```




    'gi|45478711|ref|NC_005816.1| Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence'




```python
# We can only access information that is actually within the record
record.dbxrefs
```




    []




```python
record.annotations
```




    {}




```python
record.features
```




    []



## Part 2
```python
# We can upload a Gen Bank file
record = SeqIO.read("NC_005816.gb.txt","genbank")
```


```python
# We can access the information in the Gen Bank file
record
```




    SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=['Project:58037'])




```python
record.seq
```




    Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG')




```python
record.id
```




    'NC_005816.1'




```python
record.name
```




    'NC_005816'




```python
record.description
```




    'Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence'




```python
record.letter_annotations
```




    {}




```python
# We can look at the length of the record's annotations
len(record.annotations)
```




    13




```python
# We can access the annotations of the file
record.annotations["source"]
```




    'Yersinia pestis biovar Microtus str. 91001'




```python
record.dbxrefs
```




    ['Project:58037']




```python
# We can put all of the records in a features table
len(record.features)
```




    41




```python
# We can add concrete information on the genes
from Bio import SeqFeature
```


```python
# We can say around where the start position starts
start_pos = SeqFeature.AfterPosition(5)
```


```python
# We can say the end position is in a certain area
end_pos = SeqFeature.BetweenPosition(9, left = 8, right = 9)
```


```python
# We can create a location based on the info above
my_location = SeqFeature.SimpleLocation(start_pos, end_pos)
```


```python
print(my_location)
```

    [>5:(8^9)]



```python
# We can pull out just the start location
my_location.start
```




    AfterPosition(5)




```python
# We can pull out just the end location
my_location.end
```




    BetweenPosition(9, left=8, right=9)




```python
# We can make it just print the number of the end or start position
int(my_location.end)
```




    9




```python
int(my_location.start)
```




    5




```python
# We can just pass the numbers to the verbs to print just the numbers
exact_location = SeqFeature.SimpleLocation(5,9)
```


```python
print(exact_location)
```

    [5:9]



```python
exact_location.start
```




    ExactPosition(5)



## Part 3
```python
from Bio.SeqRecord import SeqRecord
```


```python
# We can create a SeqRecord
record = SeqRecord(Seq("MMYQQGCFAGGTVLRLAKDLAENNRGARVLVVCSEITAVTFRGPSETHDSMVGQALFGDGAGAVIVGSDPDLSVERPLYELVWTGATLLPOSEGAIDGHLREVGLTFHLLKDVPGLISKNIEKSLKEAFTPLGISDWVSTFWIAHPGGPAILDQVEAKLGLKEEKMRATREVLSEYGNMSSAC"), id= "gi|14150838|gb|AAK54648.1|AF376133_1", description= "chalcone synthase [Cucumis sativus]",)
```


```python
# We can print our SeqRecord as a fasta file
print(record.format("fasta"))
```

    >gi|14150838|gb|AAK54648.1|AF376133_1 chalcone synthase [Cucumis sativus]
    MMYQQGCFAGGTVLRLAKDLAENNRGARVLVVCSEITAVTFRGPSETHDSMVGQALFGDG
    AGAVIVGSDPDLSVERPLYELVWTGATLLPOSEGAIDGHLREVGLTFHLLKDVPGLISKN
    IEKSLKEAFTPLGISDWVSTFWIAHPGGPAILDQVEAKLGLKEEKMRATREVLSEYGNMS
    SAC
    



```python
# We can print the SeqRecord not as a fasta file
print(record)
```

    ID: gi|14150838|gb|AAK54648.1|AF376133_1
    Name: <unknown name>
    Description: chalcone synthase [Cucumis sativus]
    Number of features: 0
    Seq('MMYQQGCFAGGTVLRLAKDLAENNRGARVLVVCSEITAVTFRGPSETHDSMVGQ...SAC')



```python
from Bio import SeqIO
```


```python
# We can import our GenBank file
record = SeqIO.read("NC_005816.gb.txt", "genbank")
```


```python
record
```




    SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=['Project:58037'])




```python
# We can print the length of the record
len(record)
```




    9609




```python
# We can determine how many features are in the record
len(record.features)
```




    41




```python
# We can determine a specific feature
print(record.features[20])
```

    type: gene
    location: [4342:4780](+)
    qualifiers:
        Key: db_xref, Value: ['GeneID:2767712']
        Key: gene, Value: ['pim']
        Key: locus_tag, Value: ['YP_pPCP05']
    



```python
print(record.features[21])
```

    type: CDS
    location: [4342:4780](+)
    qualifiers:
        Key: codon_start, Value: ['1']
        Key: db_xref, Value: ['GI:45478716', 'GeneID:2767712']
        Key: gene, Value: ['pim']
        Key: locus_tag, Value: ['YP_pPCP05']
        Key: note, Value: ['similar to many previously sequenced pesticin immunity protein entries of Yersinia pestis plasmid pPCP, e.g. gi| 16082683|,ref|NP_395230.1| (NC_003132) , gi|1200166|emb|CAA90861.1| (Z54145 ) , gi|1488655| emb|CAA63439.1| (X92856) , gi|2996219|gb|AAC62543.1| (AF053945) , and gi|5763814|emb|CAB531 67.1| (AL109969)']
        Key: product, Value: ['pesticin immunity protein']
        Key: protein_id, Value: ['NP_995571.1']
        Key: transl_table, Value: ['11']
        Key: translation, Value: ['MGGGMISKLFCLALIFLSSSGLAEKNTYTAKDILQNLELNTFGNSLSHGIYGKQTTFKQTEFTNIKSNTKKHIALINKDNSWMISLKILGIKRDEYTVCFEDFSLIRPPTYVAIHPLLIKKVKSGNFIVVKEIKKSIPGCTVYYH']
    



```python
# We can subdivide the record
sub_record = record[4300:4800]
```


```python
len(sub_record)
```




    500




```python
# We can determine how many features are in the subrecord
len(sub_record.features)
```




    2




```python
# We can look at each of the features individually
sub_record.features[0]
```




    SeqFeature(SimpleLocation(ExactPosition(42), ExactPosition(480), strand=1), type='gene', qualifiers=...)




```python
sub_record.features[1]
```




    SeqFeature(SimpleLocation(ExactPosition(42), ExactPosition(480), strand=1), type='CDS', qualifiers=...)




```python
# We can look at the feature in the context of our subrecord
print(sub_record.features[0])
```

    type: gene
    location: [42:480](+)
    qualifiers:
        Key: db_xref, Value: ['GeneID:2767712']
        Key: gene, Value: ['pim']
        Key: locus_tag, Value: ['YP_pPCP05']
    



```python
print(sub_record.features[1])
```

    type: CDS
    location: [42:480](+)
    qualifiers:
        Key: codon_start, Value: ['1']
        Key: db_xref, Value: ['GI:45478716', 'GeneID:2767712']
        Key: gene, Value: ['pim']
        Key: locus_tag, Value: ['YP_pPCP05']
        Key: note, Value: ['similar to many previously sequenced pesticin immunity protein entries of Yersinia pestis plasmid pPCP, e.g. gi| 16082683|,ref|NP_395230.1| (NC_003132) , gi|1200166|emb|CAA90861.1| (Z54145 ) , gi|1488655| emb|CAA63439.1| (X92856) , gi|2996219|gb|AAC62543.1| (AF053945) , and gi|5763814|emb|CAB531 67.1| (AL109969)']
        Key: product, Value: ['pesticin immunity protein']
        Key: protein_id, Value: ['NP_995571.1']
        Key: transl_table, Value: ['11']
        Key: translation, Value: ['MGGGMISKLFCLALIFLSSSGLAEKNTYTAKDILQNLELNTFGNSLSHGIYGKQTTFKQTEFTNIKSNTKKHIALINKDNSWMISLKILGIKRDEYTVCFEDFSLIRPPTYVAIHPLLIKKVKSGNFIVVKEIKKSIPGCTVYYH']
    



```python
# We can look at the annotations of the subrecord
sub_record.annotations
```




    {'molecule_type': 'DNA'}




```python
sub_record.dbxrefs
```




    []




```python
# We can add to the annotations of the subrecord
sub_record.annotations["topology"] = "linear"
```


```python
sub_record.annotations
```




    {'molecule_type': 'DNA', 'topology': 'linear'}




```python
# We can look at the id of the subrecord
sub_record.id
```




    'NC_005816.1'




```python
# We can look at the name of the subrecord
sub_record.name
```




    'NC_005816'




```python
# We can look at the description of the subrecord
sub_record.description
```




    'Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence'




```python
# We can change the description of the subrecord
sub_record.description = "Versinia pestis biovar Microtus str. 91001 plasmid pPCP1, partial sequence"
```


```python
sub_record.description
```




    'Versinia pestis biovar Microtus str. 91001 plasmid pPCP1, partial sequence'




```python
# We can format the subrecord format
print(sub_record.format("genbank")[:200] + "...")
```

    LOCUS       NC_005816                500 bp    DNA     linear   UNK 01-JAN-1980
    DEFINITION  Versinia pestis biovar Microtus str. 91001 plasmid pPCP1, partial
                sequence.
    ACCESSION   NC_00581...


## Part 4
```python
# We are going to go back to using our downloaded Gen bank file
record = SeqIO.read("NC_005816.gb.txt", "genbank")
```


```python
record
```




    SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=['Project:58037'])




```python
# We can do the length of the record
len(record)
```




    9609




```python
len(record.features)
```




    41




```python
# We can see the project number of the record
record.dbxrefs
```




    ['Project:58037']




```python
# We can look at the annotation keys for the record
record.annotations.keys()
```




    dict_keys(['molecule_type', 'topology', 'data_file_division', 'date', 'accessions', 'sequence_version', 'gi', 'keywords', 'source', 'organism', 'taxonomy', 'references', 'comment'])




```python
# We can shift the start point of the record
shifted = record[2000:] + record[:2000]
```


```python
shifted
```




    SeqRecord(seq=Seq('GATACGCAGTCATATTTTTTACACAATTCTCTAATCCCGACAAGGTCGTAGGTC...GGA'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=[])




```python
# We can look at the length of the shifted genome
len(shifted)
```




    9609




```python
# We can look at the features of the shifted genome
len(shifted.features)
```




    40




```python
# We can look at the annotation keys of the shifted genome
shifted.annotations.keys()
```




    dict_keys(['molecule_type'])




```python
# We can look at the reference for the shifted genome
shifted.dbxrefs
```




    []




```python
# We have to tell the program to keep the annotations keys and everything else the same on the shifted genome
shifted.dbxrefs = record.dbxrefs[:]
```


```python
shifted.dbxrefs
```




    ['Project:58037']




```python
shifted.annotations = record.annotations.copy()
```


```python
shifted.annotations.keys()
```




    dict_keys(['molecule_type', 'topology', 'data_file_division', 'date', 'accessions', 'sequence_version', 'gi', 'keywords', 'source', 'organism', 'taxonomy', 'references', 'comment'])




```python
record
```




    SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=['Project:58037'])




```python
# We can tell the record to print our genome as both a string value and integers for certain values and to show the properties it has
print("%s %i %i %i %i" % (record.id, len(record), len(record.features), len(record.dbxrefs), len(record.annotations)))
```

    NC_005816.1 9609 41 1 13



```python
# We can run the reverse complement of the record sequence
rc = record.reverse_complement(id = "Testing")
```


```python
# We can test the properties of the reverse complement
print("%s %i %i %i %i" % (rc.id, len(rc), len(rc.features), len(rc.dbxrefs), len(rc.annotations)))
```

    Testing 9609 41 0 0






