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




# Sequence I/O
## Part 1
```python
from Bio import SeqIO
```


```python
# We can use to parse function on our file and print certain parts of the parsed file
for seq_record in SeqIO.parse("ls_orchid.fasta.txt", "fasta"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))
```

    gi|2765658|emb|Z78533.1|CIZ78533
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGG...CGC')
    740
    gi|2765657|emb|Z78532.1|CCZ78532
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAACAG...GGC')
    753
    gi|2765656|emb|Z78531.1|CFZ78531
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGCAG...TAA')
    748
    gi|2765655|emb|Z78530.1|CMZ78530
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAAACAACAT...CAT')
    744
    gi|2765654|emb|Z78529.1|CLZ78529
    Seq('ACGGCGAGCTGCCGAAGGACATTGTTGAGACAGCAGAATATACGATTGAGTGAA...AAA')
    733
    gi|2765652|emb|Z78527.1|CYZ78527
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAG...CCC')
    718
    gi|2765651|emb|Z78526.1|CGZ78526
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAG...TGT')
    730
    gi|2765650|emb|Z78525.1|CAZ78525
    Seq('TGTTGAGATAGCAGAATATACATCGAGTGAATCCGGAGGACCTGTGGTTATTCG...GCA')
    704
    gi|2765649|emb|Z78524.1|CFZ78524
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATAGTAG...AGC')
    740
    gi|2765648|emb|Z78523.1|CHZ78523
    Seq('CGTAACCAGGTTTCCGTAGGTGAACCTGCGGCAGGATCATTGTTGAGACAGCAG...AAG')
    709
    gi|2765647|emb|Z78522.1|CMZ78522
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGCAG...GAG')
    700
    gi|2765646|emb|Z78521.1|CCZ78521
    Seq('GTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAGAATATATGATCGAGT...ACC')
    726
    gi|2765645|emb|Z78520.1|CSZ78520
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGCAG...TTT')
    753
    gi|2765644|emb|Z78519.1|CPZ78519
    Seq('ATATGATCGAGTGAATCTGGTGGACTTGTGGTTACTCAGCTCGCCATAGGCTTT...TTA')
    699
    gi|2765643|emb|Z78518.1|CRZ78518
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGGAGGATCATTGTTGAGATAGTAG...TCC')
    658
    gi|2765642|emb|Z78517.1|CFZ78517
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAG...AGC')
    752
    gi|2765641|emb|Z78516.1|CPZ78516
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAT...TAA')
    726
    gi|2765640|emb|Z78515.1|MXZ78515
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGCTGAGACCGTAG...AGC')
    765
    gi|2765639|emb|Z78514.1|PSZ78514
    Seq('CGTAACAAGGTTTCCGTAGGTGGACCTTCGGGAGGATCATTTTTGAAGCCCCCA...CTA')
    755
    gi|2765638|emb|Z78513.1|PBZ78513
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...GAG')
    742
    gi|2765637|emb|Z78512.1|PWZ78512
    Seq('CGTAACAAGGTTTCCGTAGGTGGACCTTCGGGAGGATCATTTTTGAAGCCCCCA...AGC')
    762
    gi|2765636|emb|Z78511.1|PEZ78511
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTTCGGAAGGATCATTGTTGAGACCCCCA...GGA')
    745
    gi|2765635|emb|Z78510.1|PCZ78510
    Seq('CTAACCAGGGTTCCGAGGTGACCTTCGGGAGGATTCCTTTTTAAGCCCCCGAAA...TTA')
    750
    gi|2765634|emb|Z78509.1|PPZ78509
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...GGA')
    731
    gi|2765633|emb|Z78508.1|PLZ78508
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...TGA')
    741
    gi|2765632|emb|Z78507.1|PLZ78507
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCCCCA...TGA')
    740
    gi|2765631|emb|Z78506.1|PLZ78506
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCAA...TGA')
    727
    gi|2765630|emb|Z78505.1|PSZ78505
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...TTT')
    711
    gi|2765629|emb|Z78504.1|PKZ78504
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTTCGGAAGGATCATTGTTGAGACCGCAA...TAA')
    743
    gi|2765628|emb|Z78503.1|PCZ78503
    Seq('CGTAACCAGGTTTCCGTAGGTGAACCTCCGGAAGGATCCTTGTTGAGACCGCCA...TAA')
    727
    gi|2765627|emb|Z78502.1|PBZ78502
    Seq('CGTAACCAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGACCGCCA...CGC')
    757
    gi|2765626|emb|Z78501.1|PCZ78501
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCAA...AGA')
    770
    gi|2765625|emb|Z78500.1|PWZ78500
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGCTCATTGTTGAGACCGCAA...AAG')
    767
    gi|2765624|emb|Z78499.1|PMZ78499
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAGGGATCATTGTTGAGATCGCAT...ACC')
    759
    gi|2765623|emb|Z78498.1|PMZ78498
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAAGGTCATTGTTGAGATCACAT...AGC')
    750
    gi|2765622|emb|Z78497.1|PDZ78497
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    788
    gi|2765621|emb|Z78496.1|PAZ78496
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...AGC')
    774
    gi|2765620|emb|Z78495.1|PEZ78495
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...GTG')
    789
    gi|2765619|emb|Z78494.1|PNZ78494
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGGTCGCAT...AAG')
    688
    gi|2765618|emb|Z78493.1|PGZ78493
    Seq('CGTAACAAGGATTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...CCC')
    719
    gi|2765617|emb|Z78492.1|PBZ78492
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...ATA')
    743
    gi|2765616|emb|Z78491.1|PCZ78491
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...AGC')
    737
    gi|2765615|emb|Z78490.1|PFZ78490
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    728
    gi|2765614|emb|Z78489.1|PDZ78489
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GGC')
    740
    gi|2765613|emb|Z78488.1|PTZ78488
    Seq('CTGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACGCAATAATTGATCGA...GCT')
    696
    gi|2765612|emb|Z78487.1|PHZ78487
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TAA')
    732
    gi|2765611|emb|Z78486.1|PBZ78486
    Seq('CGTCACGAGGTTTCCGTAGGTGAATCTGCGGGAGGATCATTGTTGAGATCACAT...TGA')
    731
    gi|2765610|emb|Z78485.1|PHZ78485
    Seq('CTGAACCTGGTGTCCGAAGGTGAATCTGCGGATGGATCATTGTTGAGATATCAT...GTA')
    735
    gi|2765609|emb|Z78484.1|PCZ78484
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGGGGAAGGATCATTGTTGAGATCACAT...TTT')
    720
    gi|2765608|emb|Z78483.1|PVZ78483
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GCA')
    740
    gi|2765607|emb|Z78482.1|PEZ78482
    Seq('TCTACTGCAGTGACCGAGATTTGCCATCGAGCCTCCTGGGAGCTTTCTTGCTGG...GCA')
    629
    gi|2765606|emb|Z78481.1|PIZ78481
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    572
    gi|2765605|emb|Z78480.1|PGZ78480
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    587
    gi|2765604|emb|Z78479.1|PPZ78479
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGT')
    700
    gi|2765603|emb|Z78478.1|PVZ78478
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCAGTGTTGAGATCACAT...GGC')
    636
    gi|2765602|emb|Z78477.1|PVZ78477
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGC')
    716
    gi|2765601|emb|Z78476.1|PGZ78476
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...CCC')
    592
    gi|2765600|emb|Z78475.1|PSZ78475
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GGT')
    716
    gi|2765599|emb|Z78474.1|PKZ78474
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACGT...CTT')
    733
    gi|2765598|emb|Z78473.1|PSZ78473
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGG')
    626
    gi|2765597|emb|Z78472.1|PLZ78472
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    737
    gi|2765596|emb|Z78471.1|PDZ78471
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    740
    gi|2765595|emb|Z78470.1|PPZ78470
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GTT')
    574
    gi|2765594|emb|Z78469.1|PHZ78469
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GTT')
    594
    gi|2765593|emb|Z78468.1|PAZ78468
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...GTT')
    610
    gi|2765592|emb|Z78467.1|PSZ78467
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    730
    gi|2765591|emb|Z78466.1|PPZ78466
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...CCC')
    641
    gi|2765590|emb|Z78465.1|PRZ78465
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGC')
    702
    gi|2765589|emb|Z78464.1|PGZ78464
    Seq('CGTAACAAGGTTTCCGTAGGTGAGCGGAAGGGTCATTGTTGAGATCACATAATA...AGC')
    733
    gi|2765588|emb|Z78463.1|PGZ78463
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGTTCATTGTTGAGATCACAT...AGC')
    738
    gi|2765587|emb|Z78462.1|PSZ78462
    Seq('CGTCACGAGGTCTCCGGATGTGACCCTGCGGAAGGATCATTGTTGAGATCACAT...CAT')
    736
    gi|2765586|emb|Z78461.1|PWZ78461
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...TAA')
    732
    gi|2765585|emb|Z78460.1|PCZ78460
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...TTA')
    745
    gi|2765584|emb|Z78459.1|PDZ78459
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TTT')
    744
    gi|2765583|emb|Z78458.1|PHZ78458
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TTG')
    738
    gi|2765582|emb|Z78457.1|PCZ78457
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...GAG')
    739
    gi|2765581|emb|Z78456.1|PTZ78456
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    740
    gi|2765580|emb|Z78455.1|PJZ78455
    Seq('CGTAACCAGGTTTCCGTAGGTGGACCTTCGGGAGGATCATTTTTGAGATCACAT...GCA')
    745
    gi|2765579|emb|Z78454.1|PFZ78454
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AAC')
    695
    gi|2765578|emb|Z78453.1|PSZ78453
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GCA')
    745
    gi|2765577|emb|Z78452.1|PBZ78452
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GCA')
    743
    gi|2765576|emb|Z78451.1|PHZ78451
    Seq('CGTAACAAGGTTTCCGTAGGTGTACCTCCGGAAGGATCATTGTTGAGATCACAT...AGC')
    730
    gi|2765575|emb|Z78450.1|PPZ78450
    Seq('GGAAGGATCATTGCTGATATCACATAATAATTGATCGAGTTAAGCTGGAGGATC...GAG')
    706
    gi|2765574|emb|Z78449.1|PMZ78449
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGC')
    744
    gi|2765573|emb|Z78448.1|PAZ78448
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGG')
    742
    gi|2765572|emb|Z78447.1|PVZ78447
    Seq('CGTAACAAGGATTCCGTAGGTGAACCTGCGGGAGGATCATTGTTGAGATCACAT...AGC')
    694
    gi|2765571|emb|Z78446.1|PAZ78446
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...CCC')
    712
    gi|2765570|emb|Z78445.1|PUZ78445
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGT')
    715
    gi|2765569|emb|Z78444.1|PAZ78444
    Seq('CGTAACAAGGTTTCCGTAGGGTGAACTGCGGAAGGATCATTGTTGAGATCACAT...ATT')
    688
    gi|2765568|emb|Z78443.1|PLZ78443
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGG')
    784
    gi|2765567|emb|Z78442.1|PBZ78442
    Seq('GTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACATAATAATTGATCGAGT...AGT')
    721
    gi|2765566|emb|Z78441.1|PSZ78441
    Seq('GGAAGGTCATTGCCGATATCACATAATAATTGATCGAGTTAATCTGGAGGATCT...GAG')
    703
    gi|2765565|emb|Z78440.1|PPZ78440
    Seq('CGTAACAAGGTTTCCGTAGGTGGACCTCCGGGAGGATCATTGTTGAGATCACAT...GCA')
    744
    gi|2765564|emb|Z78439.1|PBZ78439
    Seq('CATTGTTGAGATCACATAATAATTGATCGAGTTAATCTGGAGGATCTGTTTACT...GCC')
    592



```python
# We can also use the parse function to split up our gbk file
for seq_record in SeqIO.parse("ls_orchid.gbk.txt", "genbank"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))
```

    Z78533.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGG...CGC')
    740
    Z78532.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAACAG...GGC')
    753
    Z78531.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGCAG...TAA')
    748
    Z78530.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAAACAACAT...CAT')
    744
    Z78529.1
    Seq('ACGGCGAGCTGCCGAAGGACATTGTTGAGACAGCAGAATATACGATTGAGTGAA...AAA')
    733
    Z78527.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAG...CCC')
    718
    Z78526.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAG...TGT')
    730
    Z78525.1
    Seq('TGTTGAGATAGCAGAATATACATCGAGTGAATCCGGAGGACCTGTGGTTATTCG...GCA')
    704
    Z78524.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATAGTAG...AGC')
    740
    Z78523.1
    Seq('CGTAACCAGGTTTCCGTAGGTGAACCTGCGGCAGGATCATTGTTGAGACAGCAG...AAG')
    709
    Z78522.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGCAG...GAG')
    700
    Z78521.1
    Seq('GTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAGAATATATGATCGAGT...ACC')
    726
    Z78520.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGCAG...TTT')
    753
    Z78519.1
    Seq('ATATGATCGAGTGAATCTGGTGGACTTGTGGTTACTCAGCTCGCCATAGGCTTT...TTA')
    699
    Z78518.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGGAGGATCATTGTTGAGATAGTAG...TCC')
    658
    Z78517.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAG...AGC')
    752
    Z78516.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAT...TAA')
    726
    Z78515.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGCTGAGACCGTAG...AGC')
    765
    Z78514.1
    Seq('CGTAACAAGGTTTCCGTAGGTGGACCTTCGGGAGGATCATTTTTGAAGCCCCCA...CTA')
    755
    Z78513.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...GAG')
    742
    Z78512.1
    Seq('CGTAACAAGGTTTCCGTAGGTGGACCTTCGGGAGGATCATTTTTGAAGCCCCCA...AGC')
    762
    Z78511.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTTCGGAAGGATCATTGTTGAGACCCCCA...GGA')
    745
    Z78510.1
    Seq('CTAACCAGGGTTCCGAGGTGACCTTCGGGAGGATTCCTTTTTAAGCCCCCGAAA...TTA')
    750
    Z78509.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...GGA')
    731
    Z78508.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...TGA')
    741
    Z78507.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCCCCA...TGA')
    740
    Z78506.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCAA...TGA')
    727
    Z78505.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...TTT')
    711
    Z78504.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTTCGGAAGGATCATTGTTGAGACCGCAA...TAA')
    743
    Z78503.1
    Seq('CGTAACCAGGTTTCCGTAGGTGAACCTCCGGAAGGATCCTTGTTGAGACCGCCA...TAA')
    727
    Z78502.1
    Seq('CGTAACCAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGACCGCCA...CGC')
    757
    Z78501.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCAA...AGA')
    770
    Z78500.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGCTCATTGTTGAGACCGCAA...AAG')
    767
    Z78499.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAGGGATCATTGTTGAGATCGCAT...ACC')
    759
    Z78498.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAAGGTCATTGTTGAGATCACAT...AGC')
    750
    Z78497.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    788
    Z78496.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...AGC')
    774
    Z78495.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...GTG')
    789
    Z78494.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGGTCGCAT...AAG')
    688
    Z78493.1
    Seq('CGTAACAAGGATTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...CCC')
    719
    Z78492.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...ATA')
    743
    Z78491.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...AGC')
    737
    Z78490.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    728
    Z78489.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GGC')
    740
    Z78488.1
    Seq('CTGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACGCAATAATTGATCGA...GCT')
    696
    Z78487.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TAA')
    732
    Z78486.1
    Seq('CGTCACGAGGTTTCCGTAGGTGAATCTGCGGGAGGATCATTGTTGAGATCACAT...TGA')
    731
    Z78485.1
    Seq('CTGAACCTGGTGTCCGAAGGTGAATCTGCGGATGGATCATTGTTGAGATATCAT...GTA')
    735
    Z78484.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGGGGAAGGATCATTGTTGAGATCACAT...TTT')
    720
    Z78483.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GCA')
    740
    Z78482.1
    Seq('TCTACTGCAGTGACCGAGATTTGCCATCGAGCCTCCTGGGAGCTTTCTTGCTGG...GCA')
    629
    Z78481.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    572
    Z78480.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    587
    Z78479.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGT')
    700
    Z78478.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCAGTGTTGAGATCACAT...GGC')
    636
    Z78477.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGC')
    716
    Z78476.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...CCC')
    592
    Z78475.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GGT')
    716
    Z78474.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACGT...CTT')
    733
    Z78473.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGG')
    626
    Z78472.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    737
    Z78471.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    740
    Z78470.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GTT')
    574
    Z78469.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GTT')
    594
    Z78468.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...GTT')
    610
    Z78467.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    730
    Z78466.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...CCC')
    641
    Z78465.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGC')
    702
    Z78464.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAGCGGAAGGGTCATTGTTGAGATCACATAATA...AGC')
    733
    Z78463.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGTTCATTGTTGAGATCACAT...AGC')
    738
    Z78462.1
    Seq('CGTCACGAGGTCTCCGGATGTGACCCTGCGGAAGGATCATTGTTGAGATCACAT...CAT')
    736
    Z78461.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...TAA')
    732
    Z78460.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...TTA')
    745
    Z78459.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TTT')
    744
    Z78458.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TTG')
    738
    Z78457.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...GAG')
    739
    Z78456.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    740
    Z78455.1
    Seq('CGTAACCAGGTTTCCGTAGGTGGACCTTCGGGAGGATCATTTTTGAGATCACAT...GCA')
    745
    Z78454.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AAC')
    695
    Z78453.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GCA')
    745
    Z78452.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GCA')
    743
    Z78451.1
    Seq('CGTAACAAGGTTTCCGTAGGTGTACCTCCGGAAGGATCATTGTTGAGATCACAT...AGC')
    730
    Z78450.1
    Seq('GGAAGGATCATTGCTGATATCACATAATAATTGATCGAGTTAAGCTGGAGGATC...GAG')
    706
    Z78449.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGC')
    744
    Z78448.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGG')
    742
    Z78447.1
    Seq('CGTAACAAGGATTCCGTAGGTGAACCTGCGGGAGGATCATTGTTGAGATCACAT...AGC')
    694
    Z78446.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...CCC')
    712
    Z78445.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGT')
    715
    Z78444.1
    Seq('CGTAACAAGGTTTCCGTAGGGTGAACTGCGGAAGGATCATTGTTGAGATCACAT...ATT')
    688
    Z78443.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGG')
    784
    Z78442.1
    Seq('GTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACATAATAATTGATCGAGT...AGT')
    721
    Z78441.1
    Seq('GGAAGGTCATTGCCGATATCACATAATAATTGATCGAGTTAATCTGGAGGATCT...GAG')
    703
    Z78440.1
    Seq('CGTAACAAGGTTTCCGTAGGTGGACCTCCGGGAGGATCATTGTTGAGATCACAT...GCA')
    744
    Z78439.1
    Seq('CATTGTTGAGATCACATAATAATTGATCGAGTTAATCTGGAGGATCTGTTTACT...GCC')
    592



```python
# We can extract information as a list such as pulling out identifiers
identifiers = [seq_record.id for seq_record in SeqIO.parse("ls_orchid.gbk.txt", "genbank")]
```


```python
identifiers
```




    ['Z78533.1',
     'Z78532.1',
     'Z78531.1',
     'Z78530.1',
     'Z78529.1',
     'Z78527.1',
     'Z78526.1',
     'Z78525.1',
     'Z78524.1',
     'Z78523.1',
     'Z78522.1',
     'Z78521.1',
     'Z78520.1',
     'Z78519.1',
     'Z78518.1',
     'Z78517.1',
     'Z78516.1',
     'Z78515.1',
     'Z78514.1',
     'Z78513.1',
     'Z78512.1',
     'Z78511.1',
     'Z78510.1',
     'Z78509.1',
     'Z78508.1',
     'Z78507.1',
     'Z78506.1',
     'Z78505.1',
     'Z78504.1',
     'Z78503.1',
     'Z78502.1',
     'Z78501.1',
     'Z78500.1',
     'Z78499.1',
     'Z78498.1',
     'Z78497.1',
     'Z78496.1',
     'Z78495.1',
     'Z78494.1',
     'Z78493.1',
     'Z78492.1',
     'Z78491.1',
     'Z78490.1',
     'Z78489.1',
     'Z78488.1',
     'Z78487.1',
     'Z78486.1',
     'Z78485.1',
     'Z78484.1',
     'Z78483.1',
     'Z78482.1',
     'Z78481.1',
     'Z78480.1',
     'Z78479.1',
     'Z78478.1',
     'Z78477.1',
     'Z78476.1',
     'Z78475.1',
     'Z78474.1',
     'Z78473.1',
     'Z78472.1',
     'Z78471.1',
     'Z78470.1',
     'Z78469.1',
     'Z78468.1',
     'Z78467.1',
     'Z78466.1',
     'Z78465.1',
     'Z78464.1',
     'Z78463.1',
     'Z78462.1',
     'Z78461.1',
     'Z78460.1',
     'Z78459.1',
     'Z78458.1',
     'Z78457.1',
     'Z78456.1',
     'Z78455.1',
     'Z78454.1',
     'Z78453.1',
     'Z78452.1',
     'Z78451.1',
     'Z78450.1',
     'Z78449.1',
     'Z78448.1',
     'Z78447.1',
     'Z78446.1',
     'Z78445.1',
     'Z78444.1',
     'Z78443.1',
     'Z78442.1',
     'Z78441.1',
     'Z78440.1',
     'Z78439.1']




```python
# We can use the parse and iterator function to work with a certain record
record_iterator = SeqIO.parse("ls_orchid.fasta.txt", "fasta")
```


```python
# We can use "next" to move to the next record in the information
first_record = next(record_iterator)
```


```python
# We can print the ID of the first record
print(first_record.id)
```

    gi|2765658|emb|Z78533.1|CIZ78533



```python
# We can print the description of the first record
print(first_record.description)
```

    gi|2765658|emb|Z78533.1|CIZ78533 C.irapeanum 5.8S rRNA gene and ITS1 and ITS2 DNA



```python
# We can use the same "next" command to move to the second record
second_record = next(record_iterator)
```


```python
# We can print the ID of the second record
print(second_record.id)
```

    gi|2765657|emb|Z78532.1|CCZ78532



```python
# We can print the description of the second record
print(second_record.description)
```

    gi|2765657|emb|Z78532.1|CCZ78532 C.californicum 5.8S rRNA gene and ITS1 and ITS2 DNA



```python
# We can access the records in any order
records = list(SeqIO.parse("ls_orchid.gbk.txt", "genbank"))
```


```python
# We can print the amount of records
print("Found %i records" % len(records))
```

    Found 94 records



```python
# We can print information on a certain record such as the last record
print("The last record")
last_record = records [-1]
print(last_record.id)
print(repr(last_record.seq))
print(len(last_record))
```

    The last record
    Z78439.1
    Seq('CATTGTTGAGATCACATAATAATTGATCGAGTTAATCTGGAGGATCTGTTTACT...GCC')
    592



```python
# We can use this method to print information on any record by changing the number we ask for
print("The first record")
first_record = records[0]
print(first_record.id)
print(repr(first_record.seq))
print(len(first_record))
```

    The first record
    Z78533.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGG...CGC')
    740


## Part 2
```python
record_iterator = SeqIO.parse("ls_orchid.gbk.txt", "genbank")
```


```python
# We can use the "next" function to start at the first record
first_record = next(record_iterator)
```


```python
# We can print the first record
print(first_record)
```

    ID: Z78533.1
    Name: Z78533
    Description: C.irapeanum 5.8S rRNA gene and ITS1 and ITS2 DNA
    Number of features: 5
    /molecule_type=DNA
    /topology=linear
    /data_file_division=PLN
    /date=30-NOV-2006
    /accessions=['Z78533']
    /sequence_version=1
    /gi=2765658
    /keywords=['5.8S ribosomal RNA', '5.8S rRNA gene', 'internal transcribed spacer', 'ITS1', 'ITS2']
    /source=Cypripedium irapeanum
    /organism=Cypripedium irapeanum
    /taxonomy=['Eukaryota', 'Viridiplantae', 'Streptophyta', 'Embryophyta', 'Tracheophyta', 'Spermatophyta', 'Magnoliophyta', 'Liliopsida', 'Asparagales', 'Orchidaceae', 'Cypripedioideae', 'Cypripedium']
    /references=[Reference(title='Phylogenetics of the slipper orchids (Cypripedioideae: Orchidaceae): nuclear rDNA ITS sequences', ...), Reference(title='Direct Submission', ...)]
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGG...CGC')



```python
# We can print the annotation list of the first record
print(first_record.annotations)
```

    {'molecule_type': 'DNA', 'topology': 'linear', 'data_file_division': 'PLN', 'date': '30-NOV-2006', 'accessions': ['Z78533'], 'sequence_version': 1, 'gi': '2765658', 'keywords': ['5.8S ribosomal RNA', '5.8S rRNA gene', 'internal transcribed spacer', 'ITS1', 'ITS2'], 'source': 'Cypripedium irapeanum', 'organism': 'Cypripedium irapeanum', 'taxonomy': ['Eukaryota', 'Viridiplantae', 'Streptophyta', 'Embryophyta', 'Tracheophyta', 'Spermatophyta', 'Magnoliophyta', 'Liliopsida', 'Asparagales', 'Orchidaceae', 'Cypripedioideae', 'Cypripedium'], 'references': [Reference(title='Phylogenetics of the slipper orchids (Cypripedioideae: Orchidaceae): nuclear rDNA ITS sequences', ...), Reference(title='Direct Submission', ...)]}



```python
# We can print the annotation keys to the first record
print(first_record.annotations.keys())
```

    dict_keys(['molecule_type', 'topology', 'data_file_division', 'date', 'accessions', 'sequence_version', 'gi', 'keywords', 'source', 'organism', 'taxonomy', 'references'])



```python
# We can print the values of the first record annotations
print(first_record.annotations.values())
```

    dict_values(['DNA', 'linear', 'PLN', '30-NOV-2006', ['Z78533'], 1, '2765658', ['5.8S ribosomal RNA', '5.8S rRNA gene', 'internal transcribed spacer', 'ITS1', 'ITS2'], 'Cypripedium irapeanum', 'Cypripedium irapeanum', ['Eukaryota', 'Viridiplantae', 'Streptophyta', 'Embryophyta', 'Tracheophyta', 'Spermatophyta', 'Magnoliophyta', 'Liliopsida', 'Asparagales', 'Orchidaceae', 'Cypripedioideae', 'Cypripedium'], [Reference(title='Phylogenetics of the slipper orchids (Cypripedioideae: Orchidaceae): nuclear rDNA ITS sequences', ...), Reference(title='Direct Submission', ...)]])



```python
# We can print the source of the organism in the first record which is normally the common name
print(first_record.annotations["source"])
```

    Cypripedium irapeanum



```python
# We can also print the scientific name of the organism
print(first_record.annotations["organism"])
```

    Cypripedium irapeanum



```python
# We can create an empty variable
all_species = []
```


```python
#We can make each record in the file show the organism name
for seq_record in SeqIO.parse("ls_orchid.gbk.txt", "genbank"):
    all_species.append(seq_record.annotations["organism"])
```


```python
# We can print all species within the record
print(all_species)
```

    ['Cypripedium irapeanum', 'Cypripedium californicum', 'Cypripedium fasciculatum', 'Cypripedium margaritaceum', 'Cypripedium lichiangense', 'Cypripedium yatabeanum', 'Cypripedium guttatum', 'Cypripedium acaule', 'Cypripedium formosanum', 'Cypripedium himalaicum', 'Cypripedium macranthon', 'Cypripedium calceolus', 'Cypripedium segawai', 'Cypripedium parviflorum var. pubescens', 'Cypripedium reginae', 'Cypripedium flavum', 'Cypripedium passerinum', 'Mexipedium xerophyticum', 'Phragmipedium schlimii', 'Phragmipedium besseae', 'Phragmipedium wallisii', 'Phragmipedium exstaminodium', 'Phragmipedium caricinum', 'Phragmipedium pearcei', 'Phragmipedium longifolium', 'Phragmipedium lindenii', 'Phragmipedium lindleyanum', 'Phragmipedium sargentianum', 'Phragmipedium kaiteurum', 'Phragmipedium czerwiakowianum', 'Phragmipedium boissierianum', 'Phragmipedium caudatum', 'Phragmipedium warszewiczianum', 'Paphiopedilum micranthum', 'Paphiopedilum malipoense', 'Paphiopedilum delenatii', 'Paphiopedilum armeniacum', 'Paphiopedilum emersonii', 'Paphiopedilum niveum', 'Paphiopedilum godefroyae', 'Paphiopedilum bellatulum', 'Paphiopedilum concolor', 'Paphiopedilum fairrieanum', 'Paphiopedilum druryi', 'Paphiopedilum tigrinum', 'Paphiopedilum hirsutissimum', 'Paphiopedilum barbigerum', 'Paphiopedilum henryanum', 'Paphiopedilum charlesworthii', 'Paphiopedilum villosum', 'Paphiopedilum exul', 'Paphiopedilum insigne', 'Paphiopedilum gratrixianum', 'Paphiopedilum primulinum', 'Paphiopedilum victoria', 'Paphiopedilum victoria', 'Paphiopedilum glaucophyllum', 'Paphiopedilum supardii', 'Paphiopedilum kolopakingii', 'Paphiopedilum sanderianum', 'Paphiopedilum lowii', 'Paphiopedilum dianthum', 'Paphiopedilum parishii', 'Paphiopedilum haynaldianum', 'Paphiopedilum adductum', 'Paphiopedilum stonei', 'Paphiopedilum philippinense', 'Paphiopedilum rothschildianum', 'Paphiopedilum glanduliferum', 'Paphiopedilum glanduliferum', 'Paphiopedilum sukhakulii', 'Paphiopedilum wardii', 'Paphiopedilum ciliolare', 'Paphiopedilum dayanum', 'Paphiopedilum hennisianum', 'Paphiopedilum callosum', 'Paphiopedilum tonsum', 'Paphiopedilum javanicum', 'Paphiopedilum fowliei', 'Paphiopedilum schoseri', 'Paphiopedilum bougainvilleanum', 'Paphiopedilum hookerae', 'Paphiopedilum papuanum', 'Paphiopedilum mastersianum', 'Paphiopedilum argus', 'Paphiopedilum venustum', 'Paphiopedilum acmodontum', 'Paphiopedilum urbanianum', 'Paphiopedilum appletonianum', 'Paphiopedilum lawrenceanum', 'Paphiopedilum bullenianum', 'Paphiopedilum superbiens', 'Paphiopedilum purpuratum', 'Paphiopedilum barbatum']



```python
# We can use another method to print all of the organism names
all_species = [
    seq_record.annotations["organism"]
    for seq_record in SeqIO.parse("ls_orchid.gbk.txt", "genbank")
]
```


```python
print(all_species)
```

    ['Cypripedium irapeanum', 'Cypripedium californicum', 'Cypripedium fasciculatum', 'Cypripedium margaritaceum', 'Cypripedium lichiangense', 'Cypripedium yatabeanum', 'Cypripedium guttatum', 'Cypripedium acaule', 'Cypripedium formosanum', 'Cypripedium himalaicum', 'Cypripedium macranthon', 'Cypripedium calceolus', 'Cypripedium segawai', 'Cypripedium parviflorum var. pubescens', 'Cypripedium reginae', 'Cypripedium flavum', 'Cypripedium passerinum', 'Mexipedium xerophyticum', 'Phragmipedium schlimii', 'Phragmipedium besseae', 'Phragmipedium wallisii', 'Phragmipedium exstaminodium', 'Phragmipedium caricinum', 'Phragmipedium pearcei', 'Phragmipedium longifolium', 'Phragmipedium lindenii', 'Phragmipedium lindleyanum', 'Phragmipedium sargentianum', 'Phragmipedium kaiteurum', 'Phragmipedium czerwiakowianum', 'Phragmipedium boissierianum', 'Phragmipedium caudatum', 'Phragmipedium warszewiczianum', 'Paphiopedilum micranthum', 'Paphiopedilum malipoense', 'Paphiopedilum delenatii', 'Paphiopedilum armeniacum', 'Paphiopedilum emersonii', 'Paphiopedilum niveum', 'Paphiopedilum godefroyae', 'Paphiopedilum bellatulum', 'Paphiopedilum concolor', 'Paphiopedilum fairrieanum', 'Paphiopedilum druryi', 'Paphiopedilum tigrinum', 'Paphiopedilum hirsutissimum', 'Paphiopedilum barbigerum', 'Paphiopedilum henryanum', 'Paphiopedilum charlesworthii', 'Paphiopedilum villosum', 'Paphiopedilum exul', 'Paphiopedilum insigne', 'Paphiopedilum gratrixianum', 'Paphiopedilum primulinum', 'Paphiopedilum victoria', 'Paphiopedilum victoria', 'Paphiopedilum glaucophyllum', 'Paphiopedilum supardii', 'Paphiopedilum kolopakingii', 'Paphiopedilum sanderianum', 'Paphiopedilum lowii', 'Paphiopedilum dianthum', 'Paphiopedilum parishii', 'Paphiopedilum haynaldianum', 'Paphiopedilum adductum', 'Paphiopedilum stonei', 'Paphiopedilum philippinense', 'Paphiopedilum rothschildianum', 'Paphiopedilum glanduliferum', 'Paphiopedilum glanduliferum', 'Paphiopedilum sukhakulii', 'Paphiopedilum wardii', 'Paphiopedilum ciliolare', 'Paphiopedilum dayanum', 'Paphiopedilum hennisianum', 'Paphiopedilum callosum', 'Paphiopedilum tonsum', 'Paphiopedilum javanicum', 'Paphiopedilum fowliei', 'Paphiopedilum schoseri', 'Paphiopedilum bougainvilleanum', 'Paphiopedilum hookerae', 'Paphiopedilum papuanum', 'Paphiopedilum mastersianum', 'Paphiopedilum argus', 'Paphiopedilum venustum', 'Paphiopedilum acmodontum', 'Paphiopedilum urbanianum', 'Paphiopedilum appletonianum', 'Paphiopedilum lawrenceanum', 'Paphiopedilum bullenianum', 'Paphiopedilum superbiens', 'Paphiopedilum purpuratum', 'Paphiopedilum barbatum']



```python
# We can create an empty variable
all_species = []
```


```python
# We can also print all of the organism names within the fasta file
for seq_record in SeqIO.parse("ls_orchid.fasta.txt", "fasta"):
    all_species.append(seq_record.description.split()[1])
```


```python
print(all_species)
```

    ['C.irapeanum', 'C.californicum', 'C.fasciculatum', 'C.margaritaceum', 'C.lichiangense', 'C.yatabeanum', 'C.guttatum', 'C.acaule', 'C.formosanum', 'C.himalaicum', 'C.macranthum', 'C.calceolus', 'C.segawai', 'C.pubescens', 'C.reginae', 'C.flavum', 'C.passerinum', 'M.xerophyticum', 'P.schlimii', 'P.besseae', 'P.wallisii', 'P.exstaminodium', 'P.caricinum', 'P.pearcei', 'P.longifolium', 'P.lindenii', 'P.lindleyanum', 'P.sargentianum', 'P.kaiteurum', 'P.czerwiakowianum', 'P.boissierianum', 'P.caudatum', 'P.warszewiczianum', 'P.micranthum', 'P.malipoense', 'P.delenatii', 'P.armeniacum', 'P.emersonii', 'P.niveum', 'P.godefroyae', 'P.bellatulum', 'P.concolor', 'P.fairrieanum', 'P.druryi', 'P.tigrinum', 'P.hirsutissimum', 'P.barbigerum', 'P.henryanum', 'P.charlesworthii', 'P.villosum', 'P.exul', 'P.insigne', 'P.gratrixianum', 'P.primulinum', 'P.victoria', 'P.victoria', 'P.glaucophyllum', 'P.supardii', 'P.kolopakingii', 'P.sanderianum', 'P.lowii', 'P.dianthum', 'P.parishii', 'P.haynaldianum', 'P.adductum', 'P.stonei', 'P.philippinense', 'P.rothschildianum', 'P.glanduliferum', 'P.glanduliferum', 'P.sukhakulii', 'P.wardii', 'P.ciliolare', 'P.dayanum', 'P.hennisianum', 'P.callosum', 'P.tonsum', 'P.javanicum', 'P.fowliei', 'P.schoseri', 'P.bougainvilleanum', 'P.hookerae', 'P.papuanum', 'P.mastersianum', 'P.argus', 'P.venustum', 'P.acmodontum', 'P.urbanianum', 'P.appletonianum', 'P.lawrenceanum', 'P.bullenianum', 'P.superbiens', 'P.purpuratum', 'P.barbatum', 'C.irapeanum', 'C.californicum', 'C.fasciculatum', 'C.margaritaceum', 'C.lichiangense', 'C.yatabeanum', 'C.guttatum', 'C.acaule', 'C.formosanum', 'C.himalaicum', 'C.macranthum', 'C.calceolus', 'C.segawai', 'C.pubescens', 'C.reginae', 'C.flavum', 'C.passerinum', 'M.xerophyticum', 'P.schlimii', 'P.besseae', 'P.wallisii', 'P.exstaminodium', 'P.caricinum', 'P.pearcei', 'P.longifolium', 'P.lindenii', 'P.lindleyanum', 'P.sargentianum', 'P.kaiteurum', 'P.czerwiakowianum', 'P.boissierianum', 'P.caudatum', 'P.warszewiczianum', 'P.micranthum', 'P.malipoense', 'P.delenatii', 'P.armeniacum', 'P.emersonii', 'P.niveum', 'P.godefroyae', 'P.bellatulum', 'P.concolor', 'P.fairrieanum', 'P.druryi', 'P.tigrinum', 'P.hirsutissimum', 'P.barbigerum', 'P.henryanum', 'P.charlesworthii', 'P.villosum', 'P.exul', 'P.insigne', 'P.gratrixianum', 'P.primulinum', 'P.victoria', 'P.victoria', 'P.glaucophyllum', 'P.supardii', 'P.kolopakingii', 'P.sanderianum', 'P.lowii', 'P.dianthum', 'P.parishii', 'P.haynaldianum', 'P.adductum', 'P.stonei', 'P.philippinense', 'P.rothschildianum', 'P.glanduliferum', 'P.glanduliferum', 'P.sukhakulii', 'P.wardii', 'P.ciliolare', 'P.dayanum', 'P.hennisianum', 'P.callosum', 'P.tonsum', 'P.javanicum', 'P.fowliei', 'P.schoseri', 'P.bougainvilleanum', 'P.hookerae', 'P.papuanum', 'P.mastersianum', 'P.argus', 'P.venustum', 'P.acmodontum', 'P.urbanianum', 'P.appletonianum', 'P.lawrenceanum', 'P.bullenianum', 'P.superbiens', 'P.purpuratum', 'P.barbatum']



```python
# We can use the record iterator and parse to modify information
record_iterator = SeqIO.parse("ls_orchid.fasta.txt", "fasta")
```


```python
first_record = next(record_iterator)
```


```python
# We can print the ID of the first part of the file
first_record.id
```




    'gi|2765658|emb|Z78533.1|CIZ78533'




```python
# We can change the ID of the first file to a new name
first_record.id = "new_id"
```


```python
first_record.id
```




    'new_id'




```python
first_record
```




    SeqRecord(seq=Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGG...CGC'), id='new_id', name='gi|2765658|emb|Z78533.1|CIZ78533', description='gi|2765658|emb|Z78533.1|CIZ78533 C.irapeanum 5.8S rRNA gene and ITS1 and ITS2 DNA', dbxrefs=[])




```python
# We can change the description of the record
first_record.description = first_record.id + " " + "Mutations induced randomly"
```


```python
# We can print out the record with the new name and description
print(first_record.format("fasta"[:200]))
```

    >new_id Mutations induced randomly
    CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGGAATAAA
    CGATCGAGTGAATCCGGAGGACCGGTGTACTCAGCTCACCGGGGGCATTGCTCCCGTGGT
    GACCCTGATTTGTTGTTGGGCCGCCTCGGGAGCGTCCATGGCGGGTTTGAACCTCTAGCC
    CGGCGCAGTTTGGGCGCCAAGCCATATGAAAGCATCACCGGCGAATGGCATTGTCTTCCC
    CAAAACCCGGAGCGGCGGCGTGCTGTCGCGTGCCCAATGAATTTTGATGACTCTCGCAAA
    CGGGAATCTTGGCTCTTTGCATCGGATGGAAGGACGCAGCGAAATGCGATAAGTGGTGTG
    AATTGCAAGATCCCGTGAACCATCGAGTCTTTTGAACGCAAGTTGCGCCCGAGGCCATCA
    GGCTAAGGGCACGCCTGCTTGGGCGTCGCGCTTCGTCTCTCTCCTGCCAATGCTTGCCCG
    GCATACAGCCAGGCCGGCGTGGTGCGGATGTGAAAGATTGGCCCCTTGTGCCTAGGTGCG
    GCGGGTCCAAGAGCTGGTGTTTTGATGGCCCGGAACCCGGCAAGAGGTGGACGGATGCTG
    GCAGCAGCTGCCGTGCGAATCCCCCATGTTGTCGTGCTTGTCGGACAGGCAGGAGAACCC
    TTCCGAACCCCAATGGAGGGCGGTTGACCGCCATTCGGATGTGACCCCAGGTCAGGCGGG
    GGCACCCGCTGAGTTTACGC
    


## Part 3
```python
from Bio.Seq import Seq
```


```python
from Bio.SeqRecord import SeqRecord
```


```python
# We can input a new record along with its ID and description
rec1 = SeqRecord(
    Seq("MMYQQGCFAAGGTVLRLAKDLAENNRGARVLVVCSEITAVTFRGPSETHLDSMVGQALFGD"
       "GAGAVIVGSDPDLSVERPLYELVWTGATLLPDSEGAIDGHLREVGLTFHLLKDVPGLISK"
       "NIEKSLKEAFTPLGISDWNSTFWIAHPGGPAILDQVEAKLGKEEKMRATREVLSEYGNM"
       "SSAC",
),
    id = "gi|14150838|gb|AAK54658.1|AF376133_1",
    description = "chalcone synthase [Cucumis sativus]",
)
```


```python
# We can input a second new record along with its ID and description
rec2 = SeqRecord(
    Seq("YPDYYFRITNREHKAELKEKFQRMCDKSMIKKRYMYLTEEILKENPSMCEYMAPSLDARQ"
        "DMVVVEIPKLGKEAAVKAIKEWGQ",
),
    id = "gi|13919613|gb|AAK33142.1|",
    description = "chalcone synthase [Fragraria vesca subsp.bracteata]",
)
```


```python
# We can input a third new record along with its ID and description
rec3 = SeqRecord(
    Seq("MVTVEEFRRAQCAEGPATVMAIGTATPSNCVDQSTYPDYYFRITNSEHKVELKEKFKRMC"
        "EKSMIKKRYMHLTEEILKENPNICAYMAPSLDARQDIVVVEVPKLGKEAAQKAIKEWGQP"
        "KSKITHLVFCTTSGVDMPGCDYQLTKLLGLRPSVCRFFMMYQQGCFAGGTVLRMAKDLAEN"
        "NKGARVLVVCSEITAVTFRGPNDTHLDSLVGQALFGDGAAAVIIGSDPIPEVERPLFELV"
        "SAAQTLLPDSEGAIDGHLREVGLTFHLLKDVPGLISKNIEKSLVEAFQPLGISDWNSLFW"
        "IAHPGGPAILDQVELKLGLKQEKLKATRKVLSNYGNMSSACVLFIDLEMRKASAKEGLGT"
        "TGEGLEWGVLFGFGPGLTVETVVLHSVAT",
),
    id = "gi|13925890|gb|AAK49457.1|",
    description = "chalcone synthase [Nicotania tabacum]",
)
```


```python
# We can put all three of our peptide records into one file
my_records = [rec1, rec2, rec3]
```


```python
# We can create a new file with the three peptide records
SeqIO.write(my_records,"myexample.faa", "fasta")
```




    3




```python
records = SeqIO.parse("ls_orchid.gbk.txt", "genbank")
```


```python
# We can convert the records from the Genbank file to a fasta file
count = SeqIO.write(records, "my_example.fasta", "fasta")
print("Converted %i records" % count)
```

    Converted 94 records



```python
# We can create a loop to print all reverse complements from the genbank file
for record in SeqIO.parse("ls_orchid.gbk.txt", "genbank"):
    print(record.id)
    print(record.seq.reverse_complement())
```

    Z78533.1
    GCGTAAACTCAGCGGGTGCCCCCGCCTGACCTGGGGTCACATCCGAATGGCGGTCAACCGCCCTCCATTGGGGTTCGGAAGGGTTCTCCTGCCTGTCCGACAAGCACGACAACATGGGGGATTCGCACGGCAGCTGCTGCCAGCATCCGTCCACCTCTTGCCGGGTTCCGGGCCATCAAAACACCAGCTCTTGGACCCGCCGCACCTAGGCACAAGGGGCCAATCTTTCACATCCGCACCACGCCGGCCTGGCTGTATGCCGGGCAAGCATTGGCAGGAGAGAGACGAAGCGCGACGCCCAAGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAAGACTCGATGGTTCACGGGATCTTGCAATTCACACCACTTATCGCATTTCGCTGCGTCCTTCCATCCGATGCAAAGAGCCAAGATTCCCGTTTGCGAGAGTCATCAAAATTCATTGGGCACGCGACAGCACGCCGCCGCTCCGGGTTTTGGGGAAGACAATGCCATTCGCCGGTGATGCTTTCATATGGCTTGGCGCCCAAACTGCGCCGGGCTAGAGGTTCAAACCCGCCATGGACGCTCCCGAGGCGGCCCAACAACAAATCAGGGTCACCACGGGAGCAATGCCCCCGGTGAGCTGAGTACACCGGTCCTCCGGATTCACTCGATCGTTTATTCCACGGTCTCATCAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78532.1
    GCCTCAACTCAGCGGGTGGCCCCGCCTGACCTGGGGTCGCATCTGAATGGAAATCAACTGCCCAATGGTTATTTTAGCTCCATTGGGGTTCAATTAGGTTCTTGTGTAGGTTCGAAAAAATACAACAACATGGGGGATTCAAATAGCAGCCTTATGACTGTTAGCATTCTCCACCTCGTGCCACATTCCTACCCATCAAAGCAACAATCCTTAGACCCACCGCACCTAGGCACAAGGGGCCAATCATTCACATCCGTATAATGCCAGCTTAGCGATATGCCAAGCAAGCATTGGTAGGAGAGACGCAACACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGAGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCGCATTTAGCTGCGTTCTTCATCGATGTTAGAGCCAAGATATCCGTTGCGAGAGTCATAAAAATTCACTGGCATGCTCAGTAGCATACTGCCCCTCTGATTTTTTCTGACAATAATGTCATTCATCAGTGATGCTTTATATGACTTGGCGCAAACTGCACCGTACTAAAGTTCAAACCTGCCATGAAAGCTCTTGAGGAGGCCCAACAACAAAGCAGGGTCACGACAAAAGCAGTGCCACGACGAGCTGAGTTACCACAGGTCCTCCAGATTCACTCGATCATATATTCTGTTGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78531.1
    TTACTCACGGGTGGCCCGCCTGACTGGGGTCGCATCTGAATGGAAATCACCGCCCGAGGGCGGTTTTGCGTCCACTGGGGGTTCAAACAGGTTCTTCTGTAGGCCTGACGAGTACGACAACACGGGGGATTCGAACGGCAGCCTTGCGGCTGCCAGCATCCTCCACCTACTGCCAAGTTCCGGGCCATCAGAACACCGATGCTTAGACCCGCCGCACCTAGGCGCAAGGGGCCAATCTTTCACGTCCTCACAATGACAGACTAGCCGCATGCCAATCAAGCATTATCAGGAGAGACGCAGCACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCGCGGGATTCTGCAATTCACACCACTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCAGATATCCCGTTGCCGAGAGTCGTCATAAATTCATTGGCACGCGCAACAGACGCCGCCCCTCCGAGTTTTTCTTGACAAAAATGCCATCCATCGGTGACGCTCCATATGACTTGGCGCAAACTGCGCCGCGCTAGACGTTCAAACCCGCCATGAAAGTTCCCGAGGCGGCCCGGTCGCAAACCGGGTTCACCACGAGAGCAAAGCCACGGTGAGCCGTGTAACCACGGGTCCTCCGGATTCACTCGATCGTATGTTCTGCTGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78530.1
    ATGGGGTGGCCCCCTGACTGGGGTCGCATCTGATTGGAAATCAACCACCCAAGGATTGTTTTGCCTTAATTAGGGTCCAAACAAGTTGTTCTATAGGCCCAACAAACATGACAACATGGGGGATTCAAACAATAGCCTTGCGGCTGCCAGCATCCTCCACCTCTTGCCAAGTTTCAGACCATCAAAACACATATCCTTAGACCCACCGCACCTAAGCACAAGGGGCCAATCTTTCACATCCACACAATGATGGCCTAGCTATATGCTGGACAAGCATTGGAAGGAGAGATAAAACATACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTCGGTGCGTTCTTCATCGATGCAAGAGTCAAGATATCCGTCTGCCGAGAGTCATCAAAATTCATTGGAATGAACATTAACACACCACACCTCCGACGTTGTGTGGGGCAATAATGTCATTCATCGGTGATTCTTTATATACCTTGGTGCAAACTGGACCGTGCTAGAGGTTCAAACCCGCCATGAAAGCTCTCAATGAGGCCCAATGACAAATCATGGTCACCACAAAAGGATATCCCTAGCGAGCCAAATTACCACAAGTCCTCCAGATTCACTCAATCGTTTATTATGTTGTTTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78529.1
    TTTTAGCTCAGCTGGTGGCCCGCTGACTGGGGTCGCATCTGATTGGAAATCAACCACCCAAGGATTGTTTTGCCTTAATTAGGGTCCAAACAAGTTGTTCTATAGGCCCAACAAATATGACAACATGGGGGATTCAAACAATAGCCTTGCGGCTGCCAGCATCCTCCACCTCTTGCCAAGTTTCAGACCATCAAAACACATATCCTTAGACCCACCGCACCTAAGCACAAGGGGCCAATCTTTCACATCCACACAATGATGGCCTAGCTATATGCTGGACAAGCATTGGAAGGAGAGATAAAACATACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTCGGTGGGATTCTTCATCGATGCAAGAGCCAAGATATCCCTCGGCGAGAGTCATCAACAATTCATTGGCATGACCAATAACACACCACACCTCCGACTTTTTTGACAATAATGTCATTCATCGGTGATTCTTTATATACCTTGGTGCAAACTGCACCGTGCTAGAGGTTCAAACCCGCCATGAAAGCTCTCAATGAGGCCCAATGACAAATCATGGTCACCACAAAAGGAAATCCCTAGCGAGCCAAATAACCACAAGTCCTCCAGATTCACTCAATCGTATATTCTGCTGTCTCAACAATGTCCTTCGGCAGCTCGCCGT
    Z78527.1
    GGGTGGCCCGCCTGACCTGGGGTCGCATCTAAATGGAAATCAACCGCCAAGGGTCATTTTACAATCCATTGGGGTTCAAGCATGTTCTTTTATAGGTTCGACAAATATGACAACATGGGGGATTCGAACGACAGCCTTGCGGCTTCCAGCATCCTCCACCTCCTGCCAGGTTTCGATCCATCAAAACATCAATCCTTAGACCCACCGCACCTAGGCACAAGGGGCCAATCTTTCACATCCACCCAATGCCAGCCTAGTTGTATGCCGGGTAAGCATTGGCAAGAGAGATGCGATACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCGCGGGATTCTGCAATTCACACCACTTATCGCATTTCGCTGCGTTCCTTCATCGATGCAAAGAGCCAGATATCCGTTGGCCGAGAGTCATCAAAAATCCATTAGCGCATGCAATAACACACTGCCACTCCAACTTTTGGTGACAGTAGTGCCATTTGCAGTTATGCTTCATATGACTTGGCGCAAACTGCACCGTGGTAGAGGTTCAAACCCGCCATGTAAGCTCCTAAGGAAGCCCAACAACAAATCAGGCGAGCCGAGTAACCACAGGTCATCCAGATTCACTCGATCATATATTCTACTGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78526.1
    ACATAAACTCAGCGGGTGGCCCCGCCTGACCTGGGGTCGCATCTAAATGGAAATCAACCGCCAAGGGTCATTTTACATCCATTGGGGTTCAAGCATGTTCTTTTATAGGTTCGACAAATATGACAACATGGGGGATTCGAACGACAGCCTTGCGGCTTCCAGCATCCTCCACCTCCTGCCAGGTTTCGATCCATCAAAACATCAATCCTTAGACCCACCGCACCTAGGCACAAGGGGCCAATCTTTCACATCCGCACAATGCCAGCCTAGTTGTATGCCGGGTAAGCATTGGCAAGAGAGATGCGATACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCGCGGGATTCTGCAATTCACACCACTTATCGCATTTCGCTGGGTTCTTCATCGATGCAAGAGCCAAGATATCCGTTGGCGAGAGTCATCAAAAATCCATTAGCGCATGCAATAGCACACTGCCACTCCAACTTTTGGTGACAATAATGCCATTGTTCAGTTATGCTTCATATGACTTGGCGCAAACTGCACCGTGGTAGAGGTTCAAACCCGGCATGTAAGCTCCTAAGGAAGCCCAACAACAAATCAGGGGAGCCGAGTAACCACAGGTCCTCCAGATTCACTCGATCATATATTCTACTGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78525.1
    TGCAGCGGGTGGCCCGCCTGACCTGGGGTCGCATATGAATGGCAGTCAACCATCCAAGGGTCATTTTGCCTCCACTGGGGTTCAAACAAGTTATTATATAGGTCCGATAAATACGATAACATGGGGGATTTAAATGACAACCTTGTGGCAGCCAATGTCCTCCACCTCCTGCCAAGTTTCGGGCCATCAAAACATTCCTTAGACCCAACGCGCCAAGGCACAAGGGGCCAATCTTTCGCATCCGCACAATGCAAGCCTAGATATATGGACAGGCATTGGCAAGAGAGACACAACACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAAGACTCGATGGTTCACGGGATCTTGCAATTCACACCACTTATCGCATTTCGCTGGGTTCTTCATCCGATGCAAGAGCCAAGATATCCGCTTGTCGAGAGTCATCAAAATTATATTCGGCATGCACAGGAGGACACCGCCCCTCCAACTTTTTTGACAATAATGTCATTTGTCGGTGATGCTTCATATGACTTGGCACAAACTGTGCCGTGCTAGAGTTTGAAACCTGCCATGAAAGCTCCCGAGGAGGCCCAACGACAAATTTGGGTCACCACAAAAGCAAAGCCCTCGGCAAGCCGAATAACCACAGGTCCTCCGGATTCACTCGATGTATATTCTGCTATCTCAACA
    Z78524.1
    GCTTAATCTCAGGGTGGCCACCTGACCTGGGGTTGCATCTGAATGAAAATCAACCGCCCAAGGGTCATTTTTCCTTCATTGGGGTTCAAACAAGTTCTTTTATAGGTTCAACAAATACGACAACATGGGGAATTCAAATGGCAGCCTTGCAGCTACCAACATTCTCCACCTCCTGCCGAGTTCCGAGCTATCAAAATACTAATCCTTAGACCCACCGCACCTAAGCACAAGGGGCCAATATTTCATCCGCACAATGCTAACCTAGCTATATGCTGGGCAAGCATTGGNAGGAGAGATGCAGCACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTTGGGCGCAACTTGCGTTCAAATACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCGCATTTCGCTGCGTTCTTCATCGGATGTAAAGAGCCAAGTTCCGTTGCGAGAGTCATCAAAAAATTCATTAGCATGCACAACAGCATGCCGCCCCTCCGACTTTTTTAACAACAATGCCATTCATCAGGGATGCTTATATGACTTGGCACAAACTGCGCCGTGCTAAAGGTTTGAACCCGCCATGAAAGCTCCTAAGGAGGCCCAATTAAGGTCACCACAAAAGCAAAGCACTGACGAGCCAAGTAACCACATGTCCTCCATATTCACTCAATCATATATTCTACTATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78523.1
    CTTAAACTCAGCGGTGGCCCCGCCTGACCTGGGTCGCAATATGAATGGCAATCAACCGCCCGAGGGTTGTTTGCCTCCATTGGGATTCAAACAGGTTCTTCTCTAGGTCCGACAAACACGACAATATGGGGGGTTCAAACGATAGCCTTATAGCTGCCAACATCCTCCACCTCCTGCCAGGTTTCGAACCATCAAAACACTAGTCCTTAGACCCACCGCACCTAGGCACAAGGGGCCAATCTTTCGCATCCACGCAATGCAAACCAATCTATATGATAGGAAAGCATTGATAGGAGAGACGCAGCACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCATCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACATATGTCGAGAGTCAATAAAAAATTCATTATGCATGCATGCAAGAGCACGCTTCCACTCTAGCTTTTGGGACAATAATGTCATTCCTCAGTGATGCTTCATATGACTTGGCGCATACTGCACCGTGCTAGAGGTTCACACCCACAAGGAAAGCTTGGGAGGAGGCCCAATGACAAGTTAGGGTCACCGCAAAAGCAAAGCCTATGTCGAGCTGAGTAACCACAAGTCCACCGGATTCACTCGATCATATATTCTGCTGTCTCAACAATGATCCTGCCGCAGGTTCACCTACGGAAACCTGGTTACG
    Z78522.1
    CTCAGCGGGTGGCCGCCTGACCTGGGGTCGCATATGAATGGCAATCAACCGCCCGAGGGTTGTTTGCCCCCATTGGGATTCAAACATGTTCTTCTCTAGGTCCGACAAACACGACAATATGGGGGATTCAAATGATAGCCTTATAGCTGCCAACATCCTCCACCTCCTGCCAGGTTTCGAACCATCAAAACACTAGTCCTTAGACCCACCGCACCTAGGCACAAGGGGCCAATCTTTCACATCCACGCAATGCAAACCTGTCTATATGACGTATAACCATTGACAGGAGAGACGCAGCACACGACGCCCAGGCAGGCGTTCCCTTAGCCTGATGGCATCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCGATCACACCACGATATGTCGAGAGTCATAAAAAATTCATTTGCATGTATGCAATAGCACGCTTCCACTCCAACTTTTTGGACAATAATGTCATTCATCAGTGATGCTTCATATGACTTGGCGCATACTGCACCGTGCTAGAGGTTCAAACCCACAAGGAAAGCTTGGGGGGAGGCCCAATGACAAATTAGGGTCACCGCAAAAGCAAAGCCTATGTCGAGCTGAGTAACCACAAGTCCACCGGATTCACTCGATCATATATTCTGCTGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78521.1
    GGTGGCCCCGCCTGACCTGGGGTCGCATATGAATGGCAATCAACCGCCCGAGAGTTGTTTGCCTCCATTGGGATTCAAACATGTTTTTCTCTAGGTCCGACAAACACGACAATATGGGGGATTCAAATGATAGCCTTATAGCTGCCAACATCCTCCACCTCCTGCCAGGTTTCGAACCATCAAAACACTAGTCCTTAGACCCACCGTACCTAGGCACAAGGGGCCATTCTTTCACATCCACGCAATGCAAACCTATCTATATGACATGAAAGCATTGACAGGAGAGACGCAGCACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCATCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCGCATTTTGCTGCGTTCTTCATCGATGCAAGAGCCAAGATATCCGTTGTCGAGAGTCATAAAAAATTCATTTGCATGCATGCAATAGCACGCTTCCACTCCAACTTTTTGACAATAATGTCATTCATCGGTGATGCTTCATATGACTTGGCGCATACTGCACCGTGCTAGAGGTTCAAACCCACAAGGAAAGCTTTGGAGGAGGCCCAATGGCAAATTAGGGTCACCGCCAAAGCCAAGCCTATGTCGAGCTGAGTAACCACAAGTCCACCGGATTCACTCGATCATATATTCTACTGTCTCAACAATGATCCTTCCGCAGGTTCACCTAC
    Z78520.1
    AAACTCAGCGGGTGGCCCGCCTGACCTGGGGTCGCATATGAATGGCAATCAACCGCCCGAGGGTTGTTTGCCTCCATTGGGATTCAAACAGGTTCTTCTCCAGGTCCGACAAACACGACAATATGGGGGATTCACATGACAACCTTATAGCTGCCAACATCCTCCACCTCCTGCCAGGTTTCGAACCATCAAAACACTAGTCCTTAGACCCACCGCACCTAGGCACAAGGGGCCAATCTTTCACATCCACACAATGCAAACCTATCGATATGACAGGAAAGCATTGACAGGAGAGACGCAGCACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCATCGGGCGCAACTTGCGTTCACATACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACAATTTCGCTGCGTTCTTCATCCGATGCAAGAGCCAAGATATCCCTTGTCGAGAGTCATAAAATATTCTATTGCATGCATGCAATAGCACGCTTCTACTCCAACTTTTTGGACAATAATGTCATTCATCGGTGATGATTCATATGACTTGGCGCATACTGCACCGTGCTAGAGGTTCAAACCCACAAGGAAAGCTTTGGAGGAGGCCCAATGACAAATTAGGGTCACCGCAAAAGCAAAGCCTATGGCGAGCTGAGTAACCACAAGTCCACCGGATTCACTCGATCATATATTCTGCTGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78519.1
    TAAACTCAGCGGGTGGCCCCGCCTGACCTGGGGTCGCATATGAATGGCAATCAACCGCCCGAGGGTTGTTTGCCTCCATTGGGATTCAAACATGTTCTTCTCTAGGTCCAACAAACGCGACAATATGGGGGATTCAAATGATAGCCTTATATAGCTGCCAACATCCTCCACCTCCTGCCAGGTTTCGAACCATCAAAACATTAAGTTCTTAGACCCACCGCACCTAGGCACAAGGGGCCAATCTTTCACATCCACACAATGCAAACCTATCTATATGATGGGAAAGCATTGACAGGAGAGACGCAGCACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCATCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCAAGATATCCGTTGTCGAGAGTCATAAAAAAATTCATTTGCATGCATGCAGTAGCACGCTTCCACTCCAACTTTTGGACAATAATGTCATTCATCGGCAATGCTTCATATGACTTGGTGCATACTGCACCGTGCTAGAGGTTCAAACCCACAAGGAAAGCTTGGGAGGAGGCCCAATGACAAATTAGGGTCACCGCAAAAGCAAAGCCTATGGCGAGCTGAGTAACCACAAGTCCACCAGATTCACTCGATCATAT
    Z78518.1
    GGAAAACTCAGCGGGTGGCCCCGCCTGACCTGGGGTCGTATCTGAATGGCAATCAATCGCTCAAGGGTTGTTTTGCCTCCATAGGGGTTCAAATAGGTTCTTCTGTGGGTCCGGCAAATACGACAACATGTGGGATTTGAACGACAGCCTTGCGACTATCAACATGCTCCACCTCCTGCCGGGGTTTGGGCCATCAAAACACCAATCCTTAGACCCACCGCACCTAAGCACAAGGGGTCAATCTTTCATATCTATGCAATACTGGCCTAATTGTATGTCGAGCAAGCATTGGCAAAAGAGATGCAGCACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTTGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGTATGCGTAACACACACCGTCCATCCGAGTTTGGGGGACAATAATGTAATTCGTCGGCGATGCTTCATATGACTTGGCGCAAACTGCGCCGTGATAGAGGCTCAAACCCGCCATAAAAGCTCTCGAGGATGCCCAACTACAAATCAGGGTCACCACAAAAGTTAAGCCTTCGACGAGCCGAGTAACCACAAGTCCTCCGGATTCACTCGATCGAATATTCTACTATCTCAACAATGATCCTCCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78517.1
    GCTTAAACTCAGCGGGTGGCCCCGCCTGACCTGGGGTCGTATCTGAATGGCAATCAATCGCTCGAGGGTTGTTTTGCCTCCATTGGGGTTCAAATAGGTTCTTCTGTGGGTCCGACAAATACGACAACATGTGGGATTTGAATGACAGTCTTGCGACTATCAACATGATCCACCTCCTGCCGGGGTTCGGGCCACCAAAACACCAATCCTTAGACCCACCGCACCTAAGCACAAGGGGTCAATCTTTCATATCTATGCAATACTGTCCTAATTGTATGTCGGGCAAGCATTGGCAAAAGAGATGCAGCACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTTGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCGCATTTCGCTGCGTTCTTCATCGATGGAAGAGCAAGATATCCGTTGGCCGAGAGTCATCAAAAATTTATTGGCATGCACAACAGCACACCGCCCATCCGACTTTTGTGACAATAATGTCATTCATCGGTGATGCTTCATATGACTTGGCGCAAACTGCACCGTGATAGAGGTTCAAACCCGCCATAAAATCTCCTGAGGATGCCCAACGACAAATCAGGGCCACAAAAGCTAAGACTTGGACGAGCCGAGTAACCACAAGTCCTCTGGATTCACTCGATCGTATATTCTACTGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78516.1
    TTAAACTCAGCGGGTGGCCCCGCCTGACCTGGGGTCATATCTGAATGGCAATCAATCACTCGAGGGTTGTTTTGCCTCCATTGGGGTTCAAATAGGTTCTTCTGTGGGTCCGACAAATACGACAACATGTGGGATTTGAACGATAGTCTTGCGACTATCAACATGCTCCACCTCCTGCCGGGGTTCGGGACATCAAAACACCAATCCTTAGACCCACCGCACCTAAGCACAAGGGGTCAATCTTTCATATCTATGCAATACTGGCCTAATTGTATGTCGGGCAAGCATTGGCAAAAGAGATGCAGCACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTTGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCGCATTTCGCTAAGATATCCGTTGACGAGAGTCATCAAAAAATTATTGGTATGCACAATAGCACACCGCCCATCCGACTTTGGTGACAATAATGTCATTCGTCGGTGATGCTTCATATGACTTGGCGCAAACTGCGCCGTGATAGAGGTTCAAACCCGCCATAAAAGCTCCTGAGGATGCCCAACGACAAATCAGGGCCACAAAAGCTAAGACTTGGACGAGCCGAGTAACCACAAGTCCTCCGGATTCACTCGATCATATATTATACTGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78515.1
    GCTTAAACTCAGCGGGTGGCCTCACCTGACCTGGGGTTGCTTCCAAATGGCAATCAACCGCCCATGGGTTGGGTTGATGTGGCACCAATGGGGTTCCGAAGGGTCCTACGAATTATGATGGCCCATCAAATCGGACAACATACGGCATTCGCACAACAGCTTTGGCTCAGCCATTGCCCATCCACCTCTTGTCGAATTTTGGACCATCAAAAGCCCGAGGTCTTAGACCCATTGAACCGAAGAACAAGGGGCCAAACTCACATCCGCACCGACCGAACGGTATTCTTGGACAACCATTGGTGGGAGAGACACAGTACACAACGCCCAGGCAGGCGTGCCCTTAGCCTGACGGCCTCGAGCGCAACTTGCGTTCAAAGACTCGATGTTCACGGGATTCTGCAATTCACACCACTTATCGCATTTCGCTGCGTCCTCCATCGATGCAAGAGCCGAGATATCCGTTGCCGAGAGTTGTAAAAATGTAACTTGGGGGGCACGAATCAAGTGCCACCCCCACCGGCTACAAGGGAACACCCCTATGCCGATTCTTGTCTCAAAAGACTTGGCGCGAAGCTGCGCCGGGCTAGAAGTTTGATCGCCTTCAAGAACACTCCCAAGGAGGCCCGATGGCAGATCAAAGGCCACCACAGTAGCGAAAGCCCCTGGGAGATCGAAGTACACCGGTCCTCCAGATTAACTCGATCGTTTATTCTACGGTCTCAGCAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78514.1
    TAGCGGGTAGTCTCACCTGACCTGGGGTTGCATCCTAATGGCCGTCAACCGCCCATGGGTTGGGTCGACGTGCCTCCGACGGGGTTCGAAAGGGATCATTCCGGCCCATCAAATACGACAACATAGGGCATTCGCACAACAGCTTTGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTGTGGACCGTCACAAGCCCAAGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACTCTCATCCTCGCCGGCCGAACAGTATGGCTGTTCAACCTTAGCATTGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTCGCTGGGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGACGAGAGTTGTCAGACGAATCTTAAGGGGGCACAGCAGCACGCCGGCCCTCCCCGGTCCATCCATGGAGGACGAAGACCCTTGTCGGTTCTATCTCATTTGACTTGGAGCGAAACTGCGCCGGGCCAAGGGTTCAATTGCCATAGAGAACTCTCCCAAGGAGGTTCGATGGGAAATCATGGTCACCACAGCGGCGGAAGCCCCTGGGAGACCAAAGTACACCGGTCCTCCGGATTATCTCGATCGTATATTTGGGGGCTTCAAAAATGATCCTCCCGAAGGTCCACCTACGGAAACCTTGTTACG
    Z78513.1
    CTCAGCGGTAGCTCACTGACTGGGTTGCATCCTAATGGCGTCACCGCCCATGGGTTGGGTCGACGTGCCTCCGAAGGGGTTCGAAAGGGATCGTTCTGGCACATCAAATACGACAACATAGGGCATTCGCACAAAAGCTTCGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTGTGGACCGTCACAAGCCCAAGTCTTGGACCCATCGCACAGAAGAACAAGGGGCCAAACTCTCATCCTCGCCGGCCGAACAGTATGGCTGTTCAACCTTAGCATTGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCCGAGAGTTGTCGAAAATTCTTTCGGGCACAGCAGCATGCCGGCCTCCCCGGCTCCATGGAGGACGATGCCCTTGCCGGTTCTATCTCATTTGACTTGGCGCGAAACTGCGCCGGGCCAAGGGTTCAATTGCCATGGAGAAATCTCCCAATGAGGCTCGATGGCAAATCATGGTCACCACAGCGGCGGAAGCCCCTGGGAGACCAAAGTACACCGGTCCTCCGGATTAACTCGATCGTGTATTTGGCGGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78512.1
    GCTTAAACTCAGCGGGTGGCCTCACCTGACCTGGGGTTGGATCCAAATGACCGTCGACCGCCCACGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCCGGCCCATCAAATACGACAACGCAGGGCATTCGCACGACAGCTTCGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTCCGGCCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACACTCATCCGCGCCGGCCGAACAGTATGCCTGTTCAGCCTTAGCAATGGGAGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCAGCTTATCACATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCGAGAGTTGTCAAAAAATAATTTCATTGGGGGAACGGAAGAACGCCGGCCCTCCCCGGTTCCATGGGGGACGATGCCCTCCGCCGGTTCCATCTCATTTGATTTGGGGCGAAACTGCGCCGGGCCAAGGGTTTCATTTGCCATCAAGAATTCCTCCCAGAGGGTTCGACGGAAATTCACGGTCACCACAGTAGCCAAGGCCCCTGGGAGACCAAACTACACCGGCCCTCCGGATTAATTCGATCGTTTTTTTGGGGGCTTCAAAAATGATCCTCCCGAAGGTCCACCTACGGAAACCTTGTTACG
    Z78511.1
    TCCTCACCTGACCTGGGGTTGCATCCAAATGGCCGTCGACCGCCCACGGGTTGACGTGCCTCCAATGGGGCTCAAAAGGGATTTATTCCGGCCCATCAAATACGACAGCGTAGGGCATTCGCACGACAGCTTCGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTCCGGCCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACACTCATCCGCGCCGGCCGAACAGTATGCCTGTTCAGCCTTAGCAATGGGAGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTCGCTGGGTTCTTTCATCGGATGCAAAGAGCCGAGATATCCGTTGCCGAGAGTTGTCAAAAAAAAATTCATGGGGGGAACGGAAGAACGCCGGCCCTCCCCGGTTCCATGGAGGACGATGCCCTCCGCCGGTTCCATCTCATTTGACTTGGGGCGAAACTGCGCCGGGCCAAGGGTTCAATTGCCATCAAGAATTCTCCCAAGGAGGTTCGACGGAAATTCACGGCCACCACAGTAGCCAAAGCCCCTGGGAGACCAAACTACACCGGTCCTCCGGATTAACTCGATCGTTTTTTTGGGGGTCTCAACAATGATCCTTCCGAAGGTTCACCTACGGAAACCTTGTTACG
    Z78510.1
    TAAACTCAGCGGGTGGCTCACCTGACCTGGGGTTGCATCCAAATGGCCGTCAACCGCCCATGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCCGGCCCATCAAAAACGACAACGCAGGGCATTCGCACGACAGCTTTGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTTGGACCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAATCTCTCATCCGCGCCGGCCGAACAGTATGCCCGTTCAACCTTAGCAATGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTCGCTGCGTTCTTCCATCCGATGCAAAGAGGCGAGATATCCGTTGCGAGAGTTGTCAAAAAATCCGAGGGGGGAACGGAAGAATGCCGGCCCTCCCCGGTTCCATGGGGGGGGACGCCCCCTGCCGGTCCTATCTCATTTGATTGGGGGGGAAATTGCGCCGGGCCAAGGGTCCATTGGCCACCAAGAATTCTCCCAGGGGGGTTCGATGGAAATTCACGGCCACCAAGGGGGGAAAGCCCCTGGGGAGACCAAATTAAACCGGTCCTCCGGTTTAATTCGATCGTTTTTTCGGGGGCTTAAAAAGGAATCCTCCCGAAGGTCACCTCGGAACCCTGGTTAG
    Z78509.1
    TCCTAAAATTAACGGGTGCCTTAACTTGCCGGGGTTAACCAAATGCCCGTCACCCCCCATTGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCTGGCCCATCAAAAACGACAACGTAGGGCATTCGCACGACAGCTCTGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTTGGACCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACTCTCATCCGCGCCGGCCGAACAGTATGCCCGTTCAACCTTAGCAATGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAACCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCTGCTTATCACATTTAAACCTTGTTACGTCCGTTGCCGAGAGTTGTCAGAAAAATTCATGGGGGGGCACGGAGCATGCCGGCCCTCCCCGCTCCATGGAGGACGACGCCCTCCGCCGGTTCTATCTCATATGACTTGGCGCGAAACTGCGCCGGGCCAAGGGTTCAATTGCCATCAAGAAATCTCCCAAGGAGGCTCGATGGCAAATCACGGTCACCACAGCGGCGAAAGCCCCTGGGAGACCAAACTACACCGGTCCTCCGGATTAACTCGATCGTATATTTGGCGGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78508.1
    TCAGCGGTGGCTCACTGACTGGGTTGCATCCAAGTGGCCGTCACCGCCCATGGGGTTGACGTGCCTCCAATGGGGTTCAAAGGGATTTATTCCGGCCCATCAAATACGACAACGTAGGGCATTCGCACGACAGCTCTGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTTGGACCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACTCTCATCCGCGCCGGCCGAACAGTATGCCCGTTCAACCTTAGCAATGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCCGAGAGTTGTCAAAAAATTCATTGGGGGCACGGCAGCATGCCGGCCCTCCCCGGCTCCATGGAGGACGACGCCCCCTGCCGGTTCTATCTCATATGACTTGGCGCGAAACTGCGCCGGGCCAAGGGTTCAATTGCCATCAAGAAATCTCCCAAGGAGGCTCGATGGCAAATCACGGTCACCACAGCGGCGAAAGCCCCTGGGAGACCAAACTACACCGGTCCTCCGGATTAACTCGATCGTATATTTGGCGGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78507.1
    TCAGCGGTGGCTCACTGACTGGGTTGCATCCAAATGGCCGTCGACCGCCACGGGTTGACGTGCCTCCAATGGGGTTCAAAGGGATTTATTCCGGCCCATCAAATACGACAACGCAGGGCATTCGCACGACAGCTTCGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTCCGGCCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACACTCATCCGCGCCGGCCGAACAGTATGCCTGTTCAGCCTTAGCAATGGGAGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTCGCTGTGTTCTTCATCGATGCAAGAGCGAGATATCCGTTGCGAGAGTTGTCAAAAAAATTCATTGGGGGACGGAAGCACGCCGGCCCTCCCCGGTTCCATGGAGGACGATGCCCTCCGCCGGTTCCATCTCATTTGACTTGGGGCGAAACTGCGCCGGGCCAAGGGTTCAATTGCCATCAAGAAATCTCCCAAGGAGGCTCGACGGAAATTCACGGTCACCACAGTAGCCAAAGCCCCTGGGAGACCAAACTACACCGGTCCTCCGGATTAACTCGATCGTTTTTTTGGGGGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78506.1
    TCAGCGGGTGGCCTCACCTGACTGGGTTGGCATCCAATGGCCGTCAAACCGCCCATGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCTGGCCCATCAAATACGACAACGTAGGGCATTCGCACGACAGCTTTGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTTTTACCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACTCTCATCCGCGCCGGCCGAACAGTATATGCCTGTTCGACCTTAGCAACGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTATGAGCTGAGGGCTCCTTGAACGAGAGTTGTCAGCCCGATTCAAAGGGGGTACGGCGGCATGCCGGTCCTCCCCGGTTCCATGGAGGACGAAGCCCTCTGCCGGTTCTATCTCATACGACTTGGCGCGAAACTGCGCCGGGCCAAGGGTTCAATTGCCATCAAGAAATCTCCCGAGGAGGCTCGATGGAAAATCATGGTCACCGCAGTAGCCAACGCCCCTGGGAGACCAAACTACACCAGTCCTCCGGATTAACTCGATCGTATATTTTGCGGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78505.1
    AAAACTCAGCGGGTGGCTCACCTGACCTGGGGTTGCATCCAAATGGCCGTCAACCGCCCATGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCCGGCCCATCAAATACGACAACGCAGGGCATTCGCACGACAGCTTTGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTTTTACCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACTCTCATCCGCGCCGACCGAACAGTATGCCTGTTCGACCTTAGCAACGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTGACGAGAGTTGTCACAAAGATTCATTGGGGGTATGGTGGTATGCCGGGCCTCCCCGGTTCCATGGAGGACGAAGACCTCTGCCGGTTCTATCTCATACGACTTGGCGCGAAACTGCGCCGGGCCACGGGTTCAATTGCCATCAAGAAATCTCCCGAGGAGGCTCGATGGAATATCATGGTCACCGCAGTAGCCACCGCCCCTGGGAGACCAAACTACACCAGTCCTCCGGATTAACTCGATCGTATATTTGGCGGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78504.1
    TTACTCAGCGGTGGCTCACTGACTGGGGTTGCATCCAATGGCCGTCACCGCCCATGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCTGGCCCATCAAATACGACAACGTAGGGCATTCGCACGACAGTTTTGTCGCAGCCACCGTCCGTCCACCTCTTGTCGGATTTGGTACCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACTCTCATCCGCGCCGGCCGAACAGTATGTCTGTTCGACCTTAGCAACGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTCGCTGCGTTCTTCATCGATGAAGAGCGAGATATCCGTTGCGAGAGTTGTCAAAAAGATTCATTGGGGGCACGGCGGGCATGCCGGCCCTCCCCGGCTCCATGGAGGGACGATGCCCTCTGACGGTTCTATCTCATACGACTTGGCGCGAAACTGCGCCGGGCCAAGGGTTCAATTGCCATCAAGAAATCTCCCGAGGAGGCTCGATGGAAAATCATGGTCACCGCAGTAGCCAACGCCCCTGGGAGACAAACTACACCAGTCCTCCGGATTAACTCGATCGTATATTTTGCGGTCTCAACAATGATCCTTCCGAAGGTTCACCTACGGAAACCTTGTTACG
    Z78503.1
    TTATCTCAGCGGGTGGCTCACCTGACCTGGGGTTGCATCCAAATGGCCGTCAACCGCCCATGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCTGGCCCATCAAATACGACAACGTAGGGCATTCGCACGACAGCTTTGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTTGGACCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACTCTCGTCCGCGCCGGCCGAACAGTATGCCCGTTCAACCTTAGCAATGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACAAAGAGCTGAGACTCCGTTGCCGAGAGTTGTCAAAATAATTCATTGGGGGTACGGTAGCATGCCGGCCCTCCCCGGTTCCATGGAGGACGACGCCCCCTGCCGGTTCTATCTCATATGACTTGGCGCGAAACTGCGCCGGGCCATGGGTTCAATTGCCATCAAGAAATCTCCCAAGGAGGCTCGATGGTAAATCACGGACACCACAGCGGCGAAAGCCCCTGGGAGACCAAACTACACCGGTCCTCCGGATTAACTCGATCGTATATTTGGCGGTCTCAACAAGGATCCTTCCGGAGGTTCACCTACGGAAACCTGGTTACG
    Z78502.1
    GCGTAAACTCACGGGGTGGCCTCACCTGACCTGGGGTTGCATCCAAATGGCCGTCAACCGCCCATGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCTGGCCCATCAAATACGACAACGTAGGGCATTCGCACGACAGCTTTGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTTGGACCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACTCTCGTCCGCGCCGGCCGAACAGTACGCCCGTTCAACCTTAGCAATGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTTCGCTGCGTCCTTCCATCGATGCAAGAGCCGAGATTCCGGCCGAGAGTTGTCAATATAAAATTCATTGGGGGGACGGTAGTATGCCGGCCCTCCCCGGCTCCATGGAGGACGACGACCCCTGCCGGTTCTATCTCATATGACTTGGCGCGAAACTGCGCCGGGTCATGGGTTCAATTGCCATCAAGAAATCTCCCAAGGAGGCTCGATGGCAAATCACGGTCACCACAGGGGCGAAAGCCCCTGGGAGACCAAACTACACCGGTCCTCCGGATTAACTCGATCGTATATTTGGCGGTCTCAACAATGATCCTTCCGGAGGTTCACCTACGGAAACCTGGTTACG
    Z78501.1
    TCTTAAACTCAGCGGGTGGCCTCACCTGACCTGGGGTTGCATCCAAATGGCCGTCGACCGCCCACGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCCGGCCCATCAAATACGACAGCGTAGGGCATTCGCTCGACAGCTTCGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTCCGGCCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACACTCATCCGCGCCGGCCGAACAGTATGCCTGTTCAGCCTTAGCAATGGGAGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAAGACTCGATGGTTCACGGGATTCTGCAATTCACAACACTTATCACAATTTCGCTGCGTCCTTCCATCCGATGCAAGGAGCCGAGATTTCCCGTTTGGCCGAGAGTTGCCCAAAAAAAAATTCATTGGGGGCACGGCAACACGCCGGCCCTCCCCGGCTCCATGGAGGACGATTCCCTCCGCCGGTTCCATCTCATATGACTTGGCGCGAAACTGCGCCGGGCCAAGGGTTCAATTTCCATCAAGAAATCTCCCAAGGAGGCTCGACGGCAAAGTCACGGTCACCACAGTAACCAAAGCCCCTGGGAGACCAAACTACACCGGTCCTCCGGATTAACTCGATCGTATATTTTGCGGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78500.1
    CTTAAACTCAGCGGGTGGCCTCACCTGACCTGGGGTTGCATCCAAATGGCCGTCGACCGCCCACGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCCGGCCCATCAAATACGACAGCGTAGGGCATTCGCACGACAGCTTCGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTCCGGCCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACACTCATCCGCGCCGGCCGAACAGTATGCCTGTTCAGCCTTAGCAATGGGAGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACAATTCGCTGCGTCCTTCCATCCGATGCAAAGAGCCGAGATTCCCGTTGGCCGAGAAGTTTTCCCAAAGAAAAATTCATTGGGGGCACGGCAGGACGCCGGCCCTCCCCGGCTCCATGGAGGGCGATGCCCTCCGCCGGTTCCATCTCATATGACTTGGCGCGAAACTGCGCCGGGCCAAGGGTTCAATTGCCATCAAGAAATCTCCCAAGGAGGCTCGACGGCAAAGTCACGGTCACCACAGTAACCAAAGCCCCTGGGAGACCAAACTACACCGGTCCTCCGGATTAACTCGATCAATTATTTTGCGGTCTCAACAATGAGCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78499.1
    GGTTGCTCACCTGACCTGGGGTCGCATCCAAAAGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCACAAGGGTCTCTAGATTATGATGGCCCATCTAATACGACAACCTGGAGCATTCGCCCAACAGTTTTGTTGTAGCATGATTCGTCCACCTCTTGCCGAATTTAGGACCATCAAAAGCCCCTGCAGCTCTTAGACCCCCAGCACCAAAGAACAAGGGGGCAAACTCACATCCGCACCGGCTGAACAGTATGCATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCTTGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATCTTGCAATTCACACACTTATCGCATTTCGCTGCGTCCTCCATCCGATGCAAAGAGCCGAGATATCCGTTGGCTGAGAGTTGTCAAAAAAATAATTGGGAGAGGTCAAGGCCACATGCCGTCCCCTCCTCGAATTCAACGAAGGGGGTATGTCCTTCACCATTTATGTGTCATATGACTTGGTGCAAAACTGCGACGGGCAAGGGTTAGATCGCCGGCAAGAAAGCTCCCAAGGAGGCTCCCTCCATGGAAAATCTAGGTCACCGCAATAGCAGGCGCCCAAGGGTGACCAAAGTAAACCGATCCTCTGGATTATCTCGATCAATTATTATGCGATCTCAACAATGATCCCTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78498.1
    GCTTAAACTCAGCGGGTTGCTCACTGACTGGGGTCGCAACAAATGGCCATCAACTGATCACGGGTTGATGTGCCTCCAGTGGGGTTCACAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCTGGAGCATTCGCACAACAGTTTTGTTGTAATTCGTCCACCTCTTGCCGAATTTAGGACCATCAAAAGCCCCTGCAGCTCTTAGACCCCCAGCACCGAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGCATGTACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTTGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATCCTGCAATTCACACCACTTATCGCATTTCGCTGGGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTCCAAAAATAATTTGGAGAGGGGGAGTGATTGTAGGGTCAAGGGCAAAATGCCGCCCCCTCCTTCGCCATTATTGTGTCATATGACTTGGGGCAAAACTGCCCCGGGCAAGGGTAAGATCGCCGGCAAGAAAACTCCCAAGGAGGCTCCATGGCAAATCTAGGTCACCGCAATAGCAAGCACCCATGGGTGACCGGAGTAAACTGATCCTCTGGAGTAACTCGATCAATTATTATGTGATCTCAACAATGACCTTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78497.1
    GCTTAAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCACAAGGGTGTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCGCGTTGCGTCCACCTCTTGCCGAATTTAGGACCATCAAAAGCCCCTGCAGCTCTTAGACCCCCAGCACCAAAGAACAAGGGGCCAAACTCACACCCGCACCGGCTGAACAGTATGCCTGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTTGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCGCATTTCGCTGCGTTCCTTCATCGATGCAAGAGCTGAGATCTCCGTTGCTGAGAGTTGTTTCACAAAAAATAATTTGGAGAGGGGGAGGGAGTGTAGGCTCAAGGCCACATGCCCTTCACCTCCTTGAATTCACAAAGGGCCGTGCCCTTCACCATTTATGTGTCATGTGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTAGATCGCCGGCAAGAAAGCTCCCAAGGAGGCTCCATGGCAAATCTAGGTCACCGCAATAGCAAGCGCCCATGGGCGACCAAAGTAAACTGATCCTCTGGATTCACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78496.1
    GCTTAAACTCAGCTGGTTGCCTCACCTGACCTGGGGTCGCAACCCAAAAGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCACAAGGGTCTCTAGATTATGATGGCCCATCTAATACGACAACCTGGAGCATTCGCCCAACAGTTTTGTTGTAGCATGATTCGTCCACCTCTTGCCGAATTTAGGACCATCAAAAGCCCCTGCAGCTCTTAGACCCCCAGCACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGCATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTTGGGCGCAACTTGCGTTCAAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCGCATTTCGGTGCGTCCCTCCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTCAAAAAATTAATTGGAGAGGTCAAGGCCACATGCCGGCCCCTCCTCGAATTCACGAAGGGATATGCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTAGATCGCCGGCAAGAAAGCTCCCAAGGAGGCTCCCTCCATGGCAAATCTAGGTCACCGCAATAGCAGGCGCCCAAGGGTGACCAAAGTAAACCGATCCTCTGGATTAACTCGATCAATTATTATGCGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78495.1
    CACTCAGCGGGTTGCTCACTGACCTGGGGTCGCAACCGAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCACAAGGGTCTCTGGATTATGCCGGCCCATCTAATACGACAACCTGGAGCATTCGCACAACAACAGTGTTGTTGTAGCCTGATTCGTCCACCTCTTGCCGAATTTAGGACCATCAAAAGCCCCTGCTGCTCTTAGACCCCCAGCACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGCATGGACAGCCTCATTGAGGGAGAGAGATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTCATGGCCTCGGGCGCAACTTGCGTTCAAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCGGATTTCGCTGCGTTCCTCCATCCGATGCAAGAGCTGAGATATCCGTTGCTGAGAGTTGTTCAGAAAAATACTTTGGAGAGGGGGAGAGCGAGTGCAGTGTCAAGGCCACATGCCGCCCCCCCTCCTTGAATACACGAAGGGCTATGCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACCGCGCCGGACAAGGGTTAGATCGCCGGCAAGAAAGCTCCCAAGGAGGCTCTCCATGGCAACTCTAGGTCACCGCAATAGCAAGCGCCCATGGGTGACCAAAGTAAACCGATCCTCTGGATTATCTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGGAGGTTCACCTACGGAAACCTTGTTACG
    Z78494.1
    CTTAAACTCAGCGGGTTGCCTCACCTGACCTGGGATCGCAACCAAATGGCCATTCAACTGATCATGGGTCGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCTGATTAAACACAACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCGAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTACGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGACAGGCGTTCCCTTGGCCTGTTGGCCTCGGGCCCAACTTGCGTTCAAAGACTCGATGTTAACGGGATTCTGCAATTCACACCTGAGATATCCGAAGTTGAGAGTGTTACCTCAATAATGGGGAGTTGGTCAAGGCAGTATGCCGTCCCCCTTCATTAATTATGGGTCATTTGATTTGGCGCATTACTGCGCCGGGCAAGGGTATAGATCGCCAGCAGGAAAGTTCCCAAGGAGGCCCGATGGTAAATCTAGGTCACTGCAACAGTAAATGTCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGCGACCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78493.1
    GGGTTGCCTCACCTGACCTGGGATCGCAACCAAATGGCCATTCAACTGATCATGGGTCGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCTGATTAAACACAACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCGAAGAACAAGGGGCCAAACTCACATCCACACCGGATGACCAGTACGTTTGGACAGCCTCATTAAGGGAGAGATATGAGTCGCAATGTCCAGGCAGGCGTCCCCTTGGGCTGATGGCCTCGCGCGCAACTTGCGTTCAAAGATTCGATGGTTCACGGGATTCTGCAAGGCACACCATTTATCGAATCTCGCTGCGTTCCTTCATCGATGCATGAGCCGAGATATCCGTTGCTGAGAGTTGTTACCTAAATACTTGGGGGCTGGTCAAGGAAGAATGCCGCCCCCCTTCACCACTTATGTAAAATTTGACTTGGCGCAAAACTGCGCCGGGCAAGGGTATAGATCGCCAGCAGGAAAGCTCCCAGGGAGGCCCGATGGAAATTCTAGGTCACTGCAACAGAAAATGCCCATGGGTGACCAAAGTAAACCGATCCTCCAGATTAACTCGATCAATTATTATGCGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAATCCTTGTTACG
    Z78492.1
    TATTAAACTCAGCGGTGTGCCTCACCTGACTGGGATCGCAACCAAATGGCCAATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTCCAAAAGGGTCTTCTGATTAAGCTGGCCCATGTAACACAACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAATTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCGTGCAGCTCTTAGACCCCCCGTCCCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATTTGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTCCACGGGATTCTGCAAATCACACCAGTATCGCATTTCGCTGGGTTCTACATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGAGTGGGTCAAGGCAGTACGCTCCCCCCCTTCACCAATTATGTGACATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTATAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATTCTAGGTCACTGCAACAGCAAATGCCCGTGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGCGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78491.1
    GCTTATCTCAGCGGGTTGCTCACCTGACTGGGATCGCAACAAATGGCCATTCAACTGATCATGGGTTGATGGGCCTCTAATGGGGTCCAAAGGGTCTTCTGATTATGCTGGTCCATGTAACACAACAACCCGGGGCATTCGCACAACATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAATATAATTTGGAGCGGGTCAAGGCAGCACGCCTCCCCCCAACACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTATAGATCGCCAGCAAGAAAGCTCCTAAGGAGGCCCGATGGCAAATCTAGGTCACTGCAACAGCAAATGCCCGTGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGCGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78490.1
    TCAGCTGGTTGCTCACCTGACTGGGGTCGCAACCAAATGGCCGTCAACTGATCATGGGTTGACGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAAATTATGCTGGCCCATCTAATACGACAACGCGGGGCATTCGCACAACACTTGTTGCAGAACAGTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCGTTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTTTGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAATAATTTGGGGAGGGTCGAGGCAGCATGCCGCCCCCTTCACCAATTATGTGCCATATGACTTGGCGCAAAACTGCGCCGGACTAGGGTTCAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTCGGTCACTGCAACAGCAAATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78489.1
    GCCTAAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCATGGGGTTCAAAAGGGTCTTCAGATATGCTGGGCCCATCTAATACGACAACTCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAATCCCATGCAGCTCTTAGACCCCCCGTACCAATGAACAAGGGGCCAAACCCACATCCCCACCGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAATAATCTGGGGAGGGTCAAGGCAGCATGCCGCCCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCTCGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTCGGTCACTGCAATAGCAAATGCCCTCGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78488.1
    AGCCGGTTGCTCACCTGACTGGGGTCGCAACAAATGTCCATCAACTGATCATGGGTTGATGGGCCTCCGATGGGGTTCAAAAGGGTCTTCAGATTATGATGGCCCATCTAATACGACAACCCGGGGGATCCGCACAACAGTTTGGTTGTGGCGTTGTTCGTCCACCTCTTGCCGTTTTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAATGAACAAGGGGGCAAACTCACATCCGCACCGGTTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATTTGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGGGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTCCATCCCGTTGATGAGAGTTGTTAGAAAATAATTAGGGGAGGGTCAAGGCACCATGCCGCCCCCTTCACCAATTATGTGTCGTATGACTTGGCGCAAAACTGCGCCGGGCTAGGGTTAAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTCGGTCACTGCAATAACAAATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTGCGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACAG
    Z78487.1
    TTACTCAGCGGTTGCTCACCTGACTGGGGTCGCAACAAGTGGCCCCAACTGATCATGGGTTGATGGGCCTCTAATGGGGTTCAAAAGGGTCTTTAGATTACGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGCATTTAGGACCATCGATAGCCCATGCTGCTCTTAGACCCCCCGTACCGAAGAACAAGGGGCCGAACTCACATCCGCACCGGCTGAACAATATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAATAATTTGGGGAGGGTAAAGGCAGCATGCCGCCCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCTAGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCCCGGTCACTGCAATAGCAAATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGAGTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78486.1
    TCAGCGGGTTGCCTCACCTGACCTGGGGTCGCTACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCACGCAGCTCTTAGACCCCCCGTACCAATGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCGAGATATCCGTTGCTGAGAGTTGTTAAGAAATAATTTGGGGAGGGTCAGGCACCATGCCGCACCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCTAGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTCGGGCACTGCAATAGCAAATGTCCATGGGTGACCAAAGTAAACTGATCCTCCAGATAATCTCGATCAATTATTATGTGATCTCAACAATGATCCTCCCGCAGATTCACCTACGGAAACCTCGTGACG
    Z78485.1
    TACTAAATCAGGGGTTGCTCAGCTGACTGGGGTCGCAACACATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAGGGGTCTTTAGATTGTGCTGGCCCGTCTAATACGACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGGCGTATTTAGGACCATCAAAAGCCCACGCAGCTCTTAGACCCCCCGTACCAATGAACAAGGGCCCAAACTCACATCCGCACCGGCTGCACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGTTCACGGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAGAAATAATTTGGGGAGGGTCAAGGCACCATGCCGCACCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCCCCGGGCTAGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTCGGTCACTGCAATAGCAAATGCCCATGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGATATCTCAACAATGATCCATCCGCAGATTCACCTTCGGACACCAGGTTCAG
    Z78484.1
    AAAACTCAGGGAGTTGCTTCACCTGACTTGGGGTCGCAACCCAATGGACCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTGTCAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCACGCAGCTCTTAGACCCCCCGTACCAATGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCCAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCGAGATATCCGTTGCTGAGAGTTGGGAGGGTAGGAAACATGCCGCCCCTTCACAATTATGTGTCATATGACTTGGGCCAAAACTGCGCCGGGCTAGGGTTTAGATCGCCAGCAAGAAAGCTCCCATGGAGGCTCGATGGTAAACTCTCGGTCACTGCAATAGCCAAATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCCCAGGTTCACCTACGGAAACCTTGTTACG
    Z78483.1
    TGCTAAACTCAGGGGTTGCTTCACTTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCAACAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAATGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGATATGTGTGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATTGTTCACGGGATTCTGCAATTCACACCATTTATGATTTCGGTTTGGTCTTCATCCATGCAAGAGGCCAGATTTCCGTTGCTGAGAGTTGTTAAAAAATAATTTGGGGGAGGGTCAAGGCAGCATGCCGCCCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCTAGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTCGGTCACTGCAATAGCAAATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78482.1
    TGCTAAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCACGCAGCTCTTAGACCCCCCGTACCAATGAACAAGGGGCCAAACGCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGNTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGAAGAGGCGAGATATCCGTTGCNGAGAGTTGTTAAAAAATAATTTGGGGAGGGTCAAGGCAGCATGCCGCCCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCTACGGTTTAGATCGCCAGCAAGAAAGCTCCCAGGAGGCTCGATGGCAAATCTCGGTCACTGCAGTAGA
    Z78481.1
    TCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCAACAAAAGCCCACGCAGCTCTTAGACCCCCCGTACCAATGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTGTGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTCCGTTGCTGAGAGTTGTTAAAAAAATGATTTGGGGAGGGTCAAGGCACCATGCCGCCCCCTTCGCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCTAGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTCGGTCACTGCAATAGCAAATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78480.1
    TCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTATACGACAACCCGGGGCATTCGCACAACATTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGTACCAACAAAAGCCCACGCAGCTCTTAGACCCCCCGTACCAATGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTGTGGACAGCCTCATTGAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTATCCGTTGCTGAGAGTTGTTAAAAAATAATTTGGGGAGGGTCAAGGCAGCATGCCGCCCCCTTCGCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCTAGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTCGGTCACTGCAATAGCAAATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78479.1
    ACTCAGAGGGTTGCCTCACCTGACCTGGGGTCGCATCCAAATGGCCGTCAACTGATCTTGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAGTTCTGTTGTCGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACGACAGCCTCATTAAGGGAGAGATATGTCTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAGACTCGATGGTTCACGGGATCTCTGTATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGGCGAGATATCCGTTGCTGAGAGTTGTTAAAAAATAATTTGGGGGGCAAGGCAAGGCAGCATCCCTCCCCTTCCAGTAATGTCTCATATAAATTGGCGCAAATCTGCGCCGGGCAAGGGTTTAGTTCGGCTCGATGGCAAATCTAGGTCACTGAAACAGCAAATGCCCATGGGTGACCAAAGTAACCTGATCCTCCTGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78478.1
    GCCTAAACTCAGCGGGTGCCTCATCTGACCTGGGTCGCAACCAAATGTCCTCAACTGATAATGGGTTGATGGGCCTCCAATGGGGGTCAAAAGGGTCTTTAGATTATGCTGACCCATCTAAACGACAACCGGGGCATTCGCACAACAGTTCTGTTGTCACATAGCTCGTCCACCTCTTGCCGTATTTAGGACCATCCAACCCATACAGCTCTTAGACCCCCGTACCAAAGGACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGCAGCCTCATTAAGGGAGAGATATGTCTCACAATTTCCAGGCAGGCGTCCCCTTGACCTGATGGCCTCGGACGCAACATGCGTTCAAAACTCGATGTTCAGGGATCTGCTATCACGCCATTATCAATATTTGGGGGAGCCATGGCAGCATGCCGCCCCTTCCAGTTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCGGCTCGATGGCAAACTCTAGGCCACTGCAACAGCAAATGCCCATGGGTGACCAAAGTAAACTGATCCCCCAGATTAACTCGATCAATTATTATGTGATCTCAACACTGATCCTTCCGGAGGTTCACCTACGGAAACCTTGTTACG
    Z78477.1
    GCATAAACTCAGCGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAGTTCTGTTGTCGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTGTGGACAGCCTCATTAAGGGAGAGATATGTCTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGCCGCAACTTGCGTTCCAAAGACTCGATGGTTCAGGGATTCTTGCAATTCACACCATTTATCGCATTTCGCTGCGTCCTTCCATCCGATGCAAGAGCCGAGATACTCGTTGCTGAGAGTTGGTTAACAAATATTTGGGGAGGGCAAGGCAGCATGCCGCCCCTTCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCGGCTCGATGGCAAATCTAGGTCACTGCAACAGCAAATGCCCATGGGTGACCAAAGCAAACTGATCCTCCAGATTAACTCGATCAATTATCATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78476.1
    GGGTCGCAACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGATTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCGCAACAGTTCTGTTGTCGCATAGTTCGTCCACCTCTTGCCGTATTTAGGTGGATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACTGTATGTATGGACAGCCTCATTAAGGGAGAGATATGTCTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCGGATGGCCTAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGCCAAGGCAGCATGCCGCCCTTCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCGGCTCGATGGCAAATCTAGGTCACTACAACAGCAAATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATCATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78475.1
    ACCTGACCTGGGTCGCAACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTTGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCATTGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNAAGAGCCGAGATATCCGTTGCTGAGAGTTGTCAAAAAAATAATTTGGGGAGGGTCAAGGCAGCATGCCACCCCCTTCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGACGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCCAGGTCACTGCAAGAGCAGATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78474.1
    AAGCTCAACGGGTTGCTTCACCTGACCTGGGGTCGCAACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAGAAGGGTCTGTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTTGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGTGTCTTCATCGATGCAAGAGGCGAGATATCCGTTGCTGAGAGTTGTCAAAAATATAATTTGGGGAGGGTCAAGGCAGGATGCCGCCCTTCCAATTATGTGTCATATGATTGGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCCAGGTCACTGAAAGAGCAGATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTACGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78473.1
    CCTGGGTCGCAACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGTACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTAGGCAAGAACAAGGGGCCAAACTCACATCCGCACCGACTGAACAGTATGTATGGACAGCCTCATTGAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAACAAAATAATTTGGGGAGGGTCAAGGCAGCATGCCGCCCCCTTCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCGCCACCAAGAAAGCTCCCATGGAGGCTCGATGGCAAATCCAGGTCACTACAAGAGCAGATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACCCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78472.1
    GCTTAAACTCAGCGGTTGCCTCACCTGACTTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATTGGCCTCTAAAGGGGTTCAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCTGGGCATTCGAACTACAATTTTGTTGTAGCATACTTCGTCCACCTCTTGCCGTATTTAAGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCTGTGGCTGAACAGTATGTATAGACAGCCTCATTAAAGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTCCTCCATCGATGCAAGAGCCGAGATTCCGTTGCTGAGAGTTATCAAAAAATAATTTGGGGAGGGTCTTAGACTGCATGCCGCCCCCTTCCAATTATGTGTGAAATGACTTGGCGCAAGACTGCGCCGGGCAACGATTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCCAGGTCACTGCAAGAGCAGATGCCCATGGGTGACTAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78471.1
    GCTTAAACTCAGCGGGTTGCCTCACCTGACTTGGGGTCGCAACCAAATGGCCATCAACTGATCAGGGTTGATGGGCCTCTAATGGGGTTCAAAAGGGTCTTTAAATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACTACAATTTTGTTGTAGCATACTTTGTCCACCTCTTGCCGTATTTAGGACCAGCAAAAGCCCATGCAGCTTTTAGACCCCCCGTACTAAAGAACAAGGGGCCAAACTCACATCCGCACTGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAAGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGATGTTGAGAGTTGTCAAAAAATTAATTTGGGGAAGGGTCAAGGCAGCATGCCGCCCCCTTCCATTATGTGTCGGGGTGACTTGACGCAAGGCTGCGTCGGGCAACGATTTAGATCGCCAGCAAGACAGCTCCCAAGGAGGCTCGATGGCAAATCCAGGTCACTGTAAGAGCAGATGCCCATGGGTGACCAAAGTAAACCGATCCTCTAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78470.1
    AACTGATCAGGGTTGATGGGCCTCTAATGGGGTTCAAAAGGGTCTTTAAATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACTACAATTTTGTTGTAGCATACTTTGTCCACCTCTTGCCGTATTTAGGACCAGCAAAAGCCCATGCAGCTTTTAGACCCCCCGTAGTAAAGAACAAGGGGCCAAACTCACATCCGCACTGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAAGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTAAAAAAATAATTTGGGAAGGGTCAAGGCAGCATGCCGCCCCCTTCCAATTATGTGTCATATGACTTGACGCAAGGCTGCGTCGGGCAACGATTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCCAGGTCACTGTAAGAGCAGATGCCCATGGGTGACCAAAGTAAACTGATCCTCTAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78469.1
    AACTGATCATGGGTTGATGGGCCTCTAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCTGGGCATTCGCACTATCATTTTGTTGTAGCATACTTCGTCCACCTCTTGCCGTATTTAGGTACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCTGTGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGAAGAGCCGAGATATCCGTTGCTGAGAGTTGTCAAAAAAATAATTTGGGGAGGGTCTAGACCGCATGCCGCCCCCTTCCAATTATGTGTGATATGACTTGGCGCAAGACTGCGCCGGGCAACGAATTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCCAGGTCACTGCAAGAGCAGATGCCCATGGGTGACTAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78468.1
    AACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTATACGACAACCCGGGGCATTCACACAATAATTTTGTTGTAGCATACTTCGTCCACCTCTTGCCGTATTTAGGTCCATCAAAAGCCCAGCTGCTCTTAGACCCCCCGTACCGAAGAACAAGGGGCCAAACTCACATCCGCACTGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGAGTTCCGGAGCTGAGAGTCGTTAAAAAAATGATATGGGGAGGGTCAAGGCAGCATGCTGCCCCTCCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTCAGATCGCCATCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCCAGGTCACTACAAGAGCAGATTCCCATGAGTGACCAAAGTAAACTGATCCTCCAGATTAACTCAATCAATTATTATGCGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78467.1
    TCAGCGGGTTGCCTCACCTGACCTGGGTCGCAAACATATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAGAGGGTCTTTAGATTATGCTGGCCCAGCTAATACGACAACCCGGGTCATTCGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTTGGACCACCAAAGGCACATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTACCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTCAAAAAATAATTTGGGGAGGGTCAAGGCAGCATGCCACCCCCTTCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCCAGGTCACTGCAAGAGCAGATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78466.1
    GGGTTGCCTCACCTGACCTGGGTCGCAACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACGCAAGAACAAGGGGCCAAACTCACATCCGCACAGGCTGAACTGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGCAAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTACAAAAATAATTTGGGGAGGGTCAAGGCAGCATGCCGCCCCCTTCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCGCCAGCAAGAAAGCTCCCATGAAGGCTCGATGGCAAATCCAGGTCACTACAAGAGCAGATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78465.1
    GCAACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACTCGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTTCCCTTGGCCTGATGGCCTCGGGCGCAACTTTCGTTCAAAGACTCGATGGTTTACGGGATTCTTCAATTCACACCAATTATCGCATTTCGCTGCGTTCCTCATCGATTCAAGAACCGAGAAATCCGTTGCTGAGAGTTGTCAAAAAAATAATTTGGGGAGGGTCAAGGCAGCATGCCGCCCCCTTCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCCAGGTCACTGCAAGAGCAGATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78464.1
    GCTTAAACTCAGCGGGTTGCTCACCTGACCTGGGGTCGCACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTATTACGAAAGCCCGGGGAATCCGAACAGAAATTTGGTTGAAGAATAGTCCGCCCACCTCTTGCCGTTTTTAGGACCACCAAAAGCCCATGAAGTTCTTAGACCCCCCGCCCCAAAGAAAAAGGGGCCAAACTCACATCCGAACCGGTTGAACAGTATGTTTGGACAGCCTCATTAAGGGAGAGATTTGATTCGAAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGAACTTGGCGTTCAAAGACTCGATGGTTCACGGGATTCTGAAATTCACACCATTTATCGGATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAGAAAATAATTTGGGGAGGGTCAAGGCAGCATGCCGCCCCCTTCCAATTATGTGTCATATGACTTGGCGCAAAACTGCTCCGGGCAAGGGTTTAGATCTCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAATCCAGGTCACTACAAGAGCAGATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGACTCAACCGATCAATTATTATGTGATCTCAACAATGACCCTTCCGCTCACCTACGGAAACCTTGTTACG
    Z78463.1
    GCTTAAACTCAGCGGGTTGCTCACCTGATCTGGGGTCGCAACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGTATTCGCACAGCAATTTTGTTGCAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGCCCCAAAGAAAAAGGGGCCAAACTCACATCCGCACCGGTTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAGAAAATAATTTGGGGAGGGTCAAGGCAACATGCCGCCCCCTCCCAATTATGTGTCATATGACTTGGCGCAAAACTCCCCCGGGCGAGGGTTTAGATCCCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAATCCAGGTCACTACAAGAGCAGATGCCCATGGTGACCAAAGTAAACTGATCCTCCAGACTCAACCGAACAATGATTATGTGATCTCAACAATGAACCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78462.1
    ATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAATCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCCCATGGGTGACCAAAGTAAACAGATCCTGCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGGTCACATCCGGAGACCTCGTGACG
    Z78461.1
    TTAACTCAGCGGCTTGCTCACTGACTGGGTCGCAACAATGGCATCAACTGATCATGGGTTGATGGGCTCCAATGGGGTTCAAAGGGTCTTCAGATTACGGTGGCCCATCTTATACGACAACCCGGGGGCTTTGGACAACAATTTTGGTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGGTCTTAGACCCCCCGTACCAAAGAACAAGGGGGCAAACTCACATCCGCACCGGCTGAACAGGATGTATGGACAGGCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTTCCCTTGGCCTGATGGCCTCGGGCGCAACTTGGGTTCAAAGACTCGATGGTTCACGGGATTCTGGAATTCACACCATTTATCGGATTTCGGTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTGAAAAAATTATTTGGGGGAGGGTCAGGGAAGAATGCCACCCCCTCCACCAATTATGTGTCATTTGACTGGGGGAAAAATTGCGCCGGGTAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGAAATTCTAGGTCACTCCAGCAGCAACTGCCCATGGGTGACCAAAGTAAACAGATCCTCCAGATTACCTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGGAGGTTCACCTACGGAAACCTTGTTACG
    Z78460.1
    TAAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTGGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCCTCCATCCGATGCAAGAGCCGAGATATCCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGTCAGGGCAGGATGCCACCCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGGAGGTTCACCTACGGAAACCTTGTTACG
    Z78459.1
    AAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCAATTTATCGCAATTTCGCTGCGTTCCTCCATCCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGTCAGGGCAGGATGCCACCCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78458.1
    CAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTCTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATCTGCAATTCACACCATTTATCGCATTTCGCTGCGTCCTCCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAATAATTTGGGGAGGGTCAGGGAAGGATGCCACCCCCTTCACCAATTATGTGTCGTACGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78457.1
    CTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACGATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGGTCCTTCATCCGATGAAAGAGCCGAGATATCCGTTGTTGAGAGTTGTTAAAAAAATATTTGGGGGAGGGTCAGGGAAGAATGCCACCCCCTCCACCAATTATGTGCCATTTGACTGGGGGAAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGAAATTCTAGGTCACTACAACAGAAAATGCCCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGGAGGTTCACCTACGGAAACCTTGTTACG
    Z78456.1
    GCTAAACTCAGGGTTGCCTCACCTGACCTGGGTCGCAACCCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTCCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGTTCACGGGATTCTGCAATCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCCAGACCGGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGTCAGGGCAGGATGCCACCCCCTCCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATTCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78455.1
    TGCTAAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGAAGCTCTTAGACCCCCCGTGCCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGGTTCTTCAACCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAATTATTTGGGGGAGGGTAAGGGAAGGATGCCACCCCCTCCACCATTTTTGTGTAATTTGATTGGGGGAAAAACTGCGCCGGGAAAGGGTTTAGTTCTCGCCAAAAAGAAAGTTCCCAAGGGGGTTCGATGGAAATTCTAGGTCACTCCAAAAGAAATTGTTCATGGGTGACCAAAGTAAAAAGTTCCTCCAGATTAACTCGATCAATTTTTATGTGATCTCAAAAATGATCCTCCCGAAGGTCCACCTACGGAAACCTGGTTACG
    Z78454.1
    GTTGCCTCACCTGACCTGGGGTCGCATCCAAATGGCCATCAACTGATCATGNGGTTGATGGGCCTCCAATGGGTTCAAAAGGGGCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAGAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGTCTGATGGCCTCGGCCGCAACTTGCGTTCACAGACTCGATGGTTCACGGGATTCTGCAAATCACACCATTTATCGCATGAGTTCCGTTGATGAGAGTTGTTAACAATAGTGGGGGGAGGGTCAGGGCAGGATGCCACCCCCTTCACCGATTATGTGTCAAATGACTTGGAGAAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGTAAATCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTATCTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78453.1
    TGCTAAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCTGAGATATCCGGTTGCTGAGAGTTGTTAAAAAATAATTTGGGGGAGGGTCAGGGAAGGATGCCACCCCCTTCACCAGTTATGTGTCAAATGAGTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAACTCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78452.1
    TGCTAAACTCAGCGGTTGCCTCACCTGACCTGGGGTCGCATCCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACAGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCATCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATTCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGTCAGGGCAGGATGCCACCCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAAAAAGAAAGCTCCCAAGGAGGCTTGATGGCAAATCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAATAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78451.1
    GCTCAGCTGGTTGCTCACTGACTGGTGTCGCATCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCTATAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGGCTTCGCACAACAATTTTATTGTAGGGATAGTTCGTCCACCTCTTGCCGTTTTTAGGACCATCAAAATCCCATGCAGGTCTTAGACCCCCCGTACCAAAGAACAAGGGGGCAAACTCACATCCGCACCGGGTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATTTGACTCGCGATGCCCAGGGAGGGGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGGGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGATTTCGCTGGTTCTTCATCGATGAAAGAGCGAGATTTCCGTTGTTGAGAGTTGCAAAAAATACTTTGGGGAGGGTCAGGGCAGGATGCCCGCCCCCTACACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCACGGGTTTAGATCTCGCCAGCAAGAAAGCTCCCAGGGGGGCTCCATGGAAATCTAGGTCACTACAACAGCAAATGCCCATGGGTGACCAAAGTAAACAGATCCTCCAGATTATCTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGGAGGTACACCTACGGAAACCTTGTTACG
    Z78450.1
    CTCAGGGGTTGCTCAGCTGATCTGGGGTCGCAACACATGGTCATCAACTGATCATGGGTTGATGGTCCTCCAATGGGGTTCAACAGGGTCTACAGGTTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTGGCACAACAATTTGGTTGTAGCATAGTTCGTCCACCTCTTGCCAGTTTTGAGGACCATCATATGCCCATGCAGTTCTTAGACCCCCCGGACCAAAGAACAAGGGGCCACACTCACATCCGCACCGGATGAACAGTATGTATGGACAGCCTCGTTAGGGGAGAGATATGACTCGGAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGTCAGGGCAGGATGCCACCCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGGCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGTAAATGCTCATGGATGACCAAAGTAAACAGATCCTCCAGCTTAACTCGATCAATTATTATGTGATATCAGCAATGATCCTTCC
    Z78449.1
    GCATAAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTAGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTTCGCTGCGTTCTTCATCGATGCAAGAGCTGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGGAGGGTCAGGGCAGGATGCCACCCTTCACCAATTATGTGTCAAATGAGTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78448.1
    CCTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCANCCAAATGGCCATCAACTGATCATGGGCTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTNTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGGTGCGTTCTTCCATCCGATGCAAGAGCTGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGTCAAGGGAAGGATGCCACCCCCCTTCACCAATTATGTGTCATACGACTTGGCGCAAAACTGCGCCGGGAAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAATTCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78447.1
    GCTTAAACTCAGCGGGTTGCCTCACCTGACTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGGTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGGGTTCTTCATCGATGCAAGAGCCGAGATAGGATGCCACCCCCTTCACCACTTATGTGTCATTTGACTGGGGGCAAAATTGCGCCGGGTAAGGGTTTAGTTCTCGCCAACAAGAAAGCTCCCAGGGAGGCTCGATGGCAACTCTAGGTCACTTCAACAGGAAATGCCCATGGGTGACCAAAGTAAACAGATCCTCCAGATTACTCGATCAATTATTATGTGATCTCAACAATGATCCTCCCGCAGGTTCACCTACGGAATCCTTGTTACG
    Z78446.1
    GGGTTGCCTCACCTGACCTGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGAGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGGCCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCGAAGAACAAGGGGCCAAACTCACATCTGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGTTCACGGGATTCTGCAATCACACCATTATCGCAATTTCGCTGGTCCTCCATCGATGAAGAGCCGAGTATCCGTTGCTGAGAGTTGTTAAAAAATTATTTGGGGAGGATGCCACCCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGGAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGGAGGTTCACCTACGGAAACCTTGTTACG
    Z78445.1
    ACAGCTGTTGCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTCTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTGCACGGGATTCTGCAATTCACACCATTTATCGCATTTGCCGAAGACTGAGTCTCCGTTGCTGAGAGTTGTTAAAAAAATAGTTGGGGGAGGGTCAGGGAAGGATGCCACCCCCTTCACCAATTATGTGTCAAACGACTTGGAGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAATTCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACCCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78444.1
    AATGGCCATCAACTGACATGGGTTGATGGGCTCCAATGGGGTCCAAAAGGGTCTTCAGATATCGGTGGCCCATCTTATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCATCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGGTGGGTTCTTCATCGATGCAAGAGCTGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATATTTGGGGGAGGGTAAGGGAAGAATGCCACCCCCTCCACCAATTTTGTGTAATTTGATTGGGGGAAAAATTGAGGGTTTAGATATCGCCAACAAGGATGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCCCATGGGTGACCCAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGTTCACCCTACGGAAACCTTGTTACG
    Z78443.1
    CCTCACCTGACTGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTATGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGGTGTGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTATCGCATTTCGCTGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTGTGGCCCAAAACTCTGCCGGGCAAGGGTTTAGATATCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCCCATGGGTGACCCAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78442.1
    ACTCAGCGGGTTGCTCAGCTGACTGGGGTCGCAACAAATGGTCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGGTTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTGGGACAACAATTTTGTTGTGGCATAGTTCGTCCACCTCTTGCCGTTTTGAGGACCATCAAATGCCCATGCAGCTCTTAGACCCCCGGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGATGAACAGTATGTTTGGATAGCCTCGTTAGGGGAGAGATATGATTCCCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGTCAGGGCAGGATGCCACCCCCTTCACCAATTATGTGTCATATGACTTGGCGCCAAACTCCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCCCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTAC
    Z78441.1
    CTCAGCGCGTTGCTCAGCTGACTGGGGTCGCAACACATGGTCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGGTTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTGGGACAACAATTTGGTTGTAGCATAGTTCGTCCACCTCTTGCCGTTTTGAGGACCATCAAATGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGATGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGTCAGGGCAGGATGCCACCCCCTTCACCAATTATGTGTCATATGACTTGGCGCACAACTCCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAACTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATATCGGCAATGACCTTCC
    Z78440.1
    TGCTAAACTCAGGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGCTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTAGCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGGTGGGTTCTTCAACGATGCAAGAGCTGAGATATCCGTTGCTGAGAGTTGTTAAAAAATTATTTGGGGGAGGGTAAGGGAAGATATGCCACCCCCTCCACCATTTTTGTGTAATTTGATTGGGGGAAAAACTGCGCCGGGTAAGGGTTTAGTTCTCGCCAACAAGAATGCTCCCAAGGAGGGTCGATGGAAATTCTAGGTCACTACAACAGAAATTGCCCATGGGTGACCAAAATATGCAGATCCTCCAGATTACCTCGATCAATTATTATGTGATCTCAACAATGATCCTCCCGGAGGTCCACCTACGGAAACCTTGTTACG
    Z78439.1
    GGCCCAACTAAACGGCCAACCGGGGCATTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTTCCGTATTTAGGACCATCCAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGTTCACGGGATTCTGCAATTCACACCATTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGTCAGGGCAGGATTACCACCCCCTTCACCAATTATGTGTCATATGACTTGGCGCCCAACTCCGCCGGGCAGGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCCCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATG



```python
# We can label the reverse complements as reverse complements by annotating the label "reverse complement" into each description
records = [
    rec.reverse_complement(id = "rc_" + rec.id, description = "reverse complement")
    for rec in SeqIO.parse("ls_orchid.fasta.txt", "fasta")
]
```


```python
len(records)
```




    94




```python
# We can set conditions on when to use the label such as having a sequence less than 700 long
records = [
    rec.reverse_complement(id = "rc" + rec.id, description = "reverse complement")
    for rec in SeqIO.parse("ls_orchid.fasta.txt", "fasta")
    if len(rec) < 700
]
```


```python
# We can see the length of the records under 700 long
len(records)
```




    18




```python
# We can parse the 18 sequences under 700 long out of the file and compile them into a new file named "rev_comp"
records = [
    rec.reverse_complement(id = "rc_" + rec.id, description = "reverse complement")
    for rec in SeqIO.parse("ls_orchid.fasta.txt", "fasta")
    if len(rec) < 700
]

SeqIO.write(records, "rev_comp.fasta", "fasta")
```




    18




```python

```



# Multiple Sequence Alignment
## Part 1

```python
from Bio import AlignIO
```


```python
# We can import our file and assign it to the "alignment" variable
alignment = AlignIO.read("PF05371_seed.sth.txt", "stockholm")
```


```python
# We can print the new "alignment" variable
print(alignment)
```

    Alignment with 7 rows and 52 columns
    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRL...SKA COATB_BPIKE/30-81
    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKL...SRA Q9T0Q8_BPIKE/1-52
    DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRL...SKA COATB_BPI22/32-83
    AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA COATB_BPM13/24-72
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA COATB_BPZJ2/1-49
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA Q9T0Q9_BPFD/1-49
    FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKL...SRA COATB_BPIF1/22-73



```python
# We can format the rows how we want
alignment = AlignIO.read("PF05371_seed.sth.txt", "stockholm")
```


```python
# We can print out the specific length of our alignment and label it as such
print("Alignment length %i " % alignment.get_alignment_length())
```

    Alignment length 52 



```python
# We can shorten the printed alignments to not print the "..." before the last three letters
for record in alignment:
    print("%s - %s" % (record.seq, record.id))
```

    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA - COATB_BPIKE/30-81
    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA - Q9T0Q8_BPIKE/1-52
    DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA - COATB_BPI22/32-83
    AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA - COATB_BPM13/24-72
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA - COATB_BPZJ2/1-49
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA - Q9T0Q9_BPFD/1-49
    FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA - COATB_BPIF1/22-73



```python
# We can print the references and id of our sequences
for record in alignment:
    if record.dbxrefs:
        print("%s %s" % (record.id, record.dbxrefs))
```

    COATB_BPIKE/30-81 ['PDB; 1ifl ; 1-52;']
    COATB_BPM13/24-72 ['PDB; 2cpb ; 1-49;', 'PDB; 2cps ; 1-49;']
    Q9T0Q9_BPFD/1-49 ['PDB; 1nh4 A; 1-49;']
    COATB_BPIF1/22-73 ['PDB; 1ifk ; 1-50;']



```python
# We can see the annotations for all of the alignments
for record in alignment:
    print(record)
```

    ID: COATB_BPIKE/30-81
    Name: COATB_BPIKE
    Description: COATB_BPIKE/30-81
    Database cross-references: PDB; 1ifl ; 1-52;
    Number of features: 0
    /accession=P03620.1
    /start=30
    /end=81
    Per letter annotation for: secondary_structure
    Seq('AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA')
    ID: Q9T0Q8_BPIKE/1-52
    Name: Q9T0Q8_BPIKE
    Description: Q9T0Q8_BPIKE/1-52
    Number of features: 0
    /accession=Q9T0Q8.1
    /start=1
    /end=52
    Seq('AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA')
    ID: COATB_BPI22/32-83
    Name: COATB_BPI22
    Description: COATB_BPI22/32-83
    Number of features: 0
    /accession=P15416.1
    /start=32
    /end=83
    Seq('DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA')
    ID: COATB_BPM13/24-72
    Name: COATB_BPM13
    Description: COATB_BPM13/24-72
    Database cross-references: PDB; 2cpb ; 1-49;, PDB; 2cps ; 1-49;
    Number of features: 0
    /accession=P69541.1
    /start=24
    /end=72
    Per letter annotation for: secondary_structure
    Seq('AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA')
    ID: COATB_BPZJ2/1-49
    Name: COATB_BPZJ2
    Description: COATB_BPZJ2/1-49
    Number of features: 0
    /accession=P03618.1
    /start=1
    /end=49
    Seq('AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA')
    ID: Q9T0Q9_BPFD/1-49
    Name: Q9T0Q9_BPFD
    Description: Q9T0Q9_BPFD/1-49
    Database cross-references: PDB; 1nh4 A; 1-49;
    Number of features: 0
    /accession=Q9T0Q9.1
    /start=1
    /end=49
    Per letter annotation for: secondary_structure
    Seq('AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA')
    ID: COATB_BPIF1/22-73
    Name: COATB_BPIF1
    Description: COATB_BPIF1/22-73
    Database cross-references: PDB; 1ifk ; 1-50;
    Number of features: 0
    /accession=P03619.2
    /start=22
    /end=73
    Per letter annotation for: secondary_structure
    Seq('FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA')



```python
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
```


```python
# We can write the lines of nucleotides needed for our new file
align1 = MultipleSeqAlignment(
    [
        SeqRecord(Seq("ACTGCTAGCTAG"), id="Alpha"),
        SeqRecord(Seq("ACT-CTAGCTAG"), id="Beta"),
        SeqRecord(Seq("ACTGCTAGDTAG"), id="Gamma"),
    ]
)
align2 = MultipleSeqAlignment(
    [
        SeqRecord(Seq("GTCAGC-AG"), id="Delta"),
        SeqRecord(Seq("GACAGCTAG"), id="Epsilon"),
        SeqRecord(Seq("GTCAGCTAG"), id="Zeta"),
    ]
)
align3 = MultipleSeqAlignment(
    [
        SeqRecord(Seq("ACTAGTACAGCTG"), id="Eta"),
        SeqRecord(Seq("ACTAGTACAGCT-"), id= "Theta"),
        SeqRecord(Seq("-CTACTACAGGTG"), id="Iota"),
    ]
)
```


```python
# We can define what variables go into our "my_alignments" file
my_alignments = [align1, align2, align3]
```


```python
my_alignments
```




    [<<class 'Bio.Align.MultipleSeqAlignment'> instance (3 records of length 12) at 7f8194d79d50>,
     <<class 'Bio.Align.MultipleSeqAlignment'> instance (3 records of length 9) at 7f8194d79350>,
     <<class 'Bio.Align.MultipleSeqAlignment'> instance (3 records of length 13) at 7f8194d79450>]




```python
print(my_alignments)
```

    [<<class 'Bio.Align.MultipleSeqAlignment'> instance (3 records of length 12) at 7f8194d79d50>, <<class 'Bio.Align.MultipleSeqAlignment'> instance (3 records of length 9) at 7f8194d79350>, <<class 'Bio.Align.MultipleSeqAlignment'> instance (3 records of length 13) at 7f8194d79450>]



```python
# We can create our new file using our"my_alignments" variable
from Bio import AlignIO
AlignIO.write(my_alignments, "my_example.phy", "phylip")
```




    3




```python
# We can read in our multiple alignment files
alignments = AlignIO.parse("my_example.phy", "phylip")
```


```python
# We can load in the file we just created
for alignment in alignments:
    print(alignment)
    print()
```

    Alignment with 3 rows and 12 columns
    ACTGCTAGCTAG Alpha
    ACT-CTAGCTAG Beta
    ACTGCTAGDTAG Gamma
    
    Alignment with 3 rows and 9 columns
    GTCAGC-AG Delta
    GACAGCTAG Epsilon
    GTCAGCTAG Zeta
    
    Alignment with 3 rows and 13 columns
    ACTAGTACAGCTG Eta
    ACTAGTACAGCT- Theta
    -CTACTACAGGTG Iota
    



```python
# We can load the file in as a list
alignments = list(AlignIO.parse("my_example.phy", "phylip"))
```


```python
# We can tell it to print just the last align in our file
last_align = alignments[-1]
```


```python
print(last_align)
```

    Alignment with 3 rows and 13 columns
    ACTAGTACAGCTG Eta
    ACTAGTACAGCT- Theta
    -CTACTACAGGTG Iota



```python
# We can also just print the first align in our file
first_align = alignments[0]
```


```python
print(first_align)
```

    Alignment with 3 rows and 12 columns
    ACTGCTAGCTAG Alpha
    ACT-CTAGCTAG Beta
    ACTGCTAGDTAG Gamma


## Part 2
```python
from Bio import AlignIO
```


```python
# We can convert our file from an sth file to a clustal file
count = AlignIO.convert("PF05371_seed.sth.txt", "stockholm", "PF05371_seed.aln.txt", "clustal")
```


```python
# We can print the number of converted alignments
print("Converted %i alignments" % count)
```

    Converted 1 alignments



```python
# We can convert our file using a different method where we utilize both the "parse" and "write" functions
alignments = AlignIO.parse("PF05371_seed.sth.txt", "stockholm")
```


```python
count = AlignIO.write(alignments, "PF05371_seed.aln.txt", "clustal")
```


```python
# We can print out the converted alignments from this other method
print("Converted %i alignments" % count)
```

    Converted 1 alignments



```python
# We can also convert our file into a phylip format
AlignIO.convert("PF05371_seed.sth.txt", "stockholm", "PF05371_seed.phy", "phylip")
```




    1




```python
# We can also convert our file into a phylip-relaxed format which puts in the sample name
AlignIO.convert("PF05371_seed.sth.txt", "stockholm", "PF05371_seed.phy", "phylip-relaxed")
```




    1




```python
# We can manipulate the names of the files and number the sequences
alignment = AlignIO.read("PF05371_seed.sth.txt", "stockholm")
name_mapping = {}
for i, record in enumerate(alignment):
    name_mapping[i] = record.id
    record.id = "seq%i" % i
```


```python
print(name_mapping)
```

    {0: 'COATB_BPIKE/30-81', 1: 'Q9T0Q8_BPIKE/1-52', 2: 'COATB_BPI22/32-83', 3: 'COATB_BPM13/24-72', 4: 'COATB_BPZJ2/1-49', 5: 'Q9T0Q9_BPFD/1-49', 6: 'COATB_BPIF1/22-73'}



```python
# We can add the sequence numbers to the file
AlignIO.write([alignment], "PF05371_seed.phy.txt", "phylip")
```




    1




```python
# We can read in our stockholm file
alignment = AlignIO.read("PF05371_seed.sth.txt", "stockholm")
```


```python
# We can print the number of rows in the alignment
print("Number of rows: %i" % len(alignment))
```

    Number of rows: 7



```python
# We can print the sequences and IDs in the file
for record in alignment:
    print("%s - %s" % (record.seq, record.id))
```

    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA - COATB_BPIKE/30-81
    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA - Q9T0Q8_BPIKE/1-52
    DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA - COATB_BPI22/32-83
    AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA - COATB_BPM13/24-72
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA - COATB_BPZJ2/1-49
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA - Q9T0Q9_BPFD/1-49
    FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA - COATB_BPIF1/22-73



```python
# We can print certain numbered sequences such as the 3rd-7th
print(alignment[3:7])
```

    Alignment with 4 rows and 52 columns
    AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA COATB_BPM13/24-72
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA COATB_BPZJ2/1-49
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA Q9T0Q9_BPFD/1-49
    FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKL...SRA COATB_BPIF1/22-73



```python
# We can print a specific nucleotide from within a column and row
print(alignment[2, 6])
```

    T



```python
# We can print a specific nucleotide in a different way
print(alignment[2].seq[6])
```

    T



```python
# We can print a specific column
print(alignment[:, 6])
```

    TTT---T



```python
# We can select a range within the alignment
print(alignment[3:6, :6])
```

    Alignment with 3 rows and 6 columns
    AEGDDP COATB_BPM13/24-72
    AEGDDP COATB_BPZJ2/1-49
    AEGDDP Q9T0Q9_BPFD/1-49



```python
# We can use the ":" to tell it to print everything
print(alignment[:, :6])
```

    Alignment with 7 rows and 6 columns
    AEPNAA COATB_BPIKE/30-81
    AEPNAA Q9T0Q8_BPIKE/1-52
    DGTSTA COATB_BPI22/32-83
    AEGDDP COATB_BPM13/24-72
    AEGDDP COATB_BPZJ2/1-49
    AEGDDP Q9T0Q9_BPFD/1-49
    FAADDA COATB_BPIF1/22-73



```python
print(alignment[:, 6:9])
```

    Alignment with 7 rows and 3 columns
    TNY COATB_BPIKE/30-81
    TNY Q9T0Q8_BPIKE/1-52
    TSY COATB_BPI22/32-83
    --- COATB_BPM13/24-72
    --- COATB_BPZJ2/1-49
    --- Q9T0Q9_BPFD/1-49
    TSQ COATB_BPIF1/22-73



```python
# We can print everything after a certain point
print(alignment[:, 9:])
```

    Alignment with 7 rows and 43 columns
    ATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA COATB_BPIKE/30-81
    ATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA Q9T0Q8_BPIKE/1-52
    ATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA COATB_BPI22/32-83
    AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA COATB_BPM13/24-72
    AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA COATB_BPZJ2/1-49
    AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA Q9T0Q9_BPFD/1-49
    AKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA COATB_BPIF1/22-73



```python
# We can remove blocks of columns
edited = alignment[:, :6] + alignment [:, 9:]
```


```python
# We can print the edited version
print(edited)
```

    Alignment with 7 rows and 49 columns
    AEPNAAATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA COATB_BPIKE/30-81
    AEPNAAATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA Q9T0Q8_BPIKE/1-52
    DGTSTAATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA COATB_BPI22/32-83
    AEGDDPAKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA COATB_BPM13/24-72
    AEGDDPAKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA COATB_BPZJ2/1-49
    AEGDDPAKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA Q9T0Q9_BPFD/1-49
    FAADDAAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA COATB_BPIF1/22-73



```python
# We can sort based on the ID 
edited.sort()
```


```python
print(edited)
```

    Alignment with 7 rows and 49 columns
    DGTSTAATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA COATB_BPI22/32-83
    FAADDAAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA COATB_BPIF1/22-73
    AEPNAAATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA COATB_BPIKE/30-81
    AEGDDPAKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA COATB_BPM13/24-72
    AEGDDPAKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA COATB_BPZJ2/1-49
    AEPNAAATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA Q9T0Q8_BPIKE/1-52
    AEGDDPAKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA Q9T0Q9_BPFD/1-49


## Part 3
```python
# We can load in our required packages
from Bio.Seq import Seq
```


```python
from Bio.SeqRecord import SeqRecord
```


```python
from Bio.Align import MultipleSeqAlignment
```


```python
# We can write the lines of nucleotides needed for our new file
alignment = MultipleSeqAlignment(
    [
        SeqRecord(Seq("ACTCCTA"), id="seq1"),
        SeqRecord(Seq("AAT-CTA"), id="seq2"),
        SeqRecord(Seq("CCTACT-"), id="seq3"),
        SeqRecord(Seq("TCTCCTC"), id="seq4"),
    ]
)
```


```python
# We can print our new alignment
print(alignment)
```

    Alignment with 4 rows and 7 columns
    ACTCCTA seq1
    AAT-CTA seq2
    CCTACT- seq3
    TCTCCTC seq4



```python
# We can calculate the number of different letter pairs
substitutions = alignment.substitutions
```


```python
# We can print out the number of different letter pairs
print(substitutions)
```

        A    C    T
    A 2.0  4.5  1.0
    C 4.5 10.0  0.5
    T 1.0  0.5 12.0
    



```python
# We can add a fourth letter in our substitution matrix
m = substitutions.select("ATCG")
```


```python
print(m)
```

        A    T    C   G
    A 2.0  1.0  4.5 0.0
    T 1.0 12.0  0.5 0.0
    C 4.5  0.5 10.0 0.0
    G 0.0  0.0  0.0 0.0
    



```python
# We can change the order of our substitutions matrix
m = substitutions.select("ACTG")
```


```python
print(m)
```

        A    C    T   G
    A 2.0  4.5  1.0 0.0
    C 4.5 10.0  0.5 0.0
    T 1.0  0.5 12.0 0.0
    G 0.0  0.0  0.0 0.0
    



```python
# We can download our necessary applications
import Bio.Align.Applications
```


```python
# We can access the directory
dir(Bio.Align.Applications)
```




    ['ClustalOmegaCommandline',
     'ClustalwCommandline',
     'DialignCommandline',
     'MSAProbsCommandline',
     'MafftCommandline',
     'MuscleCommandline',
     'PrankCommandline',
     'ProbconsCommandline',
     'TCoffeeCommandline',
     '_ClustalOmega',
     '_Clustalw',
     '_Dialign',
     '_MSAProbs',
     '_Mafft',
     '_Muscle',
     '_Prank',
     '_Probcons',
     '_TCoffee',
     '__all__',
     '__builtins__',
     '__cached__',
     '__doc__',
     '__file__',
     '__loader__',
     '__name__',
     '__package__',
     '__path__',
     '__spec__']




```python
# We can import our ClustalwCommandLine
from Bio.Align.Applications import ClustalwCommandline
```


```python
from Bio import AlignIO
```


```python
# We can import our opuntia file
align = AlignIO.read("opuntia.aln", "clustal")
```


```python
# We can print our read in file
print(align)
```

    Alignment with 7 rows and 906 columns
    TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273285|gb|AF191659.1|AF191
    TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273284|gb|AF191658.1|AF191
    TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273287|gb|AF191661.1|AF191
    TATACATAAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273286|gb|AF191660.1|AF191
    TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273290|gb|AF191664.1|AF191
    TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273289|gb|AF191663.1|AF191
    TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273291|gb|AF191665.1|AF191



```python
# We need to import our necessary files
from Bio import Phylo
```


```python
# We can create a phylogenic tree from our file information
tree = Phylo.read("opuntia.dnd", "newick")
```


```python
# We can print out our phylogenic tree
Phylo.draw_ascii(tree)
```

                                 _______________ gi|6273291|gb|AF191665.1|AF191665
      __________________________|
     |                          |   ______ gi|6273290|gb|AF191664.1|AF191664
     |                          |__|
     |                             |_____ gi|6273289|gb|AF191663.1|AF191663
     |
    _|_________________ gi|6273287|gb|AF191661.1|AF191661
     |
     |__________ gi|6273286|gb|AF191660.1|AF191660
     |
     |    __ gi|6273285|gb|AF191659.1|AF191659
     |___|
         | gi|6273284|gb|AF191658.1|AF191658
    



```python

```
## Part 4
```python
# We need to import the necessary functions
from Bio import Align
```


```python
# We need to create our aligner variable
aligner = Align.PairwiseAligner()
```


```python
# We can tailor our match score
aligner = Align.PairwiseAligner(match_score = 1.0)
```


```python
# We can create our target sequence
target = "GAACT"
```


```python
# We can create our query sequence
query = "GAT"
```


```python
# We can compare the match score of the query and target sequences
score = aligner.score(target, query)
```


```python
# We can print our match score
score
```




    3.0




```python
alignments = aligner.align(target, query)
```


```python
# We can print the target and query sequences to see their alignment for calculating match score
for alignment in alignments:
    print(alignment)
```

    target            1 GAACT 6
                      0 ||||| 5
    query             0 GAACT 5
    



```python
# We can change the mode of aligner
aligner.mode = "local"
```


```python
# We can create a new target sequence
target = "AGAACTC"
```


```python
# We can create a new query sequence
query = "GAACT"
```


```python
# We can calculate the match score of the new target and query sequences
score = aligner.score(target, query)
```


```python
# We can print the match score
score
```




    5.0




```python
alignments = aligner.align(target, query)
```


```python
# We can print the alignment for calculating our match score
for alignment in alignments:
    print(alignment)
```

    target            1 GAACT 6
                      0 ||||| 5
    query             0 GAACT 5
    



```python
# We can see our scores for many different alignments
print(aligner)
```

    Pairwise sequence aligner with parameters
      wildcard: None
      match_score: 1.000000
      mismatch_score: 0.000000
      target_internal_open_gap_score: 0.000000
      target_internal_extend_gap_score: 0.000000
      target_left_open_gap_score: 0.000000
      target_left_extend_gap_score: 0.000000
      target_right_open_gap_score: 0.000000
      target_right_extend_gap_score: 0.000000
      query_internal_open_gap_score: 0.000000
      query_internal_extend_gap_score: 0.000000
      query_left_open_gap_score: 0.000000
      query_left_extend_gap_score: 0.000000
      query_right_open_gap_score: 0.000000
      query_right_extend_gap_score: 0.000000
      mode: local
    



```python
# We can see what algorithm we are using
aligner.algorithm
```




    'Smith-Waterman'




```python
# We can see a significant score
aligner.epsilon
```




    1e-06




```python
from Bio import Align
```


```python
aligner = Align.PairwiseAligner()
```


```python
# We can create a new target sequence
target = "GAACT"
```


```python
# We can create a new query sequence
query = "GAT"
```


```python
# We can see the score for our target and query sequences
alignments = aligner.align(target, query)
```


```python
# We can set our alignment to start at the first sequence
alignment = alignments[0]
```


```python
# We can see the alignment of our object
alignment
```




    <Alignment object (2 rows x 5 columns) at 0x7fc4085f3350>




```python
# We can look at the alignment score
alignment.score
```




    3.0




```python
# We can look at what the target was
alignment.target
```




    'GAACT'




```python
# We can look at what the query was
alignment.query
```




    'GAT'




```python
# We can see the alignment
print(alignment)
```

    target            0 GAACT 5
                      0 ||--| 5
    query             0 GA--T 3
    



```python
# We can see the alignment coordinates
alignment.coordinates
```




    array([[0, 2, 4, 5],
           [0, 2, 2, 3]])




```python
# We can find the length of the alignment
len(alignment)
```




    2




```python
# We can look at the shape of the alignment
alignment.shape
```




    (2, 5)




```python
# We can set the mode of the alignment
aligner.mode = "local"
```


```python
# We can create a new variable
local_alignments = aligner.align("TGAACT", "GAC")
```


```python
# We can tell our local alignment to start at the first sequence
local_alignment = local_alignments[0]
```


```python
# We can print our local alignment
print(local_alignment)
```

    target            1 GAAC 5
                      0 ||-| 4
    query             0 GA-C 3
    



```python
# We can look at the shape of the local alignment
local_alignment.shape
```




    (2, 4)




```python
# We can set the object to the global mode
aligner.mode = "global"
```


```python
# We can set the match and mismatch score
aligner = Align.PairwiseAligner(match = 1.0, mismatch_score = -10)
```


```python
# We can print the scores of our aligner
print(aligner)
```

    Pairwise sequence aligner with parameters
      wildcard: None
      match_score: 1.000000
      mismatch_score: -10.000000
      target_internal_open_gap_score: 0.000000
      target_internal_extend_gap_score: 0.000000
      target_left_open_gap_score: 0.000000
      target_left_extend_gap_score: 0.000000
      target_right_open_gap_score: 0.000000
      target_right_extend_gap_score: 0.000000
      query_internal_open_gap_score: 0.000000
      query_internal_extend_gap_score: 0.000000
      query_left_open_gap_score: 0.000000
      query_left_extend_gap_score: 0.000000
      query_right_open_gap_score: 0.000000
      query_right_extend_gap_score: 0.000000
      mode: global
    



```python
# We can set new alignment sequences
alignments = aligner.align("AAACAAA", "AAAGAAA")
```


```python
# We can get the length of our new alignments
len(alignments)
```




    2




```python
# We can print our first alignment
print(alignments[0])
```

    target            0 AAAC-AAA 7
                      0 |||--||| 8
    query             0 AAA-GAAA 7
    



```python
# We can print the alignments from the second alignment score
print(alignments[1])
```

    target            0 AAA-CAAA 7
                      0 |||--||| 8
    query             0 AAAG-AAA 7
    



```python
print(local_alignment)
```

    target            1 GAAC 5
                      0 ||-| 4
    query             0 GA-C 3
    



```python
# We can sort our local alignment
local_alignment.sort()
```


```python
# We can print our local alignment
print(local_alignment)
```

    target            0 GA-C 3
                      0 ||-| 4
    query             1 GAAC 5
    


## Part 5
```python
from Bio import Align
```


```python
# We can import the reverse_complement function
from Bio.Seq import reverse_complement
```


```python
# We can create a new target sequence
target = "AAACCC"
```


```python
# We can create a new query sequence
query = "AACC"
```


```python
# We can set our mismatch and internal gap score for our alignment
aligner = Align.PairwiseAligner(mismatch_score = -1, internal_gap_score = -1)
```


```python
# We can print our scores for the target and query
aligner.score(target, query)
```




    4.0




```python
# We can get our score between our target, reverse complement, and query
aligner.score(target, reverse_complement(query))
```




    0.0




```python
# We can get our score for our reverse complement and the negative strand
aligner.score(target, reverse_complement(query), strand = "-")
```




    4.0




```python
# We can get our aligner score for our target and query with the negative strand
aligner.score(target, query, strand = "-")
```




    0.0




```python
# We can save our alignments
alignments = aligner.align(target, query)
```


```python
# We can see the length of our alignments
len(alignments)
```




    1




```python
# We can print our first alignment
print(alignments[0])
```

    target            0 AAACCC 6
                      0 -||||- 6
    query             0 -AACC- 4
    



```python
# We can save our alignment to a different alignment format
print(alignments[0].format("bed"))
```

    target	1	5	query	4.0	+	1	5	0	1	4,	0,
    



```python
# We can align our alignments along the reverse_complement and the negative strand
alignments = aligner.align(target, reverse_complement(query), strand = "-")
```


```python
# We can format our first alignment to a bed format
print(alignments[0].format("bed"))
```

    target	1	5	query	4.0	-	1	5	0	1	4,	0,
    



```python
# We can just do the negative strand version without the negative complement
alignments = aligner.align(target, query, strand = "-")
```


```python
# We can show the length of our alignments
len(alignments)
```




    2




```python
# We can print our first alignment
print(alignments[0])
```

    target            0 AAACCC----  6
                      0 ---------- 10
    query             4 ------GGTT  0
    



```python
# We can print our second alignment
print(alignments[1])
```

    target            0 ----AAACCC  6
                      0 ---------- 10
    query             4 GGTT------  0
    



```python
# We can assign our left gap score
aligner.left_gap_score = -0.5
```


```python
# We can assign our right gap score
aligner.right_gap_score = -0.2
```


```python
# We can print the aligner score for our target and query
aligner.score(target, query)
```




    3.3




```python
# We can redefine our alignments variable
alignments = aligner.align(target,query)
```


```python
# We can print the length of our alignments
len(alignments)
```




    1




```python
# We can print our only alignment
print(alignments[0])
```

    target            0 AAACCC 6
                      0 -||||- 6
    query             0 -AACC- 4
    



```python
# We can compare our alignment to the reverse complement for the alignment score
alignments = aligner.align(target, reverse_complement(query), strand = "-")
```


```python
# We can print our alignment score aligned of the reverse complement
aligner.score(target, reverse_complement(query), strand = "-")
```




    3.3




```python
# We can print the first alignment
print(alignments[0])
```

    target            0 AAACCC 6
                      0 -||||- 6
    query             4 -AACC- 0
    



```python
# We can calculate the alignment score with the reverse complement and the positive strand
aligner.score(target, reverse_complement(query), strand = "+")
```




    -2.0




```python

```


# Challenge + Blast
## Blast
```python
# We can import our necessary functions
from Bio.Blast import NCBIWWW
```


```python
# We can connect our NCBI email
NCBIWWW.email = "dr.josh.vandenbrink@gmail.com"
```


```python
# We can blast genes in the NCBI database
result_handle = NCBIWWW.qblast("blastn", "nt", "8332116")
```


```python
# We can import our necessary files
from Bio import SeqIO
```


```python
# We can import our fasta file
record = SeqIO.read("m_cold.fasta", format = "fasta" )
```


```python
# We can print our imported fasta file
print(record)
```


```python
# We can blast our own file
result_handle = NCBIWWW.qblast("blastn", "nt", record.seq)
```


```python
# We can open our result handle
with open("m_cold.fasta", "w") as out_handle:
    out_handle.write(result_handle.read())
result_handle.close()
```


```python
# We can save it as an XML file
from Bio.Blast import NCBIXML
```


```python
# We can reopen our result handle
result_handle = open("m_cold.fasta")
```


```python
# We can read in our result handle
blast_record = NCBIXML.read(result_handle)
```


```python
# We can create an E value
E_VALUE_THRESH = 0.04
```


```python
# We can print the results of our blast record
for alignment in blast_record.alignments:
    for hsp in alignment.hsps:
        if hsp.expect < E_VALUE_THRESH:
            print("****ALIGHTMENT****")
            print("sequence:", alignment.title)
            print("length:", alignment.length)
            print("e value:", hsp.expect)
            print(hsp.query[0:75] + "...")
            print(hsp.match[0:75] + "...")
            print(hsp.sbjct[0:75] + "...")
```

## Challenge 1
```python
# We can import our necessary NCBI files
from Bio.Blast import NCBIWWW
```


```python
# We can connect our NCBI email
NCBIWWW.email = "dr.josh.vandenbrink@gmail.com"
```


```python
from Bio import SeqIO
```


```python
# We can upload our file of the KRAS gene
record = SeqIO.read("chek2.fasta.txt", format = "fasta")
```


```python
# We can print the record of the gene
print(record)
```

    ID: lcl|NC_000022.11_cds_XP_011528144.1_1
    Name: lcl|NC_000022.11_cds_XP_011528144.1_1
    Description: lcl|NC_000022.11_cds_XP_011528144.1_1 [gene=CHEK2] [db_xref=GeneID:11200] [protein=serine/threonine-protein kinase Chk2 isoform X10] [protein_id=XP_011528144.1] [location=complement(join(28687897..28687986,28689135..28689215,28694032..28694117,28695127..28695242,28695710..28695873,28696901..28696987,28699838..28699937,28703505..28703566,28710006..28710059,28711909..28712017,28724886..28724923,28725243..28725367,28734403..28734727,28737260..28737283))] [gbkey=CDS]
    Number of features: 0
    Seq('ATGCTGAGGCTGCTAACGTCTATGGTCGTGATGTCTCGGGAGTCGGATGTTGAG...TGA')



```python
# We can blast the record of the KRAS gene using the NCBI database
result_handle = NCBIWWW.qblast("blastn", "nt", record.seq)
```


```python
# We can open and close our result handle for the KRAS file
with open("chek2.fasta.txt", "w") as out_handle:
    out_handle.write(result_handle.read())
result_handle.close()
```


    ---------------------------------------------------------------------------

    PermissionError                           Traceback (most recent call last)

    <ipython-input-66-79f9babce1d8> in <module>
          1 # We can open and close our result handle for the KRAS file
    ----> 2 with open("chek2.fasta.txt", "w") as out_handle:
          3     out_handle.write(result_handle.read())
          4 result_handle.close()


    PermissionError: [Errno 13] Permission denied: 'chek2.fasta.txt'



```python
# We can reopen our result handle for the KRAS file
result_handle = open("chek2.fasta.txt")
```


```python
from Bio.Blast import NCBIXML
```


```python
blast_record = NCBIXML.read(result_handle)
```


```python
# We can set our E value threshold
E_VALUE_THRESH = 0.04
```


```python
# We can print the blast results for the KRAS gene
for alignment in blast_record.alignments:
    for hsp in alignment.hsps:
        if hsp.expect < E_VALUE_THRESH:
            print("****ALIGHTMENT****")
            print("sequence:", alignment.title)
            print("length:", alignment.length)
            print("e value:", hsp.expect)
            print(hsp.query[0:75] + "...")
            print(hsp.match[0:75] + "...")
            print(hsp.sbjct[0:75] + "...")
```

    ****ALIGHTMENT****
    sequence: gi|2217338827|ref|XM_011529842.3| PREDICTED: Homo sapiens checkpoint kinase 2 (CHEK2), transcript variant X10, mRNA >gi|2462583673|ref|XM_054325046.1| PREDICTED: Homo sapiens checkpoint kinase 2 (CHEK2), transcript variant X10, mRNA
    length: 1674
    e value: 0.0
    ATGCTGAGGCTGCTAACGTCTATGGTCGTGATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGC...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    ATGCTGAGGCTGCTAACGTCTATGGTCGTGATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGC...
    ****ALIGHTMENT****
    sequence: gi|1675009750|ref|NM_001349956.2| Homo sapiens checkpoint kinase 2 (CHEK2), transcript variant 5, mRNA
    length: 1643
    e value: 0.0
    GGTCGTGATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGCAGTGCCTGTTCACAGCCCCATGG...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    GGTCGTGATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGCAGTGCCTGTTCACAGCCCCATGG...
    ****ALIGHTMENT****
    sequence: gi|2493398256|ref|XM_055374091.1| PREDICTED: Gorilla gorilla gorilla checkpoint kinase 2 (CHEK2), transcript variant X5, mRNA
    length: 1849
    e value: 0.0
    GGTCGTGATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGCAGTGCCTGTTCACAGCCCCATGG...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    GGTCGTGATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGCAGTGCCTGTTCACAGCCCCATGG...
    ****ALIGHTMENT****
    sequence: gi|1743222572|ref|XM_030816023.1| PREDICTED: Nomascus leucogenys checkpoint kinase 2 (CHEK2), transcript variant X4, mRNA
    length: 1688
    e value: 0.0
    GGTCGTGATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGCAGTGCCTGTTCACAGCCCCATGG...
    |||||||||||||||||||||||||||||||||||||||||| ||||||||||||||||||||||||||||||||...
    GGTCGTGATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCCCATGGCAGCAGTGCCTGTTCACAGCCCCATGG...
    ****ALIGHTMENT****
    sequence: gi|2491616672|ref|XM_055252869.1| PREDICTED: Symphalangus syndactylus checkpoint kinase 2 (CHEK2), transcript variant X4, mRNA
    length: 1642
    e value: 0.0
    GGTCGTGATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGCAGTGCCTGTTCACAGCCCCATGG...
    |||||||||||||||||||||||||||||||||||||||||| ||||||||||||||||||||||||||||||||...
    GGTCGTGATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCCCATGGCAGCAGTGCCTGTTCACAGCCCCATGG...
    ****ALIGHTMENT****
    sequence: gi|1411089960|ref|XM_025400561.1| PREDICTED: Theropithecus gelada checkpoint kinase 2 (CHEK2), transcript variant X5, mRNA
    length: 1685
    e value: 0.0
    GGTCGTGATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGCAGTGCCTGTTCACAGCCCCATGG...
    |||||||||||||||||||||||||||||||||||||||||| |||||||||||| |||||||||||||||||||...
    GGTCGTGATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCCCATGGCAGCAGTACCTGTTCACAGCCCCATGG...
    ****ALIGHTMENT****
    sequence: gi|795227094|ref|XM_011995735.1| PREDICTED: Mandrillus leucophaeus checkpoint kinase 2 (CHEK2), mRNA
    length: 1631
    e value: 0.0
    GGTCGTGATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGCAGTGCCTGTTCACAGCCCCATGG...
    |||||||||||||||||||||||||||||||||||||||||| |||||||||||| |||||||||||||||||||...
    GGTCGTGATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCCCATGGCAGCAGTACCTGTTCACAGCCCCATGG...
    ****ALIGHTMENT****
    sequence: gi|2217338824|ref|XM_047441107.1| PREDICTED: Homo sapiens checkpoint kinase 2 (CHEK2), transcript variant X8, mRNA >gi|2462583669|ref|XM_054325044.1| PREDICTED: Homo sapiens checkpoint kinase 2 (CHEK2), transcript variant X8, mRNA
    length: 1807
    e value: 0.0
    TTTGCCAATCTTGAATGTGTGAATGACAACTACTGGTTTGGGAGGGACAAAAGCTGTGAATATTGCTTTGATGAA...
    ||||| | |||  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    TTTGC-AGTCTAAAATGTGTGAATGACAACTACTGGTTTGGGAGGGACAAAAGCTGTGAATATTGCTTTGATGAA...
    ****ALIGHTMENT****
    sequence: gi|2217338824|ref|XM_047441107.1| PREDICTED: Homo sapiens checkpoint kinase 2 (CHEK2), transcript variant X8, mRNA >gi|2462583669|ref|XM_054325044.1| PREDICTED: Homo sapiens checkpoint kinase 2 (CHEK2), transcript variant X8, mRNA
    length: 1807
    e value: 8.97197e-176
    ATGCTGAGGCTGCTAACGTCTATGGTCGTGATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGC...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    ATGCTGAGGCTGCTAACGTCTATGGTCGTGATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGC...
    ****ALIGHTMENT****
    sequence: gi|2493398252|ref|XM_055374090.1| PREDICTED: Gorilla gorilla gorilla checkpoint kinase 2 (CHEK2), transcript variant X4, mRNA
    length: 1978
    e value: 0.0
    TTTGCCAATCTTGAATGTGTGAATGACAACTACTGGTTTGGGAGGGACAAAAGCTGTGAATATTGCTTTGATGAA...
    ||||| | |||  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    TTTGC-AGTCTAAAATGTGTGAATGACAACTACTGGTTTGGGAGGGACAAAAGCTGTGAATATTGCTTTGATGAA...
    ****ALIGHTMENT****
    sequence: gi|2493398252|ref|XM_055374090.1| PREDICTED: Gorilla gorilla gorilla checkpoint kinase 2 (CHEK2), transcript variant X4, mRNA
    length: 1978
    e value: 2.74697e-163
    GGTCGTGATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGCAGTGCCTGTTCACAGCCCCATGG...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    GGTCGTGATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGCAGTGCCTGTTCACAGCCCCATGG...
    ****ALIGHTMENT****
    sequence: gi|2462583677|ref|XM_054325048.1| PREDICTED: Homo sapiens checkpoint kinase 2 (CHEK2), transcript variant X12, mRNA
    length: 1509
    e value: 0.0
    ATGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCA...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    ATGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCA...
    ****ALIGHTMENT****
    sequence: gi|2217338829|ref|XM_006724114.4| PREDICTED: Homo sapiens checkpoint kinase 2 (CHEK2), transcript variant X12, mRNA
    length: 1327
    e value: 0.0
    ATGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCA...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    ATGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCA...
    ****ALIGHTMENT****
    sequence: gi|1675053443|ref|NM_001005735.2| Homo sapiens checkpoint kinase 2 (CHEK2), transcript variant 3, mRNA
    length: 1973
    e value: 0.0
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    ****ALIGHTMENT****
    sequence: gi|1675053443|ref|NM_001005735.2| Homo sapiens checkpoint kinase 2 (CHEK2), transcript variant 3, mRNA
    length: 1973
    e value: 2.74697e-163
    GGTCGTGATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGCAGTGCCTGTTCACAGCCCCATGG...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    GGTCGTGATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGCAGTGCCTGTTCACAGCCCCATGG...
    ****ALIGHTMENT****
    sequence: gi|1675053443|ref|NM_001005735.2| Homo sapiens checkpoint kinase 2 (CHEK2), transcript variant 3, mRNA
    length: 1973
    e value: 2.94854e-55
    TTTGCCAATCTTGAATGTGTGAATGACAACTACTGGTTTGGGAGGGACAAAAGCTGTGAATATTGCTTTGATGAA...
    ||||| | |||  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    TTTGC-AGTCTAAAATGTGTGAATGACAACTACTGGTTTGGGAGGGACAAAAGCTGTGAATATTGCTTTGATGAA...
    ****ALIGHTMENT****
    sequence: gi|1674986356|ref|NM_001257387.2| Homo sapiens checkpoint kinase 2 (CHEK2), transcript variant 4, mRNA
    length: 1958
    e value: 0.0
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    ****ALIGHTMENT****
    sequence: gi|1674986356|ref|NM_001257387.2| Homo sapiens checkpoint kinase 2 (CHEK2), transcript variant 4, mRNA
    length: 1958
    e value: 0.0
    GGTCGTGATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGCAGTGCCTGTTCACAGCCCCATGG...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    GGTCGTGATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGCAGTGCCTGTTCACAGCCCCATGG...
    ****ALIGHTMENT****
    sequence: gi|1402055726|ref|NM_007194.4| Homo sapiens checkpoint kinase 2 (CHEK2), transcript variant 1, mRNA
    length: 1844
    e value: 0.0
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    ****ALIGHTMENT****
    sequence: gi|1402055726|ref|NM_007194.4| Homo sapiens checkpoint kinase 2 (CHEK2), transcript variant 1, mRNA
    length: 1844
    e value: 0.0
    GGTCGTGATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGCAGTGCCTGTTCACAGCCCCATGG...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    GGTCGTGATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGCAGTGCCTGTTCACAGCCCCATGG...
    ****ALIGHTMENT****
    sequence: gi|2462583681|ref|XM_054325050.1| PREDICTED: Homo sapiens checkpoint kinase 2 (CHEK2), transcript variant X14, mRNA
    length: 1271
    e value: 0.0
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    ****ALIGHTMENT****
    sequence: gi|2217338837|ref|XM_011529845.3| PREDICTED: Homo sapiens checkpoint kinase 2 (CHEK2), transcript variant X18, mRNA >gi|2462583686|ref|XM_054325051.1| PREDICTED: Homo sapiens checkpoint kinase 2 (CHEK2), transcript variant X18, mRNA
    length: 1533
    e value: 0.0
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    ****ALIGHTMENT****
    sequence: gi|2217338837|ref|XM_011529845.3| PREDICTED: Homo sapiens checkpoint kinase 2 (CHEK2), transcript variant X18, mRNA >gi|2462583686|ref|XM_054325051.1| PREDICTED: Homo sapiens checkpoint kinase 2 (CHEK2), transcript variant X18, mRNA
    length: 1533
    e value: 2.94854e-55
    GAATGTGTGAATGACAACTACTGGTTTGGGAGGGACAAAAGCTGTGAATATTGCTTTGATGAACCACTGCTGAAA...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    GAATGTGTGAATGACAACTACTGGTTTGGGAGGGACAAAAGCTGTGAATATTGCTTTGATGAACCACTGCTGAAA...
    ****ALIGHTMENT****
    sequence: gi|2217338833|ref|XM_006724116.3| PREDICTED: Homo sapiens checkpoint kinase 2 (CHEK2), transcript variant X14, mRNA
    length: 1271
    e value: 0.0
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    ****ALIGHTMENT****
    sequence: gi|2217338821|ref|XM_024452148.2| PREDICTED: Homo sapiens checkpoint kinase 2 (CHEK2), transcript variant X6, mRNA >gi|2462583665|ref|XM_054325042.1| PREDICTED: Homo sapiens checkpoint kinase 2 (CHEK2), transcript variant X6, mRNA
    length: 3357
    e value: 0.0
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    ****ALIGHTMENT****
    sequence: gi|2217338821|ref|XM_024452148.2| PREDICTED: Homo sapiens checkpoint kinase 2 (CHEK2), transcript variant X6, mRNA >gi|2462583665|ref|XM_054325042.1| PREDICTED: Homo sapiens checkpoint kinase 2 (CHEK2), transcript variant X6, mRNA
    length: 3357
    e value: 0.0
    ATGCTGAGGCTGCTAACGTCTATGGTCGTGATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGC...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    ATGCTGAGGCTGCTAACGTCTATGGTCGTGATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGC...
    ****ALIGHTMENT****
    sequence: gi|2217338816|ref|XM_047441104.1| PREDICTED: Homo sapiens checkpoint kinase 2 (CHEK2), transcript variant X3, mRNA >gi|2462583659|ref|XM_054325039.1| PREDICTED: Homo sapiens checkpoint kinase 2 (CHEK2), transcript variant X3, mRNA
    length: 1951
    e value: 0.0
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    ****ALIGHTMENT****
    sequence: gi|2217338816|ref|XM_047441104.1| PREDICTED: Homo sapiens checkpoint kinase 2 (CHEK2), transcript variant X3, mRNA >gi|2462583659|ref|XM_054325039.1| PREDICTED: Homo sapiens checkpoint kinase 2 (CHEK2), transcript variant X3, mRNA
    length: 1951
    e value: 2.74697e-163
    GGTCGTGATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGCAGTGCCTGTTCACAGCCCCATGG...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    GGTCGTGATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGCAGTGCCTGTTCACAGCCCCATGG...
    ****ALIGHTMENT****
    sequence: gi|2217338816|ref|XM_047441104.1| PREDICTED: Homo sapiens checkpoint kinase 2 (CHEK2), transcript variant X3, mRNA >gi|2462583659|ref|XM_054325039.1| PREDICTED: Homo sapiens checkpoint kinase 2 (CHEK2), transcript variant X3, mRNA
    length: 1951
    e value: 2.94854e-55
    GAATGTGTGAATGACAACTACTGGTTTGGGAGGGACAAAAGCTGTGAATATTGCTTTGATGAACCACTGCTGAAA...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    GAATGTGTGAATGACAACTACTGGTTTGGGAGGGACAAAAGCTGTGAATATTGCTTTGATGAACCACTGCTGAAA...
    ****ALIGHTMENT****
    sequence: gi|2217338815|ref|XM_017028560.2| PREDICTED: Homo sapiens checkpoint kinase 2 (CHEK2), transcript variant X2, mRNA >gi|2462583657|ref|XM_054325038.1| PREDICTED: Homo sapiens checkpoint kinase 2 (CHEK2), transcript variant X2, mRNA
    length: 1968
    e value: 0.0
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    ****ALIGHTMENT****
    sequence: gi|2217338815|ref|XM_017028560.2| PREDICTED: Homo sapiens checkpoint kinase 2 (CHEK2), transcript variant X2, mRNA >gi|2462583657|ref|XM_054325038.1| PREDICTED: Homo sapiens checkpoint kinase 2 (CHEK2), transcript variant X2, mRNA
    length: 1968
    e value: 8.97197e-176
    ATGCTGAGGCTGCTAACGTCTATGGTCGTGATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGC...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    ATGCTGAGGCTGCTAACGTCTATGGTCGTGATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGC...
    ****ALIGHTMENT****
    sequence: gi|2217338815|ref|XM_017028560.2| PREDICTED: Homo sapiens checkpoint kinase 2 (CHEK2), transcript variant X2, mRNA >gi|2462583657|ref|XM_054325038.1| PREDICTED: Homo sapiens checkpoint kinase 2 (CHEK2), transcript variant X2, mRNA
    length: 1968
    e value: 2.94854e-55
    GAATGTGTGAATGACAACTACTGGTTTGGGAGGGACAAAAGCTGTGAATATTGCTTTGATGAACCACTGCTGAAA...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    GAATGTGTGAATGACAACTACTGGTTTGGGAGGGACAAAAGCTGTGAATATTGCTTTGATGAACCACTGCTGAAA...
    ****ALIGHTMENT****
    sequence: gi|2217338814|ref|XM_011529839.3| PREDICTED: Homo sapiens checkpoint kinase 2 (CHEK2), transcript variant X1, mRNA >gi|2462583655|ref|XM_054325037.1| PREDICTED: Homo sapiens checkpoint kinase 2 (CHEK2), transcript variant X1, mRNA
    length: 3488
    e value: 0.0
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    ****ALIGHTMENT****
    sequence: gi|2217338814|ref|XM_011529839.3| PREDICTED: Homo sapiens checkpoint kinase 2 (CHEK2), transcript variant X1, mRNA >gi|2462583655|ref|XM_054325037.1| PREDICTED: Homo sapiens checkpoint kinase 2 (CHEK2), transcript variant X1, mRNA
    length: 3488
    e value: 8.97197e-176
    ATGCTGAGGCTGCTAACGTCTATGGTCGTGATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGC...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    ATGCTGAGGCTGCTAACGTCTATGGTCGTGATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGC...
    ****ALIGHTMENT****
    sequence: gi|2217338814|ref|XM_011529839.3| PREDICTED: Homo sapiens checkpoint kinase 2 (CHEK2), transcript variant X1, mRNA >gi|2462583655|ref|XM_054325037.1| PREDICTED: Homo sapiens checkpoint kinase 2 (CHEK2), transcript variant X1, mRNA
    length: 3488
    e value: 2.94854e-55
    TTTGCCAATCTTGAATGTGTGAATGACAACTACTGGTTTGGGAGGGACAAAAGCTGTGAATATTGCTTTGATGAA...
    ||||| | |||  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    TTTGC-AGTCTAAAATGTGTGAATGACAACTACTGGTTTGGGAGGGACAAAAGCTGTGAATATTGCTTTGATGAA...
    ****ALIGHTMENT****
    sequence: gi|158259742|dbj|AK289360.1| Homo sapiens cDNA FLJ78721 complete cds, highly similar to Homo sapiens CHK2 checkpoint homolog (S. pombe) (CHEK2), transcript variant 1, mRNA
    length: 1837
    e value: 0.0
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    ****ALIGHTMENT****
    sequence: gi|158259742|dbj|AK289360.1| Homo sapiens cDNA FLJ78721 complete cds, highly similar to Homo sapiens CHK2 checkpoint homolog (S. pombe) (CHEK2), transcript variant 1, mRNA
    length: 1837
    e value: 0.0
    GGTCGTGATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGCAGTGCCTGTTCACAGCCCCATGG...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    GGTCGTGATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGCAGTGCCTGTTCACAGCCCCATGG...
    ****ALIGHTMENT****
    sequence: gi|158254943|dbj|AK290754.1| Homo sapiens cDNA FLJ75391 complete cds, highly similar to Homo sapiens CHK2 checkpoint homolog (S. pombe) (CHEK2), transcriptvariant 1, mRNA
    length: 1859
    e value: 0.0
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    ****ALIGHTMENT****
    sequence: gi|158254943|dbj|AK290754.1| Homo sapiens cDNA FLJ75391 complete cds, highly similar to Homo sapiens CHK2 checkpoint homolog (S. pombe) (CHEK2), transcriptvariant 1, mRNA
    length: 1859
    e value: 0.0
    GGTCGTGATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGCAGTGCCTGTTCACAGCCCCATGG...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    GGTCGTGATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGCAGTGCCTGTTCACAGCCCCATGG...
    ****ALIGHTMENT****
    sequence: gi|109451097|emb|CU012979.1| Homo sapiens CHEK2, mRNA (cDNA clone IMAGE:100000140), complete cds, with stop codon, in Gateway system
    length: 1795
    e value: 0.0
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    ****ALIGHTMENT****
    sequence: gi|109451097|emb|CU012979.1| Homo sapiens CHEK2, mRNA (cDNA clone IMAGE:100000140), complete cds, with stop codon, in Gateway system
    length: 1795
    e value: 1.73353e-159
    ATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGCAGTGCCTGTTCACAGCCCCATGGCAGCGTT...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    ATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGCAGTGCCTGTTCACAGCCCCATGGCAGCGTT...
    ****ALIGHTMENT****
    sequence: gi|109451097|emb|CU012979.1| Homo sapiens CHEK2, mRNA (cDNA clone IMAGE:100000140), complete cds, with stop codon, in Gateway system
    length: 1795
    e value: 2.94854e-55
    TTTGCCAATCTTGAATGTGTGAATGACAACTACTGGTTTGGGAGGGACAAAAGCTGTGAATATTGCTTTGATGAA...
    ||||| | |||  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    TTTGC-AGTCTAAAATGTGTGAATGACAACTACTGGTTTGGGAGGGACAAAAGCTGTGAATATTGCTTTGATGAA...
    ****ALIGHTMENT****
    sequence: gi|45356853|gb|AY551305.1| Homo sapiens protein kinase Chk2 transcript variant insZ (CHK2) mRNA, complete cds
    length: 1756
    e value: 0.0
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    ****ALIGHTMENT****
    sequence: gi|45356853|gb|AY551305.1| Homo sapiens protein kinase Chk2 transcript variant insZ (CHK2) mRNA, complete cds
    length: 1756
    e value: 0.0
    ATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGCAGTGCCTGTTCACAGCCCCATGGCAGCGTT...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    ATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGCAGTGCCTGTTCACAGCCCCATGGCAGCGTT...
    ****ALIGHTMENT****
    sequence: gi|45356849|gb|AY551303.1| Homo sapiens protein kinase Chk2 transcript variant del2-3 (CHK2) mRNA, complete cds
    length: 1369
    e value: 0.0
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    ****ALIGHTMENT****
    sequence: gi|45356849|gb|AY551303.1| Homo sapiens protein kinase Chk2 transcript variant del2-3 (CHK2) mRNA, complete cds
    length: 1369
    e value: 6.05061e-159
    ATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGCAGTGCCTGTTCACAGCCCCATGGCAGCGTT...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    ATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGCAGTGCCTGTTCACAGCCCCATGGCAGCGTT...
    ****ALIGHTMENT****
    sequence: gi|45356847|gb|AY551302.1| Homo sapiens protein kinase Chk2 transcript variant del4 (CHK2) mRNA, complete cds
    length: 1551
    e value: 0.0
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    ****ALIGHTMENT****
    sequence: gi|45356847|gb|AY551302.1| Homo sapiens protein kinase Chk2 transcript variant del4 (CHK2) mRNA, complete cds
    length: 1551
    e value: 0.0
    ATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGCAGTGCCTGTTCACAGCCCCATGGCAGCGTT...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    ATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGCAGTGCCTGTTCACAGCCCCATGGCAGCGTT...
    ****ALIGHTMENT****
    sequence: gi|45356845|gb|AY551301.1| Homo sapiens protein kinase Chk2 transcript variant sub3 (CHK2) mRNA, complete cds
    length: 1532
    e value: 0.0
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    ****ALIGHTMENT****
    sequence: gi|45356845|gb|AY551301.1| Homo sapiens protein kinase Chk2 transcript variant sub3 (CHK2) mRNA, complete cds
    length: 1532
    e value: 0.0
    ATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGCAGTGCCTGTTCACAGCCCCATGGCAGCGTT...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    ATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGCAGTGCCTGTTCACAGCCCCATGGCAGCGTT...
    ****ALIGHTMENT****
    sequence: gi|45356837|gb|AY551297.1| Homo sapiens protein kinase Chk2 transcript variant insX (CHK2) mRNA, complete cds
    length: 1771
    e value: 0.0
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    ****ALIGHTMENT****
    sequence: gi|45356837|gb|AY551297.1| Homo sapiens protein kinase Chk2 transcript variant insX (CHK2) mRNA, complete cds
    length: 1771
    e value: 1.73353e-159
    ATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGCAGTGCCTGTTCACAGCCCCATGGCAGCGTT...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    ATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGCAGTGCCTGTTCACAGCCCCATGGCAGCGTT...
    ****ALIGHTMENT****
    sequence: gi|45356837|gb|AY551297.1| Homo sapiens protein kinase Chk2 transcript variant insX (CHK2) mRNA, complete cds
    length: 1771
    e value: 2.94854e-55
    TTTGCCAATCTTGAATGTGTGAATGACAACTACTGGTTTGGGAGGGACAAAAGCTGTGAATATTGCTTTGATGAA...
    ||||| | |||  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    TTTGC-AGTCTAAAATGTGTGAATGACAACTACTGGTTTGGGAGGGACAAAAGCTGTGAATATTGCTTTGATGAA...
    ****ALIGHTMENT****
    sequence: gi|47678366|emb|CR456418.1| Homo sapiens CHEK2 full length open reading frame (ORF) cDNA clone (cDNA clone C22ORF:pGEM.CHEK2)
    length: 1841
    e value: 0.0
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    ****ALIGHTMENT****
    sequence: gi|47678366|emb|CR456418.1| Homo sapiens CHEK2 full length open reading frame (ORF) cDNA clone (cDNA clone C22ORF:pGEM.CHEK2)
    length: 1841
    e value: 2.74697e-163
    GGTCGTGATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGCAGTGCCTGTTCACAGCCCCATGG...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    GGTCGTGATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGCAGTGCCTGTTCACAGCCCCATGG...
    ****ALIGHTMENT****
    sequence: gi|47678366|emb|CR456418.1| Homo sapiens CHEK2 full length open reading frame (ORF) cDNA clone (cDNA clone C22ORF:pGEM.CHEK2)
    length: 1841
    e value: 2.94854e-55
    TTTGCCAATCTTGAATGTGTGAATGACAACTACTGGTTTGGGAGGGACAAAAGCTGTGAATATTGCTTTGATGAA...
    ||||| | |||  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    TTTGC-AGTCTAAAATGTGTGAATGACAACTACTGGTTTGGGAGGGACAAAAGCTGTGAATATTGCTTTGATGAA...
    ****ALIGHTMENT****
    sequence: gi|10441880|gb|AF217975.1|AF217975 Homo sapiens clone PP1425 unknown mRNA
    length: 2592
    e value: 0.0
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    ****ALIGHTMENT****
    sequence: gi|10441880|gb|AF217975.1|AF217975 Homo sapiens clone PP1425 unknown mRNA
    length: 2592
    e value: 0.0
    GGTCGTGATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGCAGTGCCTGTTCACAGCCCCATGG...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    GGTCGTGATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGCAGTGCCTGTTCACAGCCCCATGG...
    ****ALIGHTMENT****
    sequence: gi|38114706|gb|BC004207.2| Homo sapiens CHK2 checkpoint homolog (S. pombe), mRNA (cDNA clone MGC:4081 IMAGE:3532866), complete cds
    length: 1740
    e value: 0.0
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    ****ALIGHTMENT****
    sequence: gi|38114706|gb|BC004207.2| Homo sapiens CHK2 checkpoint homolog (S. pombe), mRNA (cDNA clone MGC:4081 IMAGE:3532866), complete cds
    length: 1740
    e value: 0.0
    GGTCGTGATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGCAGTGCCTGTTCACAGCCCCATGG...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    GGTCGTGATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGCAGTGCCTGTTCACAGCCCCATGG...
    ****ALIGHTMENT****
    sequence: gi|5726656|gb|AF174135.1|AF174135 Homo sapiens protein kinase CHK2 (CHK2) mRNA, complete cds
    length: 1632
    e value: 0.0
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    ****ALIGHTMENT****
    sequence: gi|5726656|gb|AF174135.1|AF174135 Homo sapiens protein kinase CHK2 (CHK2) mRNA, complete cds
    length: 1632
    e value: 0.0
    ATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGCAGTGCCTGTTCACAGCCCCATGGCAGCGTT...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    ATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGCAGTGCCTGTTCACAGCCCCATGGCAGCGTT...
    ****ALIGHTMENT****
    sequence: gi|4007565|emb|AJ131197.1| Homo sapiens mRNA for protein kinase >gi|4206720|gb|AF096279.1| Homo sapiens HuCds1 kinase mRNA, complete cds
    length: 1632
    e value: 0.0
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    ****ALIGHTMENT****
    sequence: gi|4007565|emb|AJ131197.1| Homo sapiens mRNA for protein kinase >gi|4206720|gb|AF096279.1| Homo sapiens HuCds1 kinase mRNA, complete cds
    length: 1632
    e value: 0.0
    ATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGCAGTGCCTGTTCACAGCCCCATGGCAGCGTT...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    ATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGCAGTGCCTGTTCACAGCCCCATGGCAGCGTT...
    ****ALIGHTMENT****
    sequence: gi|3982839|gb|AF086904.1|AF086904 Homo sapiens protein kinase Chk2 (CHK2) mRNA, complete cds
    length: 1735
    e value: 0.0
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    ****ALIGHTMENT****
    sequence: gi|3982839|gb|AF086904.1|AF086904 Homo sapiens protein kinase Chk2 (CHK2) mRNA, complete cds
    length: 1735
    e value: 0.0
    ATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGCAGTGCCTGTTCACAGCCCCATGGCAGCGTT...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    ATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGCAGTGCCTGTTCACAGCCCCATGGCAGCGTT...
    ****ALIGHTMENT****
    sequence: gi|649151526|gb|KJ905414.1| Synthetic construct Homo sapiens clone ccsbBroadEn_14989 CHEK2 gene, encodes complete protein
    length: 1761
    e value: 0.0
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    ****ALIGHTMENT****
    sequence: gi|649151526|gb|KJ905414.1| Synthetic construct Homo sapiens clone ccsbBroadEn_14989 CHEK2 gene, encodes complete protein
    length: 1761
    e value: 0.0
    ATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGCAGTGCCTGTTCACAGCCCCATGGCAGCGTT...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    ATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGCAGTGCCTGTTCACAGCCCCATGGCAGCGTT...
    ****ALIGHTMENT****
    sequence: gi|123999806|gb|DQ896413.2| Synthetic construct Homo sapiens clone IMAGE:100010873; FLH194094.01L; RZPDo839F0569D CHK2 checkpoint homolog (S. pombe) (CHEK2) gene, encodes complete protein
    length: 1672
    e value: 0.0
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    ****ALIGHTMENT****
    sequence: gi|123999806|gb|DQ896413.2| Synthetic construct Homo sapiens clone IMAGE:100010873; FLH194094.01L; RZPDo839F0569D CHK2 checkpoint homolog (S. pombe) (CHEK2) gene, encodes complete protein
    length: 1672
    e value: 0.0
    ATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGCAGTGCCTGTTCACAGCCCCATGGCAGCGTT...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    ATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGCAGTGCCTGTTCACAGCCCCATGGCAGCGTT...
    ****ALIGHTMENT****
    sequence: gi|123992913|gb|DQ893133.2| Synthetic construct clone IMAGE:100005763; FLH194098.01X; RZPDo839F0579D CHK2 checkpoint homolog (S. pombe) (CHEK2) gene, encodes complete protein
    length: 1672
    e value: 0.0
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    ****ALIGHTMENT****
    sequence: gi|123992913|gb|DQ893133.2| Synthetic construct clone IMAGE:100005763; FLH194098.01X; RZPDo839F0579D CHK2 checkpoint homolog (S. pombe) (CHEK2) gene, encodes complete protein
    length: 1672
    e value: 0.0
    ATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGCAGTGCCTGTTCACAGCCCCATGGCAGCGTT...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    ATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGCAGTGCCTGTTCACAGCCCCATGGCAGCGTT...
    ****ALIGHTMENT****
    sequence: gi|306921500|dbj|AB587467.1| Synthetic construct DNA, clone: pF1KB8443, Homo sapiens CHEK2 gene for CHK2 checkpoint homolog, without stop codon, in Flexi system
    length: 1775
    e value: 0.0
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    ****ALIGHTMENT****
    sequence: gi|306921500|dbj|AB587467.1| Synthetic construct DNA, clone: pF1KB8443, Homo sapiens CHEK2 gene for CHK2 checkpoint homolog, without stop codon, in Flexi system
    length: 1775
    e value: 1.73353e-159
    ATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGCAGTGCCTGTTCACAGCCCCATGGCAGCGTT...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    ATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGCAGTGCCTGTTCACAGCCCCATGGCAGCGTT...
    ****ALIGHTMENT****
    sequence: gi|306921500|dbj|AB587467.1| Synthetic construct DNA, clone: pF1KB8443, Homo sapiens CHEK2 gene for CHK2 checkpoint homolog, without stop codon, in Flexi system
    length: 1775
    e value: 2.94854e-55
    TTTGCCAATCTTGAATGTGTGAATGACAACTACTGGTTTGGGAGGGACAAAAGCTGTGAATATTGCTTTGATGAA...
    ||||| | |||  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    TTTGC-AGTCTAAAATGTGTGAATGACAACTACTGGTTTGGGAGGGACAAAAGCTGTGAATATTGCTTTGATGAA...
    ****ALIGHTMENT****
    sequence: gi|109451675|emb|CU013267.1| Homo sapiens CHEK2, mRNA (cDNA clone IMAGE:100000044), complete cds, without stop codon, in Gateway system
    length: 1792
    e value: 0.0
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    ****ALIGHTMENT****
    sequence: gi|109451675|emb|CU013267.1| Homo sapiens CHEK2, mRNA (cDNA clone IMAGE:100000044), complete cds, without stop codon, in Gateway system
    length: 1792
    e value: 1.73353e-159
    ATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGCAGTGCCTGTTCACAGCCCCATGGCAGCGTT...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    ATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGCAGTGCCTGTTCACAGCCCCATGGCAGCGTT...
    ****ALIGHTMENT****
    sequence: gi|109451675|emb|CU013267.1| Homo sapiens CHEK2, mRNA (cDNA clone IMAGE:100000044), complete cds, without stop codon, in Gateway system
    length: 1792
    e value: 2.94854e-55
    TTTGCCAATCTTGAATGTGTGAATGACAACTACTGGTTTGGGAGGGACAAAAGCTGTGAATATTGCTTTGATGAA...
    ||||| | |||  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    TTTGC-AGTCTAAAATGTGTGAATGACAACTACTGGTTTGGGAGGGACAAAAGCTGTGAATATTGCTTTGATGAA...
    ****ALIGHTMENT****
    sequence: gi|33303934|gb|AY335648.1| Synthetic construct Homo sapiens CHK2 checkpoint-like protein (CHEK2) mRNA, partial cds
    length: 1632
    e value: 0.0
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    ****ALIGHTMENT****
    sequence: gi|33303934|gb|AY335648.1| Synthetic construct Homo sapiens CHK2 checkpoint-like protein (CHEK2) mRNA, partial cds
    length: 1632
    e value: 0.0
    ATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGCAGTGCCTGTTCACAGCCCCATGGCAGCGTT...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    ATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGCAGTGCCTGTTCACAGCCCCATGGCAGCGTT...
    ****ALIGHTMENT****
    sequence: gi|61358956|gb|AY888704.1| Synthetic construct Homo sapiens clone FLH019619.01X CHK2 checkpoint-like (CHEK2) mRNA, complete cds
    length: 1632
    e value: 0.0
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    ****ALIGHTMENT****
    sequence: gi|61358956|gb|AY888704.1| Synthetic construct Homo sapiens clone FLH019619.01X CHK2 checkpoint-like (CHEK2) mRNA, complete cds
    length: 1632
    e value: 0.0
    ATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGCAGTGCCTGTTCACAGCCCCATGGCAGCGTT...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    ATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGCAGTGCCTGTTCACAGCCCCATGGCAGCGTT...
    ****ALIGHTMENT****
    sequence: gi|2490761732|ref|XM_034948505.2| PREDICTED: Pan paniscus checkpoint kinase 2 (CHEK2), transcript variant X7, mRNA
    length: 1879
    e value: 0.0
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    ****ALIGHTMENT****
    sequence: gi|2490761732|ref|XM_034948505.2| PREDICTED: Pan paniscus checkpoint kinase 2 (CHEK2), transcript variant X7, mRNA
    length: 1879
    e value: 0.0
    GGTCGTGATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGCAGTGCCTGTTCACAGCCCCATGG...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    GGTCGTGATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGCAGTGCCTGTTCACAGCCCCATGG...
    ****ALIGHTMENT****
    sequence: gi|2490761730|ref|XM_055105905.1| PREDICTED: Pan paniscus checkpoint kinase 2 (CHEK2), transcript variant X6, mRNA
    length: 2628
    e value: 0.0
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    ****ALIGHTMENT****
    sequence: gi|2490761730|ref|XM_055105905.1| PREDICTED: Pan paniscus checkpoint kinase 2 (CHEK2), transcript variant X6, mRNA
    length: 2628
    e value: 0.0
    ATGCTGAGGCTGCTAACGTCTATGGTCGTGATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGC...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    ATGCTGAGGCTGCTAACGTCTATGGTCGTGATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGC...
    ****ALIGHTMENT****
    sequence: gi|2490761728|ref|XM_055105904.1| PREDICTED: Pan paniscus checkpoint kinase 2 (CHEK2), transcript variant X5, mRNA
    length: 1974
    e value: 0.0
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    ****ALIGHTMENT****
    sequence: gi|2490761728|ref|XM_055105904.1| PREDICTED: Pan paniscus checkpoint kinase 2 (CHEK2), transcript variant X5, mRNA
    length: 1974
    e value: 1.16804e-161
    GGTCGTGATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGCAGTGCCTGTTCACAGCCCCATGG...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    GGTCGTGATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGCAGTGCCTGTTCACAGCCCCATGG...
    ****ALIGHTMENT****
    sequence: gi|2490761728|ref|XM_055105904.1| PREDICTED: Pan paniscus checkpoint kinase 2 (CHEK2), transcript variant X5, mRNA
    length: 1974
    e value: 2.94854e-55
    GAATGTGTGAATGACAACTACTGGTTTGGGAGGGACAAAAGCTGTGAATATTGCTTTGATGAACCACTGCTGAAA...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    GAATGTGTGAATGACAACTACTGGTTTGGGAGGGACAAAAGCTGTGAATATTGCTTTGATGAACCACTGCTGAAA...
    ****ALIGHTMENT****
    sequence: gi|2490761726|ref|XM_055105903.1| PREDICTED: Pan paniscus checkpoint kinase 2 (CHEK2), transcript variant X4, mRNA
    length: 2721
    e value: 0.0
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    ****ALIGHTMENT****
    sequence: gi|2490761726|ref|XM_055105903.1| PREDICTED: Pan paniscus checkpoint kinase 2 (CHEK2), transcript variant X4, mRNA
    length: 2721
    e value: 3.81498e-174
    ATGCTGAGGCTGCTAACGTCTATGGTCGTGATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGC...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    ATGCTGAGGCTGCTAACGTCTATGGTCGTGATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGC...
    ****ALIGHTMENT****
    sequence: gi|2490761726|ref|XM_055105903.1| PREDICTED: Pan paniscus checkpoint kinase 2 (CHEK2), transcript variant X4, mRNA
    length: 2721
    e value: 2.94854e-55
    GAATGTGTGAATGACAACTACTGGTTTGGGAGGGACAAAAGCTGTGAATATTGCTTTGATGAACCACTGCTGAAA...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    GAATGTGTGAATGACAACTACTGGTTTGGGAGGGACAAAAGCTGTGAATATTGCTTTGATGAACCACTGCTGAAA...
    ****ALIGHTMENT****
    sequence: gi|2490761725|ref|XM_034948502.2| PREDICTED: Pan paniscus checkpoint kinase 2 (CHEK2), transcript variant X3, mRNA
    length: 1959
    e value: 0.0
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    ****ALIGHTMENT****
    sequence: gi|2490761725|ref|XM_034948502.2| PREDICTED: Pan paniscus checkpoint kinase 2 (CHEK2), transcript variant X3, mRNA
    length: 1959
    e value: 1.16804e-161
    GGTCGTGATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGCAGTGCCTGTTCACAGCCCCATGG...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    GGTCGTGATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGCAGTGCCTGTTCACAGCCCCATGG...
    ****ALIGHTMENT****
    sequence: gi|2490761725|ref|XM_034948502.2| PREDICTED: Pan paniscus checkpoint kinase 2 (CHEK2), transcript variant X3, mRNA
    length: 1959
    e value: 2.94854e-55
    TTTGCCAATCTTGAATGTGTGAATGACAACTACTGGTTTGGGAGGGACAAAAGCTGTGAATATTGCTTTGATGAA...
    ||||| | |||  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    TTTGC-AGTCTAAAATGTGTGAATGACAACTACTGGTTTGGGAGGGACAAAAGCTGTGAATATTGCTTTGATGAA...
    ****ALIGHTMENT****
    sequence: gi|2490761724|ref|XM_034948501.2| PREDICTED: Pan paniscus checkpoint kinase 2 (CHEK2), transcript variant X2, mRNA
    length: 2011
    e value: 0.0
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    ****ALIGHTMENT****
    sequence: gi|2490761724|ref|XM_034948501.2| PREDICTED: Pan paniscus checkpoint kinase 2 (CHEK2), transcript variant X2, mRNA
    length: 2011
    e value: 1.16804e-161
    GGTCGTGATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGCAGTGCCTGTTCACAGCCCCATGG...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    GGTCGTGATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGCAGTGCCTGTTCACAGCCCCATGG...
    ****ALIGHTMENT****
    sequence: gi|2490761724|ref|XM_034948501.2| PREDICTED: Pan paniscus checkpoint kinase 2 (CHEK2), transcript variant X2, mRNA
    length: 2011
    e value: 2.94854e-55
    TTTGCCAATCTTGAATGTGTGAATGACAACTACTGGTTTGGGAGGGACAAAAGCTGTGAATATTGCTTTGATGAA...
    ||||| | |||  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    TTTGC-AGTCTAAAATGTGTGAATGACAACTACTGGTTTGGGAGGGACAAAAGCTGTGAATATTGCTTTGATGAA...
    ****ALIGHTMENT****
    sequence: gi|2490761722|ref|XM_055105902.1| PREDICTED: Pan paniscus checkpoint kinase 2 (CHEK2), transcript variant X1, mRNA
    length: 2757
    e value: 0.0
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    ****ALIGHTMENT****
    sequence: gi|2490761722|ref|XM_055105902.1| PREDICTED: Pan paniscus checkpoint kinase 2 (CHEK2), transcript variant X1, mRNA
    length: 2757
    e value: 3.81498e-174
    ATGCTGAGGCTGCTAACGTCTATGGTCGTGATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGC...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    ATGCTGAGGCTGCTAACGTCTATGGTCGTGATGTCTCGGGAGTCGGATGTTGAGGCTCAGCAGTCTCATGGCAGC...
    ****ALIGHTMENT****
    sequence: gi|2490761722|ref|XM_055105902.1| PREDICTED: Pan paniscus checkpoint kinase 2 (CHEK2), transcript variant X1, mRNA
    length: 2757
    e value: 2.94854e-55
    TTTGCCAATCTTGAATGTGTGAATGACAACTACTGGTTTGGGAGGGACAAAAGCTGTGAATATTGCTTTGATGAA...
    ||||| | |||  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    TTTGC-AGTCTAAAATGTGTGAATGACAACTACTGGTTTGGGAGGGACAAAAGCTGTGAATATTGCTTTGATGAA...
    ****ALIGHTMENT****
    sequence: gi|2468589408|ref|XM_016938890.2| PREDICTED: Pan troglodytes checkpoint kinase 2 (CHEK2), transcript variant X10, mRNA
    length: 1606
    e value: 0.0
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    TGGTGCCTGTGGAGAGGTAAAGCTGGCTTTCGAGAGGAAAACATGTAAGAAAGTAGCCATAAAGATCATCAGCAA...
    ****ALIGHTMENT****
    sequence: gi|2468589408|ref|XM_016938890.2| PREDICTED: Pan troglodytes checkpoint kinase 2 (CHEK2), transcript variant X10, mRNA
    length: 1606
    e value: 3.59205e-54
    GAATGTGTGAATGACAACTACTGGTTTGGGAGGGACAAAAGCTGTGAATATTGCTTTGATGAACCACTGCTGAAA...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    GAATGTGTGAATGACAACTACTGGTTTGGGAGGGACAAAAGCTGTGAATATTGCTTTGATGAACCACTGCTGAAA...



```python

```

# Open CV
## Part 1
```python
# We can import our necessary files
import numpy as np
import matplotlib.pyplot as plt
%matplotlib inline
```


```python
# We can import our CV commands
import cv2
```


```python
# We can import our Mushroom picture
img = cv2.imread("Mushroom.jpg")
```


```python
# We can determine the type of the image
type(img)
```




    numpy.ndarray




```python
# We can load an incorrect image
img_wrong = cv2.imread("wrong/path/doesnot/abcdegh.jpg")
```


```python
# We can check to see if we've downloaded the right image
type(img_wrong)
```




    NoneType




```python
# We can get our image to show
plt.imshow(img)
```




    <matplotlib.image.AxesImage at 0x7f906c8e3550>




![png](output_6_1.png)



```python
# We can transform the color of our image
fix_img = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)
```


```python
# We can show our fixed image
plt.imshow(fix_img)
```




    <matplotlib.image.AxesImage at 0x7f906c8521d0>




![png](output_8_1.png)



```python
# We can make our image gray and see the shape of the image
img_gray = cv2.imread("Mushroom.jpg", cv2.IMREAD_GRAYSCALE)
img_gray.shape
```




    (1080, 1920)




```python
# We can show our gray image
plt.imshow(img_gray)
```




    <matplotlib.image.AxesImage at 0x7f906c836d10>




![png](output_10_1.png)



```python
# We can transform the image again to get a proper gray image
plt.imshow(img_gray, cmap = "gray")
```




    <matplotlib.image.AxesImage at 0x7f906c7a3b50>




![png](output_11_1.png)



```python
# We can resize our image
fix_img.shape
```




    (1080, 1920, 3)




```python
# We can resize our image
new_img = cv2.resize(fix_img,(1000,400))
plt.imshow(new_img)
```




    <matplotlib.image.AxesImage at 0x7f906c714d50>




![png](output_13_1.png)



```python
# We can see our new image dimensions
new_img.shape
```




    (400, 1000, 3)




```python
# We can change the width and height ratios
w_ratio = 0.5
h_ratio = 0.5

new_img = cv2.resize(fix_img, (0,0), fix_img, w_ratio, h_ratio)
```


```python
# We can show our new resized image
plt.imshow(new_img)
```




    <matplotlib.image.AxesImage at 0x7f906c6f5950>




![png](output_16_1.png)



```python
# We can see the shape of our new image
new_img.shape
```




    (540, 960, 3)




```python
# We can flip images
flip_img = cv2.flip(fix_img, 0)
plt.imshow(flip_img)
```




    <matplotlib.image.AxesImage at 0x7f906c65b710>




![png](output_18_1.png)



```python
# We can also flip an image horizontally and vertically
flip_img2 = cv2.flip(fix_img, -1)
plt.imshow(flip_img2)
```




    <matplotlib.image.AxesImage at 0x7f906c5cb4d0>




![png](output_19_1.png)



```python
# We can see the type of our image
type(fix_img)
```




    numpy.ndarray




```python
# We can save the flipped image
cv2.imwrite("Mushroom_fixed_image.jpg", flip_img)
```




    True



## Part 2
```python
import cv2
import matplotlib.pyplot as plt
%matplotlib inline
```


```python
# We can import our image again
img = cv2.imread("Mushroom.jpg")
```


```python
# We can show our original image
plt.imshow(img)
```




    <matplotlib.image.AxesImage at 0x7f906c53d290>




![png](output_24_1.png)



```python
# We can fix our colors to show correctly
img1 = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)
```


```python
# We can show our fixed image
plt.imshow(img1)
```




    <matplotlib.image.AxesImage at 0x7f906c52e390>




![png](output_26_1.png)



```python
# We can convert our colors to something other than RGB
img2 = cv2.cvtColor(img, cv2.COLOR_BGR2HSV)
```


```python
plt.imshow(img2)
```




    <matplotlib.image.AxesImage at 0x7f906c48ef10>




![png](output_28_1.png)



```python
# We can also convert our image to HLS color
img3 = cv2.cvtColor(img, cv2.COLOR_BGR2HLS)
```


```python
# We can show our image in HLS colors
plt.imshow(img3)
```




    <matplotlib.image.AxesImage at 0x7f906c479a50>




![png](output_30_1.png)



```python
# We can assign different images to both img1 and img2 variables
img1 = cv2.imread("DNC.jpg")
img2 = cv2.imread("Mushroom.jpg")
```


```python
# We can show our img1
plt.imshow(img1)
```




    <matplotlib.image.AxesImage at 0x7f906c3e27d0>




![png](output_32_1.png)



```python
# We can convert our image into RGB colors
img1 = cv2.cvtColor(img1, cv2.COLOR_BGR2RGB)
img2 = cv2.cvtColor(img2, cv2.COLOR_BGR2RGB)
```


```python
# We can show the RGB version of img1
plt.imshow(img1)
```




    <matplotlib.image.AxesImage at 0x7f906c34c4d0>




![png](output_34_1.png)



```python
# We can show the RGB version of img2
plt.imshow(img2)
```




    <matplotlib.image.AxesImage at 0x7f906c339210>




![png](output_35_1.png)



```python
# We can resize our images
img1 = cv2.resize(img1, (1200, 1200))
img2 = cv2.resize(img2, (1200, 1200))
```


```python
# We can set our alpha and beta values
alpha = 0.5
beta = 0.5
```


```python
# We can blend our two images together
blended = cv2.addWeighted(img1, alpha, img2, beta, gamma=0)
```


```python
# We can show the blended image
plt.imshow(blended)
```




    <matplotlib.image.AxesImage at 0x7f906c29efd0>




![png](output_39_1.png)



```python
# We can make the pictures more or less transparent
alpha = 0.2
beta = 0.8

blended1 = cv2.addWeighted(img1, alpha, img2, beta, gamma=0)
plt.imshow(blended)
```




    <matplotlib.image.AxesImage at 0x7f906c15aa90>




![png](output_40_1.png)



```python
# We can blend using different image sizes
img1 = cv2.imread("DNC.jpg")
img2 = cv2.imread("Mushroom.jpg")

img1 = cv2.cvtColor(img1, cv2.COLOR_BGR2RGB)
img2 = cv2.cvtColor(img2, cv2.COLOR_BGR2RGB)

img1 = cv2.resize(img1, (600,600))
```


```python
# We can label our resized images as new variables and set them a certain distance from the x and y axis
large_img = img2
small_img =img1

x_offset = 0
y_offset = 0

x_end = x_offset + small_img.shape[1]
y_end = y_offset +small_img.shape[0]

large_img[y_offset:y_end, x_offset:x_end] = small_img
# We can print our newly sized pictures onto each other
plt.imshow(large_img)
```




    <matplotlib.image.AxesImage at 0x7f906c0cb450>




![png](output_42_1.png)


## Part 3
```python
# We can import our necessary programs
import cv2
import matplotlib.pyplot as plt
%matplotlib inline
```


```python
# We can import our rainbow image
img = cv2.imread("Rainbow.jpg")
```


```python
# We can print our rainbow image
plt.imshow(img)
```




    <matplotlib.image.AxesImage at 0x7f906c0b3810>




![png](output_45_1.png)



```python
# We can cancel out background colors
img = cv2.imread("Rainbow.jpg", 0)
```


```python
# We can print our reduced color image
plt.imshow(img, cmap = "gray")
```




    <matplotlib.image.AxesImage at 0x7f905f46c410>




![png](output_47_1.png)



```python
# We can threshold this image
ret1, thresh1 = cv2.threshold(img, 127, 255, cv2.THRESH_BINARY)
```


```python
# We can print our threshold
ret1
```




    127.0




```python
# We can use this threshold to create a black and white image
plt.imshow(thresh1, cmap = "gray")
```




    <matplotlib.image.AxesImage at 0x7f905f3d7090>




![png](output_50_1.png)



```python
# We can inverse the limits of the image
img2 = cv2.imread("Rainbow.jpg", 0)
ret1, thresh1 = cv2.threshold(img2, 127, 255, cv2.THRESH_TRUNC)
plt.imshow(thresh1, cmap = "gray")
```




    <matplotlib.image.AxesImage at 0x7f905f3bc1d0>




![png](output_51_1.png)



```python
# We can make the image appear almost 3D
img3 = cv2.imread("Rainbow.jpg", 0)
ret1, thresh1 = cv2.threshold(img3, 127, 255, cv2.THRESH_TOZERO)
plt.imshow(thresh1, cmap = "gray")
```




    <matplotlib.image.AxesImage at 0x7f905f319d50>




![png](output_52_1.png)



```python
# We can create a new variable for our crossword image
img_r = cv2.imread("Crossword.jpg", 0)
plt.imshow(img_r, cmap = "gray")
```




    <matplotlib.image.AxesImage at 0x7f905f124d90>




![png](output_53_1.png)



```python
# We can create a function to show our images
def show_pic(img):
    fig = plt.figure(figsize = (15,15))
    ax = fig.add_subplot(111)
    ax.imshow(img, cmap = "gray")
```


```python
# We can show our crossword image
show_pic(img_r)
```


![png](output_55_0.png)



```python
# We can use our new "showpic" function to adjust the image quicker
ret1, th1 = cv2.threshold(img_r, 127, 255, cv2.THRESH_BINARY)
show_pic(th1)
```


![png](output_56_0.png)



```python
# We can try adjusting the threshold in order to produce a clearer image
ret1, th1 = cv2.threshold(img_r, 200, 255, cv2.THRESH_BINARY)
show_pic(th1)
```


![png](output_57_0.png)



```python
# We can keep adjusting the threshold to see if we can produce a clearer image
th2 = cv2.adaptiveThreshold(img_r, 255, cv2.ADAPTIVE_THRESH_MEAN_C, cv2.THRESH_BINARY, 11, 8)
```


```python
show_pic(th2)
```


![png](output_59_0.png)



```python
# We can blend the two images to show different thresholds
blended = cv2.addWeighted(src1 = th1, alpha = 0.6,
                         src2 = th2, beta = 0.4, gamma = 0)

show_pic(blended)
```


![png](output_60_0.png)



```python
# We can continue to play with blending the images in order to produce better results
th3 = cv2.adaptiveThreshold(img_r, 255, cv2.ADAPTIVE_THRESH_MEAN_C, cv2.THRESH_BINARY, 11, 8)

blended = cv2.addWeighted(src1 = th1, alpha = 0.6,
                         src2 = th2, beta = 0.4, gamma=0)

show_pic(blended)
```


![png](output_61_0.png)



```python

```

# Aspect Detection
## Corner Detection

```python
# We can import our necessary files
import cv2
import numpy as np
import matplotlib.pyplot as plt
%matplotlib inline
```


```python
# We can import our empty chess board images and convert it to RGB
flat_chess = cv2.imread("Chess Board.jpg")
flat_chess = cv2.cvtColor(flat_chess, cv2.COLOR_BGR2RGB)
plt.imshow(flat_chess)
```




    <matplotlib.image.AxesImage at 0x7f081fc14790>




![png](output_1_1.png)



```python
# We can turn our empty chess board into a gray scale
gray_flat_chess = cv2.cvtColor(flat_chess, cv2.COLOR_BGR2GRAY)
plt.imshow(gray_flat_chess, cmap = "gray")
```




    <matplotlib.image.AxesImage at 0x7f07eaadde50>




![png](output_2_1.png)



```python
# We can import our real chess board image and convert it to RGB color
real_chess = cv2.imread("SetupChessBoard.jpg")
real_chess = cv2.cvtColor(real_chess, cv2.COLOR_BGR2RGB)
plt.imshow(real_chess)
```




    <matplotlib.image.AxesImage at 0x7f07ea25c2d0>




![png](output_3_1.png)



```python
# We can convert our real chess image into a gray scale image
gray_real_chess = cv2.cvtColor(real_chess, cv2.COLOR_BGR2GRAY)
plt.imshow(gray_real_chess, cmap = "gray")
```




    <matplotlib.image.AxesImage at 0x7f07e89cf450>




![png](output_4_1.png)



```python
# We can set up Harris corner detection to detect corners
gray = np.float32(gray_flat_chess)
dst = cv2.cornerHarris(src = gray, blockSize = 2, ksize = 3, k = 0.04)

dst = cv2.dilate(dst, None)
```


```python
# We can assign certain corner detection values to a color and print out the results
flat_chess[dst>0.01*dst.max()] = [255,0,0]

plt.imshow(flat_chess)
```




    <matplotlib.image.AxesImage at 0x7f07e88ca110>




![png](output_6_1.png)



```python
# We can use Harris corner detection on our real chess board
gray = np.float32(gray_real_chess)
dst = cv2.cornerHarris(src = gray, blockSize = 2, ksize = 3, k = 0.04)
dst = cv2.dilate(dst, None)

# We can print out the corner detection results
real_chess[dst>0.01*dst.max()] = [255, 0, 0]

plt.imshow(real_chess)
```




    <matplotlib.image.AxesImage at 0x7f07e8834450>




![png](output_7_1.png)



```python
# We can also use the Shi-Tomasi corner detection method

corners = cv2.goodFeaturesToTrack(gray_flat_chess, 64, 0.01, 10)
```


```python
# We can assign the color red to the Shi-Tomasi method and print the results
corners = np.int0(corners)

for i in corners:
    x,y = i.ravel()
    cv2.circle(flat_chess, (x,y),3,(255,0,0), -1)
    
plt.imshow(flat_chess)
```




    <matplotlib.image.AxesImage at 0x7f07e8826710>




![png](output_9_1.png)



```python
# We can use the Shi-Tomasi method on the real chess board and use a different color for the corner detection
corners = cv2.goodFeaturesToTrack(gray_real_chess, 100, 0.01, 10)

corners = np.int0(corners)

for i in corners:
    x, y = i.ravel()
    cv2.circle(real_chess, (x,y), 3, (0,255,0), -1)
    
plt.imshow(real_chess)
```




    <matplotlib.image.AxesImage at 0x7f07e87908d0>




![png](output_10_1.png)



```python

```

## Edge Detection
```python
# We can import our functions
import cv2
import numpy as np
import matplotlib.pyplot as plt
%matplotlib inline
```


```python
# We can import and show our image
img = cv2.imread("Mushroom.jpg")
plt.imshow(img)
```




    <matplotlib.image.AxesImage at 0x7fae42a55490>




![png](output_1_1.png)



```python
# We can create a function to detect edges
edges = cv2.Canny(image = img, threshold1 = 127, threshold2 = 127)

plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7fae4298aa50>




![png](output_2_1.png)



```python
# We can find the median color value of the image
med_value = np.median(img)

med_value
```




    104.0




```python
# We can change our threshold limits
lower = int(max(0, 0.7*med_value))
upper = int(min(255, 1.3*med_value))

edges = cv2.Canny(img, threshold1 = lower, threshold2 = upper)

plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7fae4288fa10>




![png](output_4_1.png)



```python
# We can keep adjusting the threshold value to improve the clarity
edges = cv2.Canny(image = img, threshold1 = lower, threshold2 = upper + 100)

plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7fae427fb390>




![png](output_5_1.png)



```python
# We can try to blur the picture to increase edge detection
blurred_img = cv2.blur(img, ksize = (5,5))

edges = cv2.Canny(image=blurred_img,
                 threshold1 = lower,
                 threshold2 = upper)

plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7fae4276e610>




![png](output_6_1.png)



```python
# We can increase the kernel size of the pixels
blurred_img = cv2.blur(img, ksize = (7,7))

edges = cv2.Canny(image=blurred_img,
                 threshold1 = lower,
                 threshold2 = upper)

plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7fae4275b410>




![png](output_7_1.png)



```python
# We can increase the upper threshold limit
blurred_img = cv2.blur(img, ksize = (7,7))

edges = cv2.Canny(image=blurred_img,
                 threshold1 = lower,
                 threshold2 = upper +50)

plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7fae400b3210>




![png](output_8_1.png)



```python
# We can increase the upper threshold limit even more
blurred_img = cv2.blur(img, ksize = (7,7))

edges = cv2.Canny(image=blurred_img,
                 threshold1 = lower,
                 threshold2 = upper +100)

plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7fae40011f50>




![png](output_9_1.png)



```python
# We can change the upper threshold limit
blurred_img = cv2.blur(img, ksize = (7,7))

edges = cv2.Canny(image=blurred_img,
                 threshold1 = lower,
                 threshold2 = upper +60)

plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7fae3a7edcd0>




![png](output_10_1.png)



```python
# We can increase the upper threshold limit
blurred_img = cv2.blur(img, ksize = (8,8))

edges = cv2.Canny(image=blurred_img,
                 threshold1 = lower,
                 threshold2 = upper +50)

plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7fae3a75c750>




![png](output_11_1.png)



```python

```

# Feature Detection
## Feature Matches

```python
# We can import our necessary files
import cv2
import numpy as np
import matplotlib.pyplot as plt
%matplotlib inline
```


```python
# We can set up our display function
def display(img,cmap="gray"):
    fig = plt.figure(figsize = (12,10))
    ax = fig.add_subplot(111)
    ax.imshow(img,cmap="gray")
```


```python
# We can upload our Cap'N crunch image
Capn_crunch = cv2.imread("Cap'N Crunch.jpg", 0)
display(Capn_crunch)
```


![png](output_2_0.png)



```python
# We can upload our cereal box image
cereals = cv2.imread("Cereal Boxes.jpg", 0)
display(cereals)
```


![png](output_3_0.png)



```python
# We can use Brute Force matching to find a match for our Cap'N crunch
orb = cv2.ORB_create()

kp1,des1 = orb.detectAndCompute(Capn_crunch, mask=None)
kp2,des2 = orb.detectAndCompute(cereals, mask=None)
```


```python
# We can run the results of Brute force matching
bf = cv2.BFMatcher(cv2.NORM_HAMMING, crossCheck = True)
matches= bf.match(des1, des2)
```


```python
# We can sort the matches of our brute force matching
matches = sorted(matches, key = lambda x:x.distance)
```


```python
# We can calculate the first 25 matches
Capn_crunch_matches = cv2.drawMatches(Capn_crunch, kp1, cereals, kp2, matches[:25], None, flags = 2)
```


```python
# We can display the matches
display(Capn_crunch_matches)
```


![png](output_8_0.png)



```python
# We can use the sift technique to match
sift = cv2.SIFT_create()
```


```python
# We can continue setting up the sift technique
kp1, des1 = sift.detectAndCompute(Capn_crunch, None)
kp2, des2 = sift.detectAndCompute(cereals, None)
```


```python
# We can use brute force along with the sift method
bf = cv2.BFMatcher()
matches = bf.knnMatch(des1, des2, k=2)
```


```python
# We can determine is a match is good or bad
good = []

for match1, match2 in matches:
    if match1.distance < 0.75*match2.distance:
        good.append([match1])
```


```python
# We can print the length of total matches and good matches
print("length of total matches:", len(matches))
print("length of good matches:", len(good))
```

    length of total matches: 13686
    length of good matches: 517



```python
# We can print out the sift matches
sift_matches = cv2.drawMatchesKnn(Capn_crunch, kp1, cereals, kp2, good, None, flags =2)
display(sift_matches)
```


![png](output_14_0.png)



```python
# We can set up sift based matching again
sift = cv2.SIFT_create()

kp1, des1 = sift.detectAndCompute(Capn_crunch, None)
kp2, des2 = sift.detectAndCompute(cereals, None)
```


```python
# We can set up a new type of matching called Flanbased matching
flann_index_KDtree = 0
index_params = dict(algorithm=flann_index_KDtree, trees = 5)
search_params = dict(checks=50)
```


```python
# We can run our flan matcher
flann = cv2.FlannBasedMatcher(index_params, search_params)

matches = flann.knnMatch(des1, des2, k=2)

good = []

for match1, match2, in matches:
    if match1.distance < 0.75*match2.distance:
        good.append([match1])
```


```python
# We can print our results for the flan matcher
flann_matches = cv2.drawMatchesKnn(Capn_crunch, kp1, cereals, kp2, good, None, flags=0)
display(flann_matches)
```


![png](output_18_0.png)



```python
# We can set up our sift match again
sift = cv2.SIFT_create()

kp1, des1 = sift.detectAndCompute(Capn_crunch, None)
kp2, des2 = sift.detectAndCompute(cereals, None)
```


```python
# We can also set up our flan index
flann_index_KDtree = 0
index_params = dict(algorithm= flann_index_KDtree, trees = 50)
search_params = dict(checks = 50)
```


```python
# We can continue setting up our flann parameters
flann = cv2.FlannBasedMatcher(index_params, search_params)

matches = flann.knnMatch(des1, des2, k = 2)
```


```python
# We can create a mask for our matches
matchesMask = [[0,0] for i in range(len(matches))]
```


```python
# We can reduce the distance between our features for better matches
for i, (match1, match2) in enumerate(matches):
    if match1.distance <0.75*match2.distance:
        matchesMask[i] = [1,0]

draw_params = dict(matchColor = (0,255,0),
                  singlePointColor = (255, 0,0),
                  matchesMask = matchesMask,
                  flags = 0)
```


```python
# We can display our flann matches with masks
flann_matches = cv2.drawMatchesKnn(Capn_crunch, kp1, cereals, kp2, matches, None, **draw_params)

display(flann_matches)
```


![png](output_24_0.png)



```python

```

## Object Detection

```python
# We can import our necessary functions
import cv2
import numpy as np
import matplotlib.pyplot as plt
%matplotlib inline
```


```python
# We can import our training sunflower in
full = cv2.imread("Training Sunflower.jpg")
```


```python
# We can change the color of our sunflower to the RGB color scheme
full = cv2.cvtColor(full, cv2.COLOR_BGR2RGB)
```


```python
# We can show our sunflower
plt.imshow(full)
```




    <matplotlib.image.AxesImage at 0x7f9d82cc9910>




![png](output_3_1.png)



```python
# We can also upload our testing sunflower
test = cv2.imread("Testing Sunflower.jpg")
```


```python
# We can change the color of our sunflower to the RGB color scheme
test = cv2.cvtColor(test, cv2.COLOR_BGR2RGB)
```


```python
# We can show our testing sunflowers
plt.imshow(test)
```




    <matplotlib.image.AxesImage at 0x7f9d81d22c90>




![png](output_6_1.png)



```python
# We can print our image shapes
print("Test image shape:", full.shape)
print("Training image shape:", test.shape)
```

    Test image shape: (1555, 1404, 3)
    Training image shape: (675, 1200, 3)



```python
# We can use different methods to detect objects so we can save them to a variable
methods = ["cv2.TM_CCOEFF", "cv2.TM_CCOEFF_NORMED", "cv2.TM_CCORR_NORMED", "cv2.TM_SQDIFF", "cv2.TM_SQDIFF_NORMED"]
```


```python
# We can write a loop to use the different methods and create a copy of the image
for m in methods:
    
    test_copy = test.copy()
    method = eval(m)
    
    res = cv2.matchTemplate(test_copy, full, method)
    
    min_val, max_val, min_loc, max_loc = cv2.minMaxLoc(res)
    
#We can draw a rectangle onto our picture
    if method in [cv2.TM_SQDIFF, cv2.TM_SQDIFF_NORMED]:
        top_left = min_loc
    else:
        top_left = max_loc
# We want to add width and height to a certain point    
    height, width, channels = full.shape
    bottom_right = (top_left[0] + width, top_left[1] + height)
    
    cv2.rectangle(test_copy, top_left, bottom_right, (255, 0, 0),10)
    
# We can plot aspects of our matching pictures
    plt.subplot(121)
    plt.imshow(res)
    plt.title("Heatmap of template matching")
    plt.subplot(122)
    plt.imshow(test_copy)
    plt.title("Detection of template")
    
# We can add a title for the method used
    plt.suptitle(m)
    
    plt.show()
    print("\n")
    print("\n")
```


![png](output_9_0.png)


    
    
    
    



![png](output_9_2.png)


    
    
    
    



![png](output_9_4.png)


    
    
    
    



![png](output_9_6.png)


    
    
    
    



![png](output_9_8.png)


    
    
    
    



```python

```

