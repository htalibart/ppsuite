[TODO ajouter les versions]

# ComPotts Suite

This package contains the following tools :

* *ComPotts* aligns Potts models and corresponding sequences / alignments
* *MakePotts* builds a Potts model from a sequence or a MSA (output is a folder which also contains positions matching between model and MSA / sequence files)
* *VizPotts* allows you to visualize inferred Potts models
* *VizContacts* allows you to visualize the top N couplings of the inferred Potts model and whether the couplings are contacts in the 3D structure or not. You can visualize them around a circle thanks to Circos, or on a PDB structure using PyMOL.

## Requirements

### Necessary requirements
ComPotts was developped with Python3.6 and requires the following packages (which will normally be automatically installed by setup.py):

* numpy
* pandas
* biopython
* msgpack
* scipy
* matplotlib
* seaborn
* scikit-learn
* kneebow

and the following tools which you have to install and add to your path:

* CCMpredPy : https://github.com/soedinglab/ccmgen
```
git clone https://github.com/soedinglab/ccmgen
cd CCMgen
python3 setup.py install --user
```

* HH-suite : https://github.com/soedinglab/hh-suite
```
git clone https://github.com/soedinglab/hh-suite.git
mkdir -p hh-suite/build && cd hh-suite/build
cmake -DCMAKE_INSTALL_PREFIX=. ..
make -j 4 && make install
```
and export to your path :
```
export PATH="$(pwd)/bin:$(pwd)/scripts:$PATH"
```

* trimal : https://github.com/scapella/trimal
```
git clone https://github.com/scapella/trimal
cd trimal/source
make
```
and add trimal/source to your path 
```
export PATH="$(pwd):$PATH"
```

### If you want to use BLAST instead of HHblits to retrieve sequences
You need to have BLAST installed on your machine. Instructions are available at https://www.ncbi.nlm.nih.gov/books/NBK279671/

### If you want to use VizContacts
To visualize predicted contacts using VizContacts, you also need to install :

* PyMOL : https://pymol.org/
* Circos : http://www.circos.ca/

## Installation

### Compile the C++ solver library

```
cd apurva_compotts/
bash compile.bash
cd ..
```

### Install Python modules

```
python3 setup.py install --user
```

## Running the tests

```
python3 tests
```


## Getting started

### Use case : build a Potts model for 1CC8 (and visualize it)
Let's build a Potts model for [1CC8](https://www.rcsb.org/structure/1cc8) based on an alignment of its close homologs. We provided an alignment in tests/examples/1cc8_everything/1cc8.a3m (see section [TODO] if you want to know how to get an alignment starting from a sequence). We call MakePotts and specify where to find the sequence (in fasta format) and the alignment (in fasta or a3m format).

```
makepotts -f 1cc8_potts_folder/ -s tests/examples/1cc8_everything/1cc8.fasta -aln tests/examples/1cc8_everything/1cc8.a3m
```

This command creates a folder named 1cc8_potts_folder/ which contains the Potts model which was inferred by CCMpredPy - with carefully selected arguments - from the MSA after having processed it (i.e. trimmed it and reduced the number of sequences). The folder contains this training MSA and information on the matching between positions of the training MSA and the original MSA.

You can use it as a ComPotts input (see section [TODO lien]) to align it with another Potts model.

Optionally, you can also use the visualization tools to see what it looks like.

You can visualize the Potts model with VizPotts
```
vizpotts -f 1cc8_potts_folder/
```

<img src="tests/examples/output_examples/1cc8_potts_model.png" width="800">
[TODO rappeler les formules ?]
The first graph shows the inferred fields parameters : each column at position i is a vector vi where the higher vi(a) is, the more frequent the letter a is at position i.

The second graph shows the inferred coupling parameters : each point is the Frobenius norm of the coupling matrix wij between positions i and j. The higher the norm, the more likely there is an interaction between positions i and j (and this interaction may be a contact).

The third graph shows the Euclidean norms of the vectors vi. A higher norm suggests a more conserved position.


You can also, for example, visualize one specific coupling parameter
```
vizpotts -f 1cc8_potts_folder/ -i 18 -j 22
```

<img src="tests/examples/output_examples/1cc8_w_18_22.png" width="500">

You can visualize the top 25 predicted contacts on a 3D structure using PyMOL
```
vizpymol -f 1cc8_potts_folder/ -i 1cc8 --chain_id A --top 25
```
<img src="tests/examples/output_examples/1cc8_pymol.png" width="500">

This command displays the top 25 contacts predicted by CCMpredPy, that is to say the couplings with the highest Frobenius norm (after APC correction). A green line indicates a "true" contact (i.e. the distance between the positions in the 3D strucutre is less than 8 Angström) while a red line indicates a "false" contact.

You can visualize the same information on a Circos
```
vizcircos -f 1cc8_potts_folder/ -i 1cc8 --chain_id A --top 25
```
<img src="tests/examples/output_examples/1cc8_circos.png" width="500">



### Use case : align 1CC8 and 4YDX

Let's align 1CC8 with one of its homologs, 4YDX. Potts models were built with MakePotts in folders tests/examples/1cc8_potts_folder and tests/examples/4ydx_potts_folder.

```
compotts -f1 tests/examples/1cc8_potts_folder -f2 tests/examples/4ydx_potts_folder -o 1cc8_4ydx_compotts_output/ -oaln -osaln
```

1cc8_4ydx_compotts_output/ now contains :

* aln.csv : contains the list of aligned positions in the two Potts models
* aln_original.csv : contains the list of aligned positions of the original alignments (before trimming etc.)
* aln_sequences.csv : contains the list of aligned positions of the original sequences
* info.csv : contains information about the alignment made by ComPotts (including similarity score and computing time)
* aligned_training_sets.fasta (thanks to -oaln option) : a MSA made with the aligned columns of the two train MSAs
* aligned_sequences.fasta (thanks to -osaln option) : a fasta file with an alignment of the two original sequences derived from the alignment made by ComPotts. Since alignment was trimmed and not all positions are aligned by ComPotts, it is not possible to have an alignment of the whole sequences. Ambiguous regions are marked with an X.


You can visualize the vi parameters of the two Potts models at aligned positions using VizPotts :
```
vizpotts -f tests/examples/1cc8_potts_folder tests/examples/4ydx_potts_folder -aln 1cc8_4ydx_compotts_output/aln.csv
```

<img src="tests/examples/output_examples/1cc8_4ydx_v.png" width="500">

### More detailed presentation of the packages
ComPotts aligns two Potts models built using MakePotts. Its inputs are two folders which contain the following files :

* sequence.fasta : the sequence in fasta format
* original_aln.fasta : the alignment (a priori) derived from the sequence 
* train_aln.fasta : the alignment which was used to train the Potts model, which may be different from the original alignment due to filters and trimmings
* potts_model.mrf : the Potts model to be aligned, in a binary msgpack format as outputed by CCMpredPy
* mrf_pos_to_aln_pos.csv : a csv file which contains a mapping from positions in the Potts model to positions in the original alignment
* mrf_pos_to_seq_pos.csv : a csv file which contains a mapping from positions in the Potts model to positions in the sequence

This folder is created by invoking MakePotts with the desired arguments.


### MakePotts

#### Some use cases

* You have a sequence file, and a big alignment that you retrieved with HHblits for example. You want to filter it to remove redundancies, to trim it to remove positions with too many gaps, and you want to use only the first 1000 sequences to train the Potts model.
```
makepotts -f output_feature_folder/ -s path/to/your/sequence.fasta -aln path/to/your/alignment.fasta --hhfilter_threshold 80 --trimal_gt 0.8 -maxnb 1000
```
By default, makepotts filters the alignment with a 80% threshold using HHfilter, trims it with trimal using options -gt 0.8 -cons 0 and takes the first 1000 sequences but you can change these parameters or choose not to do any of this with options such as --dont_trim_alignment.

* You have a sequence file and you want MakePotts to call HHblits to get more sequences:
```
makepotts -f output_feature_folder/ -s path/to/your/sequence.fasta -fetch -d path/to/the/database
```

* You have a sequence file and you want to infer a Potts model only from the sequence
```
makepotts -f output_feature_folder/ -s path/to/your/sequence.fasta --inference_type one_submat
```

or
```
makepotts -f output_feature_folder/ -s path/to/your/sequence.fasta --inference_type one_hot
```

* You have a sequence file, you used BLAST to retrieve sequences (in an unaligned fasta file) and you want to use as many sequences as there are before the E-value elbow :
```
makepotts -f output_feature_folder/ -s path/to/your/sequence.fasta -ualn path/to/your/unaligned_sequences.fasta --use_evalue_cutoff --blast_xml path/to/blast_output.xml
```

* You have an alignment - and no imagination for the name of your output folder <br/>

MakePotts will use the first sequence of your alignment as the reference sequence.
```
makepotts -aln path/to/your/alignment.fasta
```
The feature folder will be created in the current folder, and will be named using the uuuid generator.


### ComPotts

Use :
```
compotts -f1 potts_folder_1/ -f2 potts_folder_2/ -o output_folder/
```

or, if you are only interested in the aligned positions of Potts models files :
```
compotts -p1 potts_model_1.mrf -p2 potts_model_2.mrf -o output_folder/
```

output_folder/ will contain at least :

* aln.csv which contains the resulting aligned positions
* aln_original.csv which contains the aligned positions in the original alignment (before trimming)
* aln_sequences.csv which contains the aligned positions in the input sequences.
* info.csv which contains information about the alignment (especially the similarity score)

and depending on the options you choose, it can also contain :

* an alignment of the train alignments in fasta format (option -oaln)
* an alignment of the input sequences in fasta format (option -osaln)


### VizPotts
VizPotts provides functions to visualize Potts models and their alignment. Though a number of functions are implemented and can be used by importing vizpotts as a Python module, the only command-line options available at the moment are :

* visualizing the parameters of one or several Potts models :
```
vizpotts -p potts_model_1.mrf potts_model_2.mrf potts_model_3.mrf
```

or
```
vizpotts -p potts_model_folder_1 potts_model_folder_2 potts_model_folder_3
```



* visualizing one specific w coupling matrix of one Potts model:
```
vizpotts -p potts_model.mrf -i 15 -j 23 
```
You can specify whether for you positions start at 0 or at 1 with option --start_at_0 or --start_at_1. 

* visualizing the fields parameters v of two Potts models aligned by ComPotts :
```
vizpotts -f potts_model_folder_1 potts_model_folder_2 -aln compotts_aln_output/aln.csv
```


### VizContacts

VizContacts allows you to visualize contacts predicted by CCMpredPy :
* around a circle thanks to Circos, using VizCircos
* on a PDB structure thanks to PyMOL
Both input a feature folder (provided by MakePotts). Green indicates that the coupling is a "true contact" and red indicates that the positions in the coupling are not in contact in the PDB structure.

#### VizPyMOL

```
vizpymol -f 5lqp_feature_folder/ -i 5lqp --chain_id AB
```


#### VizCircos
```
vizcircos -f 5lqp_feature_folder/ -i 5lqp --chain_id AB
```




## License

[TODO]

## Acknowledgments

[TODO]

