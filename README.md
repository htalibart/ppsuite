# ComPotts Suite

This package contains :
* ComPotts, a tool to align potts models and corresponding sequences / alignments
* ComFeature inputs your sequence / alignments, processes them (filtering, trimming...), infers a Potts model, and outputs a folder with all the features that ComPotts needs to align the models and the corresponding sequences / alignments
* VizPotts allows you to visualize inferred Potts models
* VizContacts allows you to visualize the top N couplings of the inferred Potts model and whether the couplings are contacts in the 3D structure or not. You can visualize them around a circle thanks to Circos, or on a PDB structure using PyMOL.


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
* kneebow

and the following tools which you have to install and add to your path:
* CCMpredPy : https://github.com/susannvorberg/CCmpredPy
* HH-suite : https://github.com/soedinglab/hh-suite
* trimal : https://github.com/scapella/trimal

or you can run install_required_tools.sh :
```
bash install_required_tools.sh
```

which will automatically download them into your home ~/, install them and add them to your path.


### If you want to use VizContacts
To visualize predicted contacts using VizContacts, you also need to install :
* PyMOL : https://pymol.org/
* Circos : http://www.circos.ca/

## Installation

### Compile the C++ solver library

```
cd apurva_compotts/
bash compile.bash
```

### Install Python modules

```
python3 setup.py install
```

## Running the tests

```
python3 tests
```


## Getting started

### ComFeature
ComPotts inputs two folders which both have to contain the following files (with the corresponding file names) :
* sequence.fasta : the sequence in fasta format
* original_aln.fasta : the alignment (a priori) derived from the sequence 
* train_aln.fasta : the alignment which was used to train the Potts model, which may be different from the original alignment due to filters and trimmings
* potts_model.mrf : the Potts model to be aligned, in a binary msgpack format as outputed by CCMpredPy
* mrf_pos_to_aln_pos.csv : a csv file which contains a mapping from positions in the Potts model to positions in the original alignment
* mrf_pos_to_seq_pos.csv : a csv file which contains a mapping from positions in the Potts model to positions in the sequence

This folder is created by invoking ComFeature with the desired arguments.


#### Review of the main use cases

* You have a sequence file, and a big alignment that you retrieved with HHblits for example. You want to filter it to remove redundancies, to trim it to remove positions with too many gaps, and you want to use only the first 1000 sequences to train the Potts model.
```
comfeature -f output_feature_folder/ -s path/to/your/sequence.fasta -aln path/to/your/alignment.fasta --hhfilter_threshold 80 --trimal_gt 0.8 -maxnb 1000
```
By default, comfeature filters the alignment with a 80% threshold using HHfilter, trims it with trimal using options -gt 0.8 -cons 0 and takes the first 1000 sequences but you can change these parameters or choose not to do any of this with options such as --dont_trim_alignment.

* You have a sequence file and you want ComFeature to call HHblits to get more sequences:
```
comfeature -f output_feature_folder/ -s path/to/your/sequence.fasta -fetch -d path/to/the/database
```

* You have a sequence file and you want to infer a Potts model only from the sequence
```
comfeature -f output_feature_folder/ -s path/to/your/sequence.fasta --inference_type one_submat
```

or
```
comfeature -f output_feature_folder/ -s path/to/your/sequence.fasta --inference_type one_hot
```

* You have a sequence file, you used BLAST to retrieve sequences (in an unaligned fasta file) and you want to use as many sequences as there are before the E-value elbow :
```
comfeature -f output_feature_folder/ -s path/to/your/sequence.fasta -ualn path/to/your/unaligned_sequences.fasta --use_evalue_cutoff --blast_xml path/to/blast_output.xml
```

* You have an alignment - and no imagination for the name of your output folder <br/>

ComFeature will use the first sequence of your alignment as the reference sequence.
```
comfeature -aln path/to/your/alignment.fasta
```
The feature folder will be created in the current folder, and will be named using the uuuid generator.


### ComPotts

Use :
```
compotts -f1 feature_folder_1/ -f2 feature_folder_2/ -o output_folder/
```

or, if you are only interested in the aligned positions of Potts models :
```
compotts -p1 potts_model_1.mrf -p2 potts_model_2.mrf -o output_folder/
```

output_folder/ will contain at least two files:
* aln.csv which contains the resulting aligned positions
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
* visualizing one specific w coupling matrix of one Potts model:
```
vizpotts -p potts_model.mrf -i 15 -j 23 
```
You can specify whether for you positions start at 0 or at 1 with option --start_at_0 or --start_at_1. 


### VizContacts (under development)
VizContacts allows you to visualize contacts predicted by CCMpredPy :
* around a circle thanks to Circos, using VizCircos
* on a PDB structure thanks to PyMOL
Both input a feature folder (provided by ComFeature). Blue indicates that the coupling is a "true contact" and red indicates that the positions in the coupling are not in contact in the PDB structure.

#### VizPyMOL

```
vizpymol -f 5lqp_feature_folder/ -i 5lqp --chain_id AB
```


#### VizCircos
```
vizpymol -f 5lqp_feature_folder/ -i 5lqp --chain_id AB
```

## License

[TODO]

## Acknowledgments

[TODO]

