# ComPotts Suite

One Paragraph of project description goes here


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

which will automatically download them, install them and put them into your home.


### If you want to use VizContacts
To visualize predicted contacts using VizContacts, you need to install :
* PyMOL : https://pymol.org/
* Circos : http://www.circos.ca/

## Installation

### Compile the C++ solver library
Go to apurva_compotts and run

```
bash compile.bash
```

### Install Python modules
run

```
python3 setup.py install
```

## Running the tests

```
cd vapotts
python3 tests
```


## Getting started

### ComFeature
ComPotts inputs two folders which both have to contain the following files (with the corresponding file name) :
* sequence.fasta : the sequence in fasta format
* original_aln.fasta : the alignment (a priori) derived from the sequence 
* train_aln.fasta : the alignment which was used to train the Potts model, which may be different from the original alignment due to filters and trimmings
* potts_model.mrf : the Potts model to be aligned, in a binary msgpack format as outputed by CCMpredPy
* mrf_pos_to_aln_pos.csv : a csv file which contains a mapping from positions in the Potts model to positions in the original alignment
* mrf_pos_to_seq_pos.csv : a csv file which contains a mapping from positions in the Potts model to positions in the sequence

This folder is created by invoking ComFeature with the desired arguments.

### ComPotts
Use :
```
compotts -f1 feature_folder_1/ -f2 feature_folder_2/ -o output_folder/
```

or, if you are only interested in the aligned positions of Potts models :
```
compotts -p1 potts_model_1.mrf -p2 potts_model_2.mrf -o output_folder/
```

output_folder/ contains at least two files:
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

