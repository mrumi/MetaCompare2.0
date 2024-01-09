# MetaCompare2.0

MetaCompare2.0 is an updated version of MetaCompare. 

## Getting Started
### System Requirements 
**Source Installation**

* git installed
* Python3 with `numpy`, `pandas` and `biopython` package installed
  * `pip install numpy`
  * `pip install pandas`
  * `pip install biopython`
  * `pip install pprodigal`

* DIAMOND installed (https://github.com/bbuchfink/diamond)
* MMseqs2 installed (https://github.com/soedinglab/MMseqs2)

**Conda Installation**

Required Python packages and alignment tools can also be installed using conda.

Python packages installation 

  * `conda install -c anaconda numpy`
  * `conda install -c anaconda pandas`
  * `conda install -c bioconda biopython`
  * `conda install -c bioconda pprodigal`
    
 Alignment tools installation 

* `conda install -c bioconda diamond=0.9.14`
* `conda install -c bioconda mmseqs2`

### Installing

**Step 1:** Change the current working directory to the location where you want the cloned `MetaCompare2.0` directory to be made.

**Step 2:** Clone the repository using git command
```
~$ git clone https://github.com/mrumi/MetaCompare2.0.git
```

**Step 3:** Download the compressed Database file using the following link and move it to the same directory as the code files.

Download link: https://drive.google.com/file/d/10Rc-Fc5gALUHZs4Yd0NUTfKu--C0fJNI/view?usp=sharing

Uncompress the database before using. 

```
~/MetaCompare2.0$ tar -zxvf metacmpDB.tar.gz
```

## Running MetaCompare2.0

### Preparing input files

MetaCompare2.0 requires one FASTA file as input and that is assembled contigs generated by a standard assembler. 

### Running

Suppose you have an assembled contigs file, `S1.fa` (*This file is already in your working directory*).

Use following command to run MetaCompare2.0.
```
~$ python metacompare.py -c S1.fa
```
The output should look like as follows:
```
Running Prodigal
Running Diamond Blastx on ARGDB
Running Diamond Blastx on ARGDB_hh
Running Diamond Blastx on MGEDB
Running mmseq2 on GTDB
```
The final output will be found in a text file ending with "_out.txt". 

You can see detailed description for command line options by using `-h` option.
```
~$ python metacompare.py -h
```

### Citation
Rumi, M. A., Oh, M., Davis, B. C., Juvekar, A. Brown, C. L., Vikesland, P. J., Pruden, A., & Zhang, L. (2024). MetaCompare 2.0: Differential ranking of ecological and human health resistome risks.

