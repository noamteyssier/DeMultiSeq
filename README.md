# DeMultiSeq

## Installation
```bash
# run installation script (will install cargo if not found in path)
# will then compile executable and symlink to bin/
./setup.sh
```

## Usage
```bash

# run help menu
./bin/demux -h

# example usage
./bin/demux \
  -i <seq1.fq.gz> \
  -I <seq2.fq.gz> \
  -c <cell barcodes> \
  -m <multiseq barcodes>


# example usage with tolerance of 1 hamming distance
./bin/demux \
  -i <seq1.fq.gz> \
  -I <seq2.fq.gz> \
  -c <cell barcodes> \
  -m <multiseq barcodes>
  -t 1
```
