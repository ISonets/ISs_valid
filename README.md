# ISs_valid
ISs_valid.py is a script to filter results of ISEScan and BLAST for rRNA, ISs and phage genes by length and coverage. This script is used in connection with Snakemake file, and can be used as a stand-alone script.
## Dependencies
  - pandas
  - numpy
## Usage
```python ISs_valid.py {path_to_assembly_depth_file} {path_to_ISEscan_res_csv}```
## How it works
It reads ```.depth``` file of assembly, stores its pathway for finding and reading BLAST results files in ```.txt```, as also ISEScan results ```.csv```. 
It merges them into 1 dataframe, and for each row in it checks, if any determinant (IS/phage protein/rRNA) found in each specific contig is **EITHER** over-covered (1.25x fold or more(default 1.25)) 
**OR** if its length exceeds 50% (default 0.5) or more of contig length, and stores its results in separate ```.csv```.

Also this script checks if median of each contig is bigger(1.5x times or more) than median of the whole assembly, and stores its results in another separate ```.csv.```
## Future plans
  - use argparse instead of sys to adjust length and coverage thresholds on the go;
  - docstring and help messages
## Authors
I.Sonets

D.Konanov
