Input file format:
```
<cut cost> <join cost> <insertion cost> <deletion cost>
<bracket sequence that represents tree>
<leaf id> <number of structures>
structure
...
structure
...
```

Usage:
```
./reconstruction --input data/input_1.txt --output out.txt
Building raw tree ...
Saving data in out.txt ...
Done!
```

Examples of input and output files can be found in data folder.
