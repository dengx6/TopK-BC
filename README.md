# TopK-BC

## Compile & Execute

Our framework requires c++17 and GCC 9.x (or later). Under the directory /src/topkbc, execute the following commands to compile the source code and execute TopKBC.

```shell
make
../bin/topkbc
```


### Basic parameter description

Please refer to our paper for a detailed explanation of the basic parameters

## Input File Format
The input graphs have edge weights and timestamps. Each edge on the data graph has a weight and timestamp. Each vertex is represented by a distinct unsigned integer (from 0 to 4294967295). There is at most one edge between two arbitrary vertices. 

For example, the file format is as follows:

```
0 0 3 1
0 1 2 2
1 0 3 3
1 1 4 4
```
Each line represents an edge, i,e, `<vertexID, vertexID, weight, timestamp>`;



