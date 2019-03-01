# ParGNG
Parallel Growing Neural Gas

ParGNG is developed at Czech national supercomputing centre IT4Innovations. The goal of this application is to contribute to the area of Artificial Neural networks by speed up computing by using High Performance Computing and improvements for sparse data with higher dimensions. 

The code has been developed and optimized for the [IT4Innovations Salomon cluster](https://docs.it4i.cz/salomon/hardware-overview/).

# Dependence
1) GCCcore/6.3.0   
2) binutils/2.27-GCCcore-6.3.0   
3) GCC/6.3.0-2.27   
4) numactl/2.0.11-GCCcore-6.3.0   
5) hwloc/1.11.7-GCC-6.3.0-2.27   
6) OpenMPI/1.10.7-GCC-6.3.0-2.27


## Build instructions
1. Clone the repository
2. Run make

After successfull compilation the `run` directory should contain one executable 'gng_linux'.

## Usage

mpirun -n 20 gng_linux -c < file > -o < file > -i < input >

**List of basic parameters:**
```
-c <name> - name of configuration file
-o <name> - name of output file
-i <name> - name of input file
-debug    - enable debug reports (Default false)
-showInput -convert input to GDF file
-withSom  - combine GNG and SOM
-s        - Change size for GDF
```

**Input configuration file example**
```
Number or records
Number of dimension
Name of data file
Output newtork configuration file
Output network data file
Number of parts
```
```
800
2
TwoDiamonds-sparseFormat.txt
diamanty1.all
diamanty1-1.parts
9
```

**Input data file example**

Each record is on a separate line. 

`<position>:<value> <position>:<value> ...`
```
0:0.222222222 1:0.625 2:0.06779661 3:0.041666667
0:0.166666667 1:0.416666667 2:0.06779661 3:0.041666667
0:0.111111111 1:0.5 2:0.050847458 3:0.041666667
0:0.083333333 1:0.458333333 2:0.084745763 3:0.041666667
```

**Configuration file example**
Each record is on a separate line.

```
gama
e_w 
e_n 
alpha
beta
a_max
NeuronMaxCount
MaxIteration
```
```
200
0.05
0.006
0.5
0.0005
300
500
200
```
**Outputs**
Default format for output is Gephi files- gdf and gn

```
nodedef> name VARCHAR,label VARCHAR,width DOUBLE,height DOUBLE,x DOUBLE,y DOUBLE,color VARCHAR
0,"n0",10.0,10.0,1971.48,1972.38,'153.153.153'
1,"n1",10.0,10.0,1213.11,4955.56,'153.153.153'
2,"n2",10.0,10.0,3776.39,4933.86,'153.153.153'
3,"n3",10.0,10.0,2002.02,2499.59,'153.153.153'
4,"n4",10.0,10.0,3479.01,2986.04,'153.153.153'
5,"n5",10.0,10.0,3780.64,39.1175,'153.153.153'
6,"n6",10.0,10.0,1219.22,48.9155,'153.153.153'
7,"n7",10.0,10.0,992.712,3007.58,'153.153.153'
```


# Licence

See the LICENSE file at the root directory.
