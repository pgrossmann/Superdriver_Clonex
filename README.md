#Foobar

Clonex is a C routine to simulate tumorigenesis with superdrivers according to clonal evolution using a Wright-Fisher model.

## Installation

Download root and execute 

```bash
make
```
to use the Makefile to produce a clonex.o from te clonex.c file.

##Usage

For a list of options, execute

```bash
./clonex -h
```

An example simulation can be run with

```bash
./clonex -s 0.01 -c 2 -u 0.00000001 -g 4500 -R 50 -N 1000000000 -n 1000000 -M 1000 -f tmp
```
for a simulation with 4500 generations, from a population of 10e6 cells to 10e9 cells using a base selection of 0.01 and a superdriver factor of 2 with 50 replicates each (and without simulating the mutator phenotype) that will be stored ina 'tmp' folder.

## Wrapper scripts

There are two wrapper scripts available. The first one is 

```bash
perform_simulations.sh
```

and wraps clonex. The second one is

```bash
perform_several_simulations.sh
```

and wraps perform\_simulations.sh
