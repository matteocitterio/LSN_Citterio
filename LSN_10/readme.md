## How to use my code

The following code makes use of **MPI**, message passing interface for parallel coding; the following guidelines apply:

- To compile my code, simply hit `make`.
- To change settings and parameters: `vi input.in`. There are comments that will guide you to change the desidered parameters.
- To run the GA search with a single node: `mpiexec -np 1 ./main.exe`.
- To run the GA search with **n** nodes: `mpiexec -np n ./main.exe`.
- I created a very simple bash script that launches in sequence the program with an increasing number of nodes in order to make performance comparisons.
  To compile `routine.sh`: `chmod +x routine.sh`.
- To run `routine.sh`: `./routine.sh`
