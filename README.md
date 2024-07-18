# LSAssigner4iEDA
## Compilation and Execution

### Compilation

To compile the project, run the following commands:

```bash
# Compile lsassigner library
make
# Compile test.cpp to produce test.o
g++ -fPIC -g -O3 -std=c++17 -fopenmp -c test.cpp -o test.o
# Link main.cpp with test.o and lsassigner library to produce executable TA
g++ -fopenmp -m64 -g -O3 -std=c++17 -o TA main.cpp test.o -L. -llsassigner

```

### Execution

To execute the program, run the following command:

```bash
# Run TA with input file ispd18_test5_metal5.input.def.ta.txt and redirect output to out.txt
./TA ispd18_test5_metal5.input.def.ta.txt

```
Replace ispd18_test5_metal5.input.def.ta.txt with your actual input file name.

The output file will be named output.txt and will contain the solution of the problem.