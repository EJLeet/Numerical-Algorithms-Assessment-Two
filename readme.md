Numerical Algorithms Assessment 2

AUTHOR

    Ethan Leet
    s5223103


ABOUT

    This assignment involved creating five programs to solve and experiment on five numerical based algorithms. As well as creating c++ source code to solve these problems, extensive experimentation and theoretical analysis was to be conducted, and is included in the report attached to this submission. Python code was used to evaluate the output of the c++ code and to create a visiual representation of the data. Question 6 is contained in a seperate repository. 

Compiling

    Each program was compiled using clang++ and the c++ standard library 14 without errors nor warnings. To compile a program, you first must be within the program directory and can then execute the following command, replacing x with the question you wish to compile:

    clang++ main.cpp -std=c++14 -o questionx -Ofast

    Alternitavely, if CMake is installed, you can run the 'make' command from either the individual question folder or the parent folder. Running the make command in the question folder will execute the above command whereas executing the make command from the parent folder will compile all six problems in this submission and place the binaries and the respective question folder.

Run Time

    After compilation, the below command can be executed from the question folder, replacing x with the question you wish to run:

    ./questionx
