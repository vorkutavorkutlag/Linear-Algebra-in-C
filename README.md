# Linear Algebra Demo in C

The naive approach, without any fancy algorithms, using Gaussian elimination as we were taught in the Linear Algebra 1 course.
Originally made to practice programming and linear algebra at the same time due to upcoming exams.

Currently supports:
- Finding the determinant of any given real square matrix.
- Solving any n-th degree real-coefficient polynomial

The original plan was to find all the eigenvalues of a matrix as well, that is left as an exercise to the reader ;]

Example usage:
```
./lac
Square matrix determinant, or solve characteristic polynomial? [M/P]
M
Enter dimension of real square matrix:
3
Enter Matrix (same as visually):
1 2 3
3 2 1
4 4 4

...

REF FORM
1 2 3
0 -4 -8
0 0 0
Determinant: 0

./lac
Square matrix determinant, or solve characteristic polynomial? [M/P]
P
Enter polynomial degree:
2
Enter real polynomial coefficients in rising degree order:
1 0 1

Root 0: -0.000000 + 1.000000i
Root 1: 0.000000 + -1.000000i
...
```
