# LandScape
A simple method to aggregate p values without a priori grouping

For a description of the method read LandScape.pdf

For help on how to run the program run the command: 
python landscape.py --help

If plink files are used as input, the package https://github.com/ymer/Plink is required.

Output files:

1. An output files containing the maximal segments and their values. This file is saved to [output]_[score_evaluation].txt
It contains the following columns:
s - The start point of the maximal segment.
e - The end point of the maximal segment.
Y - The height of the highest point in the maximal segment.
p_value - The p value of the maximal segment.
adjusted_p_value - The p value of the maximal segment adjusted for multiple testing.

2. A file containing the A_k values. This file is saved to [output]_Ak.txt
