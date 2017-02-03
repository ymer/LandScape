LandScape - A simple method to aggregate p values without a priori grouping

Paper: https://www.ncbi.nlm.nih.gov/pubmed/27269897

For help on how to run the program run the command: python landscape.py --help

If plink files are used as input, the package https://github.com/ymer/Plink is required.

Output files:

An output files containing the maximal segments and their values. This file is saved to [output]_[score_evaluation].txt.
It contains the following columns:
s - The start point of the maximal segment.
e - The end point of the maximal segment.
Y - The height of the highest point in the maximal segment.
p_value - The p value of the maximal segment.
adjusted_p_value - The p value of the maximal segment adjusted for multiple testing.

A file containing the A_k values. This file is saved to [output]_Ak.txt

Arguments:

  --infile INFILE       
                        The input file path
                        
  --outfile OUTFILE     
                        The output file path
                        
  --input_format {text,gzip,plink}
                        Default: text
                        
  --transformation {dichotomous,log,none}
                        Transformation of the input data. dichotomous applies the 1/-1 transformation shown example 1, and log applies the log transformation shown in example 2.
                        Default: log
                        
  --score_evaluation {A0,A1,A2}
                        How to evaluate the scores. (See section 4.2-4.4).
                        
  --gamma GAMMA         
                        The threshold of significance. (See Example 1 and 2).
                        Default: 0.05
                        
  --no-lattice          Set if Z is not a lattice variable. (See definition 8 and 9)
                        
  --permutations PERMUTATIONS
                        The number of permutations that are performed if score_evaluation is A1 or A2.
                        Default: If input_format is text or gzip, the number of lines in the input file.
                        
  --random_seed RANDOM_SEED
                        Random seed for permutations if input_format is plink.
                        Default: 1
                        
  --no-draw             Do not draw the landscape
