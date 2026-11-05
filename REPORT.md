# Code quality and testing report

## Repository and tested revision
Repository URL:
Final commit identifier:
Final passing Actions run URL:
Deliberate test-failure run URL:
Repair run URL:

## Design changes
The structural change I made was letting the maf and missing rate thresholds be passed in as arguments instead of 
being global variables set at the start of the analysis. In the original, the function f() used y and z without 
them being passed in. Now min_maf and max_missing are arguments to summarize_variant(), so everything is in the 
function call.The second change I made was adding input validation. The original analysis would run on any values 
in the genotypes, even ones that don’t make sense like 3 or -1, and on any threshold values. Now the functions check 
the inputs first and stop with an error if genotypes aren’t 0, 1, 2, or NA, if na_rm isn’t TRUE or FALSE, or if min_maf 
or max_missing is not a number in its allowed range. I also made it so each summary value is only computed once. In the 
original, n_called and missing_rate were calculated in f() and then calculated again in the loop. Now they are only 
calculated once, inside summarize_variant().

## Test reasoning
One test I made is calculate_maf(c(2, 2, 2, 1, 2, 0)) is set to equal 0.25. This one specifically tests the ability of
caculate_maf() to function correctly even if the counted allele is the more common allele and flips p. The counted allele
shows up 9 times out of the 12 alleles present which results in a p of 0.75, but because this is minor allele frequency
maf = 1-p which is 0.25. This would detect the mistake of forgetting to implement that correction and only returning 
counted allele frequency.

A second test I wrote is summarize_variant(c(0, 0, 0, 1), min_maf = 0.125) setting keep to TRUE. This tests the fact that
the threshold for mininmum maf is written inclusively. Because the maf resulting from these genotypes is 0.125 and min_maf
is set to the same value if mistakenly the threshold was set to be uninclusive by use of < instead of <= keep would be 
set to FALSE instead mistakenly.

## Deliberate regression
I decided to change how missingness_rate is calculated because I had originally messed this up when rewriting the calculation.
Instead of calculating correctly by (length(genotypes) - n_called)/length(genotypes), I'm instead calculating the call
rate with (n_called)/length(genotypes). The tests labeled "calculations are correct with full input" and "calculations 
are correct with some missing input" fail with this change. This needs to be an implementation change and not a test
change because this is incorrect calulation of the value this anaylis is designed to report. I restored the calculation
to (length(genotypes) - n_called)/length(genotypes) and everything passed again.

