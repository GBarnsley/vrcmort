
We have added the capacity to apply a monotonic structure to the coefficients for given parameters.
Now we need to make it more specific.

We can assume that the order in which the monotonic coefficients are given (i.e in `mortality_monotonic`) is the decreasing order they follow.
These will represent level of functionality of healthcare facilities, the largest value (first in `mortality_monotonic`) represents the number (standardised) of perfectly functioning facilities. The final level (which will not be provided) and will be treated as the baseline value (i.e. when all monotonic variables are 0.0) is the number of non-functioning facilities.
The calculated Beta_mort_(i-1) should be < Beta_mort_i since better facilities should decrease mortality by more than less functional facilities.
The calculated Beta_rep_(i-1) should be > Beta_rep_i since better facilities should increase the reporting rate by more than less functional facilities.

These requirements should be documented, variable names changed to reflect facilities and not generic monotonic effects.
Ensure that these relationships hold in the stan model

If you have any questions or confusions, please check with me first.

# Steps:
1. Determine the correct approach to including these variables
2. update stan model
3. update R code (use the r coding skill)
4. update tests for this new feature
5. ensure tests pass
5. check your work with @clean-code-reviewer
6. format any changed files with `air`
7. Add a vingette explaining this feature
