
You are an expert computation bayesian statistician with years of stan experience.

Your task is to incorporate an optional monotonic effect into two levels of the stan model defined in @inst/stan.

There are two regressions one that makes up the overall death rate (lambda) and one that informs the reporting rate (rho)

We want to be able to include co-variates with know monotonic levelling in either or both of these regressions.
The monotonic variables be real positive values that have K levels. So each entry with have K real value variables.
We want the monotonic structure to follow that beta_1 > beta 2 > ... > beta_(k-1) where the kth value is the baseline.
We should infer beta_i using B * sum_to_i(b) where b is a simplex.
See @.gemini/examples/Estimating Monotonic Effects with brms • brms.html for inspiration

# Steps:
1. Determine the correct approach to including these variables
2. update stan model
3. update R code (use the r coding skill)
4. write tests for this new feature
5. check your work with @clean-code-reviewer
6. format any changed files with `air`
