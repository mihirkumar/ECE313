#1.1a
>What differences do you see between the pmf plotted in part 1 above and the estimated pdf?

The values on the pmf vary greatly above and below with respect to the smoother estimated pdf.

> What differences do you see as the sample sizes increase from 70 to 30000?

With more samples, the curve becomes smoother and begins to looks more distinctly Gaussian.

#1.1b
> 1. Explain why tabulate(X) gives a pdf while the other gives a pmf.

Tabulate(X) makes the time into discrete values, but without the floor function, time is continuous and therefore easier displayed in a pdf.

> 2. What difference do you find for the derived min, max for each distribution?

The min and max for distribution that utilized the floor function was much wider. This appears to be a result of the floor function grouping similar values together instead of keeping each value explicitely seperae. The tabulate(x) distribution has 30,000 unique values, whereas the tabulate(floor(x0) distribution has only 76, creating the disparity in min and max shown.

> 3. Can you relate the difference to identify a distinct characteristics of a pdf?

The pdf always tends to lag behind the pmf.


#1.2c
> What differences do you see between the distribution generated in this task and the plot from Task 1.1 for 30K sample size?

The distribution generated in this task is a lot more uniform than the 30k sample size plot form task 1.1, which exhibited a left-ward skew.


#1.2d
> What differences do you see compared to the values of a and b in Task 1.1?

The value for A is slightly higher, and the value for B is lower in this task. 
This could be due to random generation, especially since samples are more sparse in the very outer edges of a Gaussian distribution.
