# volterraLMS
Simplee Volterra LMS Filter implementations.

Volterra series are widely used in nonlinear system modeling. The Volterra series expansion of a nonlinear system consists of a nonrecursive series in which the output signal is related to the input signal as

![volterraSeries](Images/volterraSeries.png)

d'(k) is the desired system output, and w<sub>oi</sub> are the **Volterra Kernels** of the system. In order to apply the LMS algorithm to a nonlinear LMS filter, we have to interpret the input vector **x** and the weight vector **w** in a different manner. For the N-th order filter with 2nd order Volterra series, we have:

![volterraLMSVectors](Images/volterraLMSVectors.png)

## References
Diniz, Paulo S. R. - Adaptive filtering: algorithms and practical implementation (2013)
