# Central Limit Theorem Demonstration

This ROOT macro serves as a demonstration of the Central Limit Theorem.
In probability theory, the Central Limit Theorem states that, if you have a random variable X which follows some arbitrary probability distribution P(X), if you get experimental values of X and then average them, and then repeat this process, the resulting distribution of average X's will converge to the normal distribution. The mean of this normal distribution will converge to the mean of our original Probability Density Function (PDF), while the standard deviation will converge to the standard deviation of our PDF divided by the square root of the population size. See [Wikipedia](https://en.wikipedia.org/wiki/Central_limit_theorem) for more information on the mathematical theory.

This program demonstrates the power of this theorem by allowing the user to define their own custom probability density function and run many experiments to recover a normal distribution. It serves as an experimental validation of the central limit theorem. The PDF is defined in `Double_t pdf(Double_t x)` in the file central_lim.C. If you wish to change the probability density function, you must edit the code there.

This program is a macro for the CERN ROOT library. You must have ROOT installed. To run, type
`root central_lim.C`
into the command line.

For an example output plot, see plot.png.
