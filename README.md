# Bayesian-Filtering

This repository provides the code necessary to simulate a variety of different single object tracking filters. 
Bayesian filtering is a method of inferring the current state from a set of observations and the entire history of the 
system. In the repository, I include implementations of the Kalman Filter, Probabilistic Data Association (PDA) Filter, 
and the Bernoulli Filter for both linear and non-linear measurement models. In all of the filters, I am assuming that the 
object to track can vary both its velocity and its position.

## Getting Started

I suggest learning the filters in the following order:
1) Kalman Filter
2) PDA Filter
3) Bernoulli Filter

Always implement the linear case before switching to the non-linear model. This will greatly aid in your ability
to debug your code! The ReadMe within each filter's folder goes into greater detail about the function of each 
MATLAB program.

### Resources

Before attempting to read the following list of resources, it is strongly recommended that the reader has a comfortable
background in probability theory, linear algebra, and calculus. Here are resources that are helpful in understanding each 
filter in greater detail:
1) Kalman Filter
* [MIT Tutorial](http://web.mit.edu/kirtley/kirtley/binlustuff/literature/control/Kalman%20filter.pdf)
* [Detailed Walkthrough](https://www.kalmanfilter.net/default.aspx)
* [YouTube Video](https://youtu.be/FkCT_LV9Syk)

2) Unscented Kalman Filter
* [UKF Paper](https://www.seas.harvard.edu/courses/cs281/papers/unscented.pdf)

3) PDA Filter
* [PDA Filter Paper](http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.212.383&rep=rep1&type=pdf)

4) Bernoulli Filter
* [Bernoulli Filter Paper](http://ba-ngu.vo-au.com/vo/RVVF_Bernoulli_TSP13.pdf)

The concepts from the unscented Kalman Filter also apply to using a nonlinear measurement model in both the PDA Filter and Bernoulli Filter.


