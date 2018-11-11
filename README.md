## : The Gaussian graphical model (GGM)
is one of the common dynamic modelling approaches in the construction of gene networks. In inference of this modelling the interaction between genes can
be detected mainly via graphical lasso (glasso) or coordinate descent-based approaches. Although these
methods are successful in moderate networks, their performances in accuracy decrease when the system
becomes sparser. We here implement a particular type of polynomial transformations, called the Bernstein
polynomials, of the network data in advance of their inference to raise the accuracy. From comparative
Monte Carlo studies and real data analyses we show that these polynomials are successful in terms of the
precision, specificity and F-measure when the scale-free networks are modelled via GGM and estimated
by glasso, and accordingly they can be used as a preprocessing step in inference of these networks. 

In this code it is suggested to perform the Bernstein and Szasz–Mirakyan polynomials in advance
of the inference via GGM to improve the accuracy of the estimates, in particular, when the
complexity of the system raises.

