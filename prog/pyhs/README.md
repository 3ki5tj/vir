# PY equation of states of hard-sphere fluids #


File       |  Description
-----------|--------------------------------------------------------
pyhs.ma    |  the original version. Intepolation along the packing fraction t. The polynomial obtained at a particular t is normalized, such that when X = 0 it matches the value of polynomial for the radius of convergence.
pyhs2.ma   |  the improved version. Intepolation along the packing fraction t. The polynomial obtained at a particular t is normalized against the leading term of X.
getvir.ma  |  handle the virial coefficients
