# ising

Multiple examples for Ising model and Metropolis algorithm

* onespin: model for a single spin flipping, see either onespin.cpp or Ising.ipynb
* two spin: two spins flipping, either twospin.cpp or Ising.ipynb
* ising1d: simulation of 1d lattice of spins, see ising1d.cpp or Ising.ipynb or ising1d.ipynb
* ising2d: simulation of 2d lattice of spins, see ising2d.cpp or ising2d.py
* ising2d_vs_T: Metropolis simulation of cooling a 1d lattice, see ising2d_vs_T.cpp or ising2d_vs_T.py
* ising1d_vs_T: Metropolis simulation of cooling a 1d lattice, see ising1d_vs_T.cpp

Complete the exercise for make plots of M,E,C in the 2d lattice using simulated annealing. <br>

This shouldn't be noticeable, but I wasn't able to connect to Rivanna, so I did this assignment in my local build of ROOT then uploaded my script and pdf to my github. I don't think that will affect anything on your end for grading, but if something goes wrong, that's what happened. 

ising.cpp - I like ROOT executable scripts, so instead of modifying ising2d_vs_T.cpp, this script does the plots for ising.pdf. It runs with .x ising.cpp in ROOT
