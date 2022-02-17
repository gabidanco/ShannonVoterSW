This program is written in fortran 90 and evaluates the Shannon entropy 
for the stationary state of the voter dynamics with external nodes in 
small world networks.

Some more detailed analysis can be found in the paper by Franco, 
Marquitti, Fernandes, Braha, and de Aguiar, 2021, "Shannon information 
criterion for low-high diversity transition in Moran and voter models".
DOI: https://doi.org/10.1103/PhysRevE.104.024315

Any suggestions, comments or questions can be sent to:
gabidanco@gmail.com.

________________________________________________________________________

To run the program:

1) Compile the file 'shannon_voter_sw2kd.f90' in your operating system.
2) Put the executable file and the input.in file in the same directory.
3) Run the executable file 'shannon_voter_sw2kd.exe'

--------Input files-----------

---> "seed.in"
- this is a file with 12 integer numbers separated by blank spaces.
- this file is updated at the end of each run.

---> "input.in"
- n: total number of dynamic nodes
- kn: degree of the network will be 2*kn
- nt: equilibrium time 
- nb: number of measures

---------Output file----------

---> "SW_Shannon_***_###.dat"
- ***: represents the rewiring probability for this run of the simulation
- ###: indicates the sample for this run - program runs 10 different networks
generated from the same rewiring probability
- First column represents external influence (symmetric)
- Second column contains Shannon entropy calculated for such influence

