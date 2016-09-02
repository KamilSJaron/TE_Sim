# TE_Sim
TE Population Simulator_Written by Elie Dolgin, Modified by Sam Groth

There are two slightly different versions of this program. 

Equilibrium allows you to drop a certain number of TE families into the genome and let the genome equilibrate.

Invasion problems will allow you to run many simulations, where a single TE family is introduce in one copy to the genome. The genome may or may not already have families in there.

To compile:

Place folders (Equilibrium or Invasion_Probs) in the appropriate location. Within the folders, type the command make.

To run:

Within the folder, provide a command file titled: input.txt.

type the command within the folder (if you don’t have your path variable pointing to the folder)

./test

 

Input.txt details

Input.txt

Line 1: N: placeholder, i.e No. MUST BE ‘N’.
Line 2: N (actual population size, not read by genome, read by main)
Line 3: ut ( transposition rate u, read by genome)
Line 4: vt (excision rate v, read by genome)
Line 5: sa (a, selection parameter, read by genome)
Line 6: sb (b, synergism within family selection parameter, read by genome)
Line 7: faf (frequency affecting fitness (1.0 is default, allowing for some to be neutral for values <1.0) ,read by genome)
Line 8: initial family count (initialTE, read by genome)
Line 9: (old run parameter, leave at 0) 
Line 10: (generations)

Not in input (hardwired) in genome.cpp
Numberofchromosomes
Ploidy
Chromlength
Rec Rate:


To change the number per previous family to equilibrium values (if different):
Population.cpp: Poisson(29) value can be changed to Poisson(x)


Regarding equilibrium vs. invasion
In population.cpp 
n = initialTE (from input.txt via Genome.cpp)

if n=0 a for loop is cutout
if n=1 or more, within whichFamily for loop# of TEs initiated for each family within the genome is set at rand.Poisson(29) or rand.Poisson(x) where x is family initial abundance

In genome is: recombination rate, chromosome number, chromosome length, ploidy













