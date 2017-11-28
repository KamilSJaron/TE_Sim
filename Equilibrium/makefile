TE_Sim_equil : main.o genome.o population.o chromosome.o locus.o random.o transposon.o mtrand.o FamilyCensus.o
	g++ -g $^ -o $@

%.o: %.cpp
	g++ -g -c $<

clean:
	rm *.o test
