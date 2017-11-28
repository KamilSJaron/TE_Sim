TE_Sim_invasion : main.o genome.o population.o chromosome.o locus.o random.o transposon.o mtrand.o FamilyCensus.o
	g++ -g $^ -o $@

%.o: %.cpp
	g++ -g -c $<

clean:
	rm *.o test
