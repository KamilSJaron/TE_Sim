test: main.o genome.o population.o chromosome.o locus.o random.o transposon.o mtrand.o FamilyCensus.o
	g++ -g  main.o genome.o population.o chromosome.o locus.o random.o transposon.o mtrand.o FamilyCensus.o -o test

main.o: main.cpp
	g++ -g -c main.cpp

genome.o: genome.cpp
	g++ -g -c genome.cpp

population.o: population.cpp
	g++ -g -c population.cpp

chromosome.o: chromosome.cpp
	g++ -g -c chromosome.cpp

locus.o: locus.cpp
	g++ -g -c locus.cpp

random.o: random.cpp
	g++ -g -c random.cpp

transposon.o: transposon.cpp
	g++ -g -c transposon.cpp

mtrand.o: mtrand.cpp
	g++ -g -c mtrand.cpp
	
FamilyCensus.o: FamilyCensus.cpp
	g++ -g -c FamilyCensus.cpp


clean:
	rm *.o test
