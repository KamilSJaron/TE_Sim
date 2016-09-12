#include "population.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <stdlib.h> 

#include "random.h"
#include "time.h"

using namespace std;

int main ()
{
    
  clock_t time;
  double seconds;
  int trials=10000;;
  int successes=0;
  int failures=0;
  int genstatsF=0;
  int genstatsS=0;
    
  cout << endl;
  time = clock();
  for (int run=1; run<=trials; run++)
  {
      

    std::cout<<"Trial #: "<<run<<std::endl;
	std::ifstream fin("input.txt");
	if (! fin.is_open())
	  { cout << "Error opening file"; exit (1); }
	
	// Initialize population size & whether to generate new population or load from file
	int N,generations;
    double invRate;
    int initialFams;

	bool fromFile = true;
	char tempChar[100];
	fin.getline(tempChar,100);
	if (tempChar[0]=='Y') fromFile = true;
	else if (tempChar[0]=='N') fromFile = false;
	fin.getline(tempChar,100);
	N= strtol(tempChar,0,10);
   
    fin.getline(tempChar,100);//ut
    fin.getline(tempChar,100);//vt
    fin.getline(tempChar,100);//sa
    fin.getline(tempChar,100);//sb
    fin.getline(tempChar,100);//faf
    fin.getline(tempChar,100);//initialTE
    initialFams = strtol(tempChar,0,10);
    fin.getline(tempChar,100);//invRate
    invRate = strtod(tempChar,0);
    fin.getline(tempChar,100);//generations
    generations = strtol(tempChar,0,10);
      
	fin.close();
	Population * pop = new Population(N);
	Population * temp;
	
	int size = pop->GetPopSize();
	bool clonal = false;//Genome::clonal;
	
	pop->Initialize(clonal, fromFile);
    

//	pop->PrintParameters("detailed.txt");
      if (run==1) {
          pop->PrintParameters("summary.txt");
          std::cout<<"\n Initial families = "<<initialFams<<std::endl;
      }
//	pop->SummaryStatistics("detailed.txt", 0);

    int invasionsPerGeneration = N*invRate;
	int gen = 1;
    int tefamily=0;
    
	for (gen; gen <= generations; gen++)
	{

/*
		if (pop->GetPopulationTECount() == 0)
		{
		  cout << "No TEs at generation [" << gen << "]." << endl << endl;
//		  pop->SummaryStatistics("summary.txt", gen);
//		  pop->SummaryStatistics("detailed.txt", gen);
          failures++;
          genstatsF+=gen;
          std::cout<<"Running totals- Failures: "<<failures<<", Successess: "<<successes<<std::endl;
		  break;
		}
    
        if((pop->GetPopulationTECount()/size) > 1){
            cout<<" Success at generation ["<<gen<<"].\n\n";
//            pop->SummaryStatistics("summary.txt", gen);
//            pop->SummaryStatistics("detailed.txt", gen);
            successes++;
            genstatsS+=gen;
            std::cout<<"Running totals- Failures: "<<failures<<", Successess: "<<successes<<std::endl;
            break;
        }
 
*/
        
 
        if(pop->isExtinct(0)){
            cout << "No TEs of invading family at generation [" << gen << "]." << endl << endl;
            failures++;
            genstatsF+=gen;
            std::cout<<"Running totals- Failures: "<<failures<<", Successess: "<<successes<<std::endl;
            break;

        }
        
        
        
        if(pop->GetPopulationMeanFamilyCount()>=(initialFams+1)){
            cout<<" Success at generation ["<<gen<<"].\n";
            successes++;
            genstatsS+=gen;
            std::cout<<"Mean family cout = "<<pop->GetPopulationMeanFamilyCount()<<std::endl;
            std::cout<<"Running totals- Failures: "<<failures<<", Successess: "<<successes<<std::endl;
            break;
        }
 
 
 
		// REPRODUCTION
        
        
        
		temp = pop->SexualReproduction();
		delete pop;
		pop = temp;
         

        
		// TRANSPOSITION & LOSS
		pop->TranspositionAndLoss();
		cout << ".";
        
        
        
/*
		if (gen%100==0)
		{
            
		  cout << endl;
		  pop->SummaryStatistics("detailed.txt", gen);

            
		}
 

		if (gen==generations)
		{
		  pop->RecordPopulation("2200WITHSELECTION.txt", gen);
          pop->FamilyStats(tefamily);
		}
*/
      if(gen%100==0)
          cout<<endl;
        
	}
	
	delete pop;
    
    
  }
    time = clock()-time;
    std::cout<<"GRAND TOTALS AND ACCOUNTING\n";
    std::cout<<"Failures: "<<failures<<std::endl<<"Successes: "<<successes<<std::endl;
    
    if(failures!=0){
        genstatsF/=failures;
        std::cout<<"Average generation at failure = "<<(genstatsF)<<std::endl;
    }else{
        std::cout<<"No failures!\n";
    }
    
    if(successes!=0){
        genstatsS/=successes;
        std::cout<<"Average generation where success was defined = "<<(genstatsS)<<std::endl;
    }else{
        std::cout<<"No successes!\n";
    }
    printf ("\nTotal runtime was: %d clicks (%f seconds).\n",time,((float)time)/CLOCKS_PER_SEC);
    
    
  return 0;
}
