// *********************************************************************
// 
// population.h
// 
// Created by: Elie Dolgin, University of Edinburgh
//
// First started: March 17, 2005
// Last edited:
//
// *********************************************************************

#ifndef POPULATION_H_EDOLGIN_TE
#define POPULATION_H_EDOLGIN_TE

#include "genome.h"
#include <vector>
#include <stdlib.h>
#include "chromosome.h"
#include "locus.h"
#include "transposon.h"
#include "random.h"
#include "FamilyCensus.h"


class Population {

public:
  Population(int);
  ~Population();
  
  Genome & GetIndividual(int);
  const Genome & GetIndividual(int) const;
  int GetPopSize() const;
  unsigned int GetPopulationTECount() const;
  unsigned int GetPopulationTECountAffectingFitness() const;
  bool isExtinct(int fam);//7/3
  double GetPopulationMeanFitness() const;
  double GetPopulationMeanFamilyCount();//7/3
  
  void Initialize(bool, bool);
  Genome MakeIndividual();
  void DeleteIndividual(int);
  Population * SexualReproduction();
  Population * AsexualReproduction();
  Population * Bottleneck(int, bool, bool);
  void TranspositionAndLoss();
  
  void ListPopulationSites() const;
  void PrintParameters(char[]);
  void SummaryStatistics(char[], int);
  void RecordPopulation(char[], int);
  
  void SummaryStatistics(int);
  void SummaryStatistics(int, int);
  void FamilyStats(int lastFam);
  
  void IntroduceElement(int);
private:
  std::vector<Genome> genoVector;
  int popSize;
  static Random rand;

};


#endif


