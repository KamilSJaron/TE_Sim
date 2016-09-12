// *********************************************************************
// 
// genome.cpp
// 
// Created by: Elie Dolgin, University of Edinburgh
//
// First started: March 11, 2005
// Last edited:
//
// *********************************************************************

#include "genome.h"
#include <math.h>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <algorithm>
#include "chromosome.h"
#include "locus.h"
#include "transposon.h"
#include "random.h"


//int Genome::N = 0;
double Genome::ut = 0;
double Genome::vt = 0;
double Genome::sa = 0;
double Genome::sb = 0;
double Genome::faf = 0;
int Genome::initialTE = 0;

//double Genome::ut = 0.01;
//double Genome::vt = 0.001;
//double Genome::sa = 0.001;
//double Genome::sb = 0.00018;
//double Genome::faf = 1.0;
//int Genome::initialTE = 50;
//int Genome::N = 0;

int Genome::numberOfChromosomes = 3;
int Genome::ploidy = 2;
int Genome::chromLengths[3] = {10000,10000,10000};
double Genome::chromRecRates[3] = {0.0001,0.0001,0.0001};
//double Genome::rGenome = 0.01;

bool Genome::clonal = false;
bool Genome::parametersSet = false;

Random Genome::rand;

void Genome::SetParameters()
{
    std::ifstream fin("input.txt");
	if (! fin.is_open())
	  {std::cout << "Error opening file"; exit (1); }
    
	char tempChar[100];
    
	while(!fin.getline(tempChar, 100).eof())
	{
    //fin.getline(tempChar,100);

    //N=strtol(tempChar,0,10);
	fin.getline(tempChar,100);
	fin.getline(tempChar,100);
	ut=strtod(tempChar,0);
	fin.getline(tempChar,100);
	vt=strtod(tempChar,0);
	fin.getline(tempChar,100);
	sa=strtod(tempChar,0);
	fin.getline(tempChar,100);
	sb=strtod(tempChar,0);
	fin.getline(tempChar,100);
	faf=strtod(tempChar,0);
	fin.getline(tempChar,100);
	initialTE=strtol(tempChar,0,10);
    fin.seekg(0,fin.end);
    
	}
	fin.close();
    
	parametersSet = true;
}

Genome::Genome():
totalTECount(0),
numberOfFamilies(0)
{
    headFamily=0;
	chromoVector.resize((numberOfChromosomes*ploidy));
	if (!parametersSet)
	  SetParameters();
	
	for (int i=1; i <= numberOfChromosomes; i++)
	{
	  if (ploidy == 1)
	  {
		chromoVector.at(i-1).SetChromNumberAndCopy(i,1);
		chromoVector.at(i-1).SetChromLengthAndRecRate(chromLengths[i-1], chromRecRates[i-1]);
	  }
	  if (ploidy == 2)
	  {
	    chromoVector.at(2*i-2).SetChromNumberAndCopy(i,1);
		chromoVector.at(2*i-2).SetChromLengthAndRecRate(chromLengths[i-1], chromRecRates[i-1]);
		chromoVector.at(2*i-1).SetChromNumberAndCopy(i,2);
		chromoVector.at(2*i-1).SetChromLengthAndRecRate(chromLengths[i-1], chromRecRates[i-1]);
	  }
    }
}

Genome::Genome(int num, int pl):
totalTECount(0),
numberOfFamilies(0)
{
    headFamily=0;

    numberOfChromosomes = num;
	ploidy = pl;
	
	chromoVector.resize((numberOfChromosomes*ploidy));
	if (!parametersSet)
	  SetParameters();
	
	for (int i=1; i <= numberOfChromosomes; i++)
	{
	  if (ploidy == 1)	
	  {
		chromoVector.at(i-1).SetChromNumberAndCopy(i,1);
		chromoVector.at(i-1).SetChromLengthAndRecRate(chromLengths[i-1], chromRecRates[i-1]);
	  }
	  if (ploidy == 2)
	  {
	    chromoVector.at(2*i-2).SetChromNumberAndCopy(i,1);
		chromoVector.at(2*i-2).SetChromLengthAndRecRate(chromLengths[i-1], chromRecRates[i-1]);
		chromoVector.at(2*i-1).SetChromNumberAndCopy(i,2);
		chromoVector.at(2*i-1).SetChromLengthAndRecRate(chromLengths[i-1], chromRecRates[i-1]);
	  }
    }
}

Genome::Genome(const Genome & rhs):
totalTECount(0),
numberOfFamilies(0)
{
    headFamily=0;

	chromoVector.resize((numberOfChromosomes*ploidy));
	if (!parametersSet)
	  SetParameters();
	
	Locus * current;
	
	for (int i=1; i <= numberOfChromosomes; i++)
	{
	  if (ploidy == 1)	
	  {
		chromoVector.at(i-1).SetChromNumberAndCopy(i,1);
		chromoVector.at(i-1).SetChromLengthAndRecRate(chromLengths[i-1], chromRecRates[i-1]);
		current = rhs.GetChromosome(i,1).GetHeadLocus();
		while (current != 0)
		{
			//chromoVector.at(i-1).Insert(current->GetTransposon());
            Insert((i-1),current->GetTransposon());
			current = current->GetNext();
		}
	  }
	  if (ploidy == 2)
	  {
	    chromoVector.at(2*i-2).SetChromNumberAndCopy(i,1);
		chromoVector.at(2*i-2).SetChromLengthAndRecRate(chromLengths[i-1], chromRecRates[i-1]);
		current = rhs.GetChromosome(i,1).GetHeadLocus();
		while (current != 0)
		{
			//chromoVector.at(2*i-2).Insert(current->GetTransposon());
            Insert((2*i-2),current->GetTransposon());
			current = current->GetNext();
		}
		chromoVector.at(2*i-1).SetChromNumberAndCopy(i,2);
		chromoVector.at(2*i-1).SetChromLengthAndRecRate(chromLengths[i-1], chromRecRates[i-1]);
		current = rhs.GetChromosome(i,2).GetHeadLocus();
		while (current != 0)
		{
			//chromoVector.at(2*i-1).Insert(current->GetTransposon());
            Insert((2*i-1),current->GetTransposon());
			current = current->GetNext();
		}
	  }
    }
}

Genome::~Genome()
{
    /*
    if(headFamily!=0){
        std::cout<<"Head present.  Family = "<<headFamily->GetFamily()<<std::endl;
        if(headFamily->GetNext()!=0)
            std::cout<<" Next present.  Family = "<<headFamily->GetNext()->GetFamily()<<".\n";
    }
    std::cout<<GetNumberOfFamilies()<<std::endl;
    */
    delete headFamily;
    
    headFamily = 0;
}

unsigned int Genome::GetNumberOfChromosomes()
{
	return numberOfChromosomes;
}

unsigned int Genome::GetPloidy()
{
	return ploidy;
}

double Genome::GetFAF()
{
	return faf;
}

unsigned int Genome::GetGenomeTECount() const
{
    return totalTECount;
    /*
     //original
	unsigned int genomeTEcount = 0;
	for (int i=1; i <= numberOfChromosomes; i++)
	{
	  genomeTEcount += GetChromosome(i, 1).GetChromTECount();
	  if (ploidy==2)
	    genomeTEcount += GetChromosome(i, 2).GetChromTECount();
	}
     //some debugging nonsense- change GenNumberOfFamilies, too.  This ain't right.
    unsigned int numberwang = GetNumberOfFamilies();
    
    bool totalequalsgenome = (genomeTEcount==totalTECount);
    bool genomeequalsfamcount = (genomeTEcount==numberwang);
    bool totalequalsfamcount = (totalTECount==numberwang);
    
    

    if(totalequalsfamcount){
        std::cout<<"Running totals agree, at least.";
    }else{
        std::cout<<"I give up.";
    }
    if(totalequalsgenome){
        std::cout<<"That's good.";
    }else{
        std::cout<<"That's bad.";
    }
    if(genomeequalsfamcount){
        std::cout<<"Weird, but workable.";
    }else{
        std::cout<<"Expected and sad.";
    }
    
        if(!(totalequalsfamcount&&totalequalsgenome&&genomeequalsfamcount)){
            for(int j=0;j<1000;j++)
                std::cout<<"Failure.\n";
            std::cout<<"\n\nmy count = "<<totalTECount<<", my fam count = "<<numberwang<<", his count = "<<genomeTEcount<<std::endl<<std::endl;
        }else{
            
        }
    
    std::cout<<std::endl;
	return genomeTEcount;
    */
     
}

unsigned int Genome::GetGenomeTECountAffectingFitness() const
{
    
    if(GetFAF()==1.0){
        return totalTECount;
    }
    
    
    unsigned int genomeTEcount = 0;
	for (int i=1; i <= numberOfChromosomes; i++)
	{
	  genomeTEcount += GetChromosome(i, 1).GetChromTECountAffectingFitness();
	  if (ploidy==2)
	    genomeTEcount += GetChromosome(i, 2).GetChromTECountAffectingFitness();
	}
	return genomeTEcount;
     
}

// input is chromosome number and copy number
const Chromosome & Genome::GetChromosome(int num, int copy) const
{
	if (ploidy == 1)
	  return chromoVector.at(num+copy-2);
	else // if (ploidy == 2)
	  return chromoVector.at((2*num)+copy-3);
}

Chromosome & Genome::GetChromosome(int num, int copy)
{
	if (ploidy == 1)
	  return chromoVector.at(num+copy-2);
	else // if (ploidy == 2)
	  return chromoVector.at((2*num)+copy-3);
}

double Genome::GetGenomeFitness() const
{
	unsigned int genomeTEcount = GetGenomeTECountAffectingFitness();
	
	//return (1 - 0.001*pow(genomeTEcount, 1.5));
	
	// assume synergistic epistatic selection
	//return exp ( -(sa * genomeTEcount) - (0.5 * sb * pow(genomeTEcount,2) ) );
    
    
    
    if(genomeTEcount==0){
        return 1.0;
    }
    
    double tbr = 1.0;
    FamilyCensus * current = headFamily;
    while(current!=0){
        tbr*= (exp ( -(sa * current->GetCount()) - (0.5 * sb * pow(current->GetCount(),2) ) ));
        current = current->GetNext();
    }
    /*
    double trouble = exp ( -(sa * genomeTEcount) - (0.5 * sb * pow(genomeTEcount,2) ) );
    std::cout<<"The new fitness is ";
    if(tbr>trouble){
        std::cout<<" GREATER THAN ";
    }else if(tbr==trouble){
        std::cout<<" EQUAL TO ";
    }else{
        std::cout<<" LESS THAN ";
    }
    std::cout<<" the old fitness.\n";
    */
    
    return tbr;
 
}

// overloaded function when mean TE count is constant
double Genome::GetGenomeFitness(int meanCount) const
{
	//return (1 - 0.001*pow(meanCount, 1.5));
    // assume synergistic epistatic selection
	return exp ( -(sa * meanCount) - (0.5 * sb * pow(meanCount,2) ) );
}

void Genome::SetChromosome(Chromosome & c)
{
	int num = c.GetChromNumber();
	int copy = c.GetChromCopy();
	
	if (ploidy == 1)
	  chromoVector.at(num+copy-2) = c;
	if (ploidy == 2)
	  chromoVector.at((2*num)+copy-3) = c;
}

//void Genome::SetGenomeParameters(int num, int pl)
//{
//	numberOfChromosomes = num;
//	ploidy = pl;
//	
//	chromoVector.resize((numberOfChromosomes*ploidy));
//	
//	for (int i=1; i <= numberOfChromosomes; i++)
//	{
//	  if (ploidy == 1)
//	  {	
//		chromoVector.at(i-1).SetChromNumberAndCopy(i,1);
//		chromoVector.at(i-1).SetChromLengthAndRecRate(chromLengths[i-1], chromRecRates[i-1]);
//	  }
//	  if (ploidy == 2)
//	  {
//	    chromoVector.at(2*i-2).SetChromNumberAndCopy(i,1);
//		chromoVector.at(2*i-2).SetChromLengthAndRecRate(chromLengths[i-1], chromRecRates[i-1]);
//		chromoVector.at(2*i-1).SetChromNumberAndCopy(i,2);
//		chromoVector.at(2*i-1).SetChromLengthAndRecRate(chromLengths[i-1], chromRecRates[i-1]);
//	  }
//    }
//}

void Genome::Transpose()
{
	unsigned int num=0, copy=0, pos=0, totalLength=0, currentLength=0;
	bool affectW = false;
	unsigned int teCount = GetGenomeTECount();
	unsigned int transposeCount = (int)rand.Poisson(ut*teCount);
    if(teCount==0){//edited
        return;
        
    }
	
	for (int i=1; i <= numberOfChromosomes; i++)
	  totalLength += GetChromosome(i, 1).GetLength();
	
	if (transposeCount > teCount)
	  transposeCount = teCount;  
	if (transposeCount > (2*totalLength - teCount))
	  transposeCount = (2*totalLength - teCount);
	  
	// for number of transpositions, randomly insert into the genome
	for (int j=0; j < transposeCount; j++)
	{
	  do {
	    pos = (int)((rand.Uniform()*totalLength) + 1);
		num = 1;
		for (int k=1; k <= numberOfChromosomes; k++)
		{
			currentLength = GetChromosome(k, 1).GetLength();
			if (pos > currentLength)
			{
				num++;
				pos -= currentLength;
			}
			else
			  break;
		}
		copy = (int)((rand.Uniform())*(ploidy) + 1);
	  } while (!GetChromosome(num, copy).TestEmpty(pos));//TestEmpty needs to be changed for families
		 
		if (faf > rand.Uniform())
		  affectW = true;
		else
		  affectW = false;
		
		//GetChromosome(num, copy).Insert(Transposon(pos, affectW,GetRandomTransposonFamily(teCount)));//getrandomtransposonfamily
        Insert(num,copy,Transposon(pos,affectW,GetRandomTransposonFamily(teCount)));
	}
}

// Overloaded Transpose method used for when at equilibrium copy number
/*
void Genome::Transpose(double rate, double meanCopyNumber)
{
	int num=0, copy=0, pos=0, currentLength=0, totalLength=0;
	bool affectW = false;
	unsigned int transposeCount = (int)rand.Poisson(rate*meanCopyNumber);
	
	for (int i=1; i <= numberOfChromosomes; i++)
	  totalLength += GetChromosome(i, 1).GetLength();
	
	if (transposeCount > meanCopyNumber)
	  transposeCount = (int)meanCopyNumber;  
	if (totalLength < (transposeCount + meanCopyNumber))
	  transposeCount = (totalLength - (int)meanCopyNumber);
	
	// for number of transpositions, randomly insert into the genome
	for (int j=0; j < transposeCount; j++)
	{
	  do {
	    pos = (int)((rand.Uniform()*totalLength) + 1);
		num = 1;
		for (int k=1; k <= numberOfChromosomes; k++)
		{
			currentLength = GetChromosome(k, 1).GetLength();
			if (pos > currentLength)
			{
				num++;
				pos -= currentLength;
			}
			else
			  break;
		}
		copy = (int)((rand.Uniform())*(ploidy) + 1);
	  } while (!GetChromosome(num, copy).TestEmpty(pos));
		
		if (faf > rand.Uniform())
		  affectW = true;
		else
		  affectW = false;
		unsigned int teCount = GetGenomeTECount();
		GetChromosome(num, copy).Insert(Transposon(pos, affectW,GetRandomTransposonFamily(teCount)));//here
	}
}
*/
void Genome::ElementLoss()
{
	if (vt == 0)
	  return;
	
	unsigned int lossCount=0, chromTEcount=0, nthTE=0;
	
	for (int i=1; i <= numberOfChromosomes; i++)
	{
		if (ploidy == 1)
		{
			chromTEcount = chromoVector.at(i-1).GetChromTECount();
			lossCount = (int)rand.Poisson(vt*chromTEcount);
			if (lossCount > chromTEcount)
			  lossCount = chromTEcount;
		
			for (int k=0; k < lossCount; k++)
			{
			  nthTE = (int)(rand.Uniform()*chromTEcount + 1);
			  //chromoVector.at(i-1).Delete(nthTE);
              Delete((i-1),nthTE);
			  chromTEcount--;
			}
		}
		if (ploidy == 2)
		{
		  for (int j=1; j <= ploidy; j++)
		  {
			chromTEcount = chromoVector.at(2*i+j-3).GetChromTECount();
			lossCount = (int)rand.Poisson(vt*chromTEcount);
			if (lossCount > chromTEcount)
			  lossCount = chromTEcount;

			for (int k=0; k < lossCount; k++)
			{
			  nthTE = (int)(rand.Uniform()*chromTEcount + 1);
			  //chromoVector.at(2*i+j-3).Delete(nthTE);
              Delete((2*i+j-3),nthTE);
			  chromTEcount--;
			}
		}
	  }
    }
}

void Genome::Recombination()
{
	if (ploidy != 2)
	  return;
	
	Locus * current1 = 0;
	Locus * current2 = 0;
	Locus * next1 = 0;
	Locus * next2 = 0;
	double site = 0.0, r = 0.0;
	int events = 0, length = 0;
	
	std::vector<int> sitesVector;
	
	  for (int i=1; i <= numberOfChromosomes; i++)
	  {
	    // if no TE's, then recombination does nothing.
		if ((GetChromosome(i, 1).GetChromTECount() + GetChromosome(i, 2).GetChromTECount()) == 0)
		  return;
		
		// number of recombination events taken from a poisson distribution
		r = GetChromosome(i, 1).GetRecRate();
		length = (GetChromosome(i, 1).GetLength() - 1);
		events = (int)rand.Poisson(r*length);
		
		sitesVector.resize(events);
		
		for (int j=0; j < events; j++)
		{
			int pos = (int)((rand.Uniform())*(length) + 1);
			sitesVector.at(j) = pos;
		}
		
		std::sort(sitesVector.begin(), sitesVector.end());		
		
		current1 = GetChromosome(i, 1).GetHeadLocus();
		current2 = GetChromosome(i, 2).GetHeadLocus();

		for(int a=0; a < events; a++)
		{
			site = sitesVector.at(a) + 0.5;
			//std::cout << site << std::endl;
			
			if ((current1 == 0) && (current2 == 0))
			  break;
			
			if (current1 != 0)
			  if (current1->GetTransposon().GetLocation() < site)
			  {
			    next1 = current1->GetNext();
				while(next1 != 0)
				{
				  if (next1->GetTransposon().GetLocation() > site)
					break;
				  current1 = next1;
				  next1 = current1->GetNext();
				}
			  }
			
			if (current2 != 0)
			  if (current2->GetTransposon().GetLocation() < site)
			  {
			    next2 = current2->GetNext();
				while(next2 != 0)
				{
				  if (next2->GetTransposon().GetLocation() > site)
					break;
				  current2 = next2;
				  next2 = current2->GetNext();
				}
			  }

			if ((current1 == 0) || (current2 == 0))
			{
			  if ((current1 == 0) && (current2->GetTransposon().GetLocation() > site))
				continue;
			
			  if ((current1 == 0) && (current2->GetTransposon().GetLocation() < site))
			  {
				current2->SetNext(current1);
				current1 = next2;
				GetChromosome(i, 1).SetHeadLocus(next2);
				continue;
			  }
			
			  if ((current2 == 0) && (current1->GetTransposon().GetLocation() > site))
				continue;
			
			  if ((current2 == 0) && (current1->GetTransposon().GetLocation() < site))
			  {
				current1->SetNext(current2);
				current2 = next1;
				GetChromosome(i, 2).SetHeadLocus(next1);
				continue;
			  }
			}
			
			if ((current1->GetTransposon().GetLocation() > site) && (current2->GetTransposon().GetLocation() > site))
			  continue;
			
			if ((current1->GetTransposon().GetLocation() > site) && (current2->GetTransposon().GetLocation() < site))
			{
			  current2->SetNext(current1);
			  current1 = next2;
			  GetChromosome(i, 1).SetHeadLocus(next2);
			  continue;
			}
			
			if ((current1->GetTransposon().GetLocation() < site) && (current2->GetTransposon().GetLocation() > site))
			{
			  current1->SetNext(current2);
			  current2 = next1;
			  GetChromosome(i, 2).SetHeadLocus(next1);
			  continue;
			}
			
			// both current1 & current2 < site
			
			if (next1 == 0 && next2 == 0)
			  break;
			else
			{
				current1->SetNext(next2);
				current2->SetNext(next1);
			}
		}
		//std::cout << "Recombination for chromo: " << i << std::endl;
		
		current1 = 0;
		current2 = 0;
		next1 = 0;
		next2 = 0;
		site = 0.0;
		events = 0;
		length = 0;
		r = 0;
		sitesVector.clear();
		
	  } // for: numberOfChromosome

} // Recombination method


//Genome Genome::MakeGamete()
//{
//	
//	// for ploidy == 2
//	
//	Recombination();
//	Genome gamete(numberOfChromosomes, 1);
//	
//	int c = 0;
//	
//	for (int i=1; i <= numberOfChromosomes; i++)
//	{
//		if (rand.Uniform() < 0.5)
//		  c = 1;
//		else
//		  c = 2;
//		  
//		Locus * current = GetChromosome(i,c).GetHeadLocus();
//		while (current != 0)
//		{
//			gamete.GetChromosome(i,1).Insert(current->GetTransposon());
//			current = current->GetNext();
//		}
//	}
//	
//	return gamete;
//}

void Genome::ListGenomeSites() const
{
	for (int i=0; i < (numberOfChromosomes*ploidy); i++)
	  chromoVector.at(i).ListChromSites();
	std::cout << std::endl;
}

int Genome::GetRandomTransposonFamily(int outOf)
{
    //int teCount = GetGenomeTECount(); //don't need this with outOf, which is quicker
    int remaining = (int) ((rand.Uniform()*outOf)+1);//picks a random te from 1 to total TEs;
    
    if(outOf==0){
        return -1;//this is strictly an error flag
    }
    int curnum=1,curcopy=1;
    
    int chromTEcount = chromoVector.at((2*curnum)+curcopy-3).GetChromTECount();
    
    while(remaining>chromTEcount){
        remaining-=chromTEcount;
        if(curnum==Genome::numberOfChromosomes){
            curnum=1;
            curcopy=2;
        }else{
            curnum++;
        }
        chromTEcount=chromoVector.at((2*curnum)+curcopy-3).GetChromTECount();
        
    }
    
    int ok;
    Locus * pointer =chromoVector.at((2*curnum)+curcopy-3).GetHeadLocus();
    for(int i=0;i<remaining-1;i++){
        pointer = pointer->GetNext();
    }
    ok = pointer->GetTransposon().GetFamily();
    return ok;
    
}

void Genome::Insert(int chrom,int cop,Transposon te)
{
    GetChromosome(chrom, cop).Chromosome::Insert(te);
    totalTECount++;
    UpdateCensus(true,te.GetFamily());
}

void Genome::Insert(int spec,Transposon te)
{
    chromoVector.at(spec).Chromosome::Insert(te);
    totalTECount++;
    UpdateCensus(true,te.GetFamily());
}

void Genome::Delete(int spec, int position)
{
    int famForUpdate = chromoVector.at(spec).Chromosome::Delete(position);
    totalTECount--;
    UpdateCensus(false,famForUpdate);
}

void Genome::UpdateCensus(bool insert,int fam)
{
    if(insert){//insert a transposon
        FamilyCensus * current = headFamily;
        
        if(current ==0){//there are no TEs yet
            FamilyCensus * newHead = new FamilyCensus(fam);
            newHead->IncCount();
            headFamily = newHead;
            numberOfFamilies++;
            return;
        }
        
        if(current->GetFamily() > fam){//new family is the lowest index in the genome, new headFamily required.
            FamilyCensus * newHead = new FamilyCensus(fam);
            newHead->IncCount();
            newHead->SetNext(headFamily);
            headFamily=newHead;
            numberOfFamilies++;
            return;
        }
        
        
        FamilyCensus * next = current->GetNext();
        
        while((next!=0) && (next->GetFamily() <= fam)){
            current = next;
            next = current->GetNext();
            
            
        }//current is now the correct spot
        
        if(current->GetFamily()==fam){//there are already elements of this family, simply increment them
            current->IncCount();
            return;
        }
        if(next==0){//this is a new family with an unprecendented family index
            FamilyCensus * newFamily = new FamilyCensus(fam);
            newFamily->IncCount();
            current->SetNext(newFamily);
            numberOfFamilies++;
            return;
        }else{//this means that there is a higher indexed family (pointed to by next), and that it needs to be replaced in order by a newFamily
            FamilyCensus * newFamily = new FamilyCensus(fam);
            newFamily->IncCount();
            newFamily->SetNext(next);
            current->SetNext(newFamily);
            numberOfFamilies++;
            return;
        }
        
    }else{//delete a transposon
        
        FamilyCensus * current = headFamily;
        if(headFamily==0){
            std::cout<<"Deletion error\n";
            return;
        }
        
        
        if (current->GetFamily()==fam){
            current->DecCount();
            if(current->GetCount()==0){
                //need to reset head
                headFamily = headFamily->GetNext();
                current->SetNext(0);
                delete current;
                current=0;
                numberOfFamilies--;
                
                
            }//reset head
            return;
        }//headFamily updated
        
        while(current->GetNext()->GetFamily()!=fam){
            current = current->GetNext();
        }//at this point, current->GetNext() is the family that needs update
        
        FamilyCensus * next = current->GetNext();
        next->DecCount();
        if(next->GetCount()==0){
            current->SetNext(next->GetNext());
            next->SetNext(0);
            delete next;
            next=0;
            numberOfFamilies--;
        }
        return;
        
    }//delete
        
}

unsigned int Genome::GetNumberOfFamilies() const
{
    /*
    //this does not do what it is titled for debugging puposes as of Jun 23
    unsigned int tbr=0;
    FamilyCensus * current = headFamily;
    
    while(current!=0){
        tbr+=current->GetCount();
        current=current->GetNext();
    }
    return tbr;
     */
    return numberOfFamilies;
}

FamilyCensus * Genome::GetHeadFamily()
{
    return headFamily;
}

void Genome::SetHeadFamily(FamilyCensus * newHead)
{
    headFamily = newHead;
}

void Genome::ResetTotalTECount()
{
    totalTECount=0;
    numberOfFamilies=0;

}