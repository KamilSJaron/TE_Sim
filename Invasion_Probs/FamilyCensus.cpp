#include "FamilyCensus.h"
#include <iostream>


FamilyCensus::FamilyCensus(int whichFam):
counter(0),
family(whichFam)
{
	next=0;
}

FamilyCensus::~FamilyCensus()
{
    
	delete next;
	next = 0;
}

int FamilyCensus::GetFamily() const
{
	return family;
}

int FamilyCensus::GetCount() const
{
	return counter;
}

FamilyCensus* FamilyCensus::GetNext()
{
	return next;
}

void FamilyCensus::SetNext(FamilyCensus* set)
{
	next = set;
}

void FamilyCensus::IncCount()
{
	counter++;
}

void FamilyCensus::DecCount()
{
	counter--;
}