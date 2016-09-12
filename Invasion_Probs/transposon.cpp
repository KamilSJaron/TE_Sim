// *********************************************************************
// 
// transposon.cpp
// 
// Created by: Elie Dolgin, University of Edinburgh
//
// First started: March 11, 2005
// Last edited:
//
// *********************************************************************

#include "transposon.h"

Transposon::Transposon():
  location(0),
  effect(false),
  family(0)
  {}

Transposon::Transposon(int site, bool e, int f):
  location(site),
  effect(e),
  family(f)
  {}

Transposon::~Transposon()
  {}

int Transposon::GetLocation() const
{
	return location;
}

bool Transposon::GetEffect() const
{
	return effect;
}

int Transposon::GetFamily() const
{
    
    return family;
}