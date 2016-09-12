#ifndef FAMILYCENSUS_H
#define FAMILYCENSUS_H


class FamilyCensus
{

	public:
		FamilyCensus(int whichFam);
		~FamilyCensus();
		int GetFamily() const;
		int GetCount() const;
		FamilyCensus * GetNext();
		void SetNext(FamilyCensus* set);
		void IncCount();
		void DecCount();
	
	private:
	  int family;
	  int counter;
	  FamilyCensus * next;
	
	
};

#endif