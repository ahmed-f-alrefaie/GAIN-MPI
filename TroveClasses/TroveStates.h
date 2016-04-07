#pragma once

#include "../BaseClasses/States.h"




class TroveStates : public States{
	void ReadJGStates(int J,int G,int j);
public:
	TroveStates(Input & pInput);
	void ReadStates();
};
