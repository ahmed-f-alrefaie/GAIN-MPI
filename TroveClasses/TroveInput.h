#pragma once
#include "../BaseClasses/Input.h"

#pragma once
class TroveInput : public Input {
private:
	void RemoveBrackets(char* string);
	void ReadIntensity();
	void HandleSymmetery();
	std::string PrepareString(std::string line);
	void ComputeIGammaPair();

public:
	TroveInput(): Input(){}
	void ReadInput(const char* filename);

};
