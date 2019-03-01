#pragma once
#include <string>
#include <cstring>
#include <stdlib.h>
class ParameterParser
{
public:
	ParameterParser(void);
	~ParameterParser(void);
	bool Parse(int numberOfitems, char *charArray[]);
	char *nameOfOutputFile;
	char *nameOfConfigFile;
	char* nameOfInputDataFile;
	int numberOfPartsLocal;
	bool showInputData;
	bool debug;
	bool optimalization;
	bool withSom;
	int SizeOfChange;
	bool StandardParalelization;
};

