#include "ParameterParser.h"


ParameterParser::ParameterParser(void)
{
	numberOfPartsLocal=0;
	showInputData=false;
	optimalization=false;
	debug=false;
	SizeOfChange=1;
	StandardParalelization=false;
	withSom=false;
}


ParameterParser::~ParameterParser(void)
{
}


bool ParameterParser::Parse(int numberOfitems, char *charArray[])
{
	for (int i = 0; i < numberOfitems; i++)
	{
		if(strcmp(charArray[i],"-c")==0)
		{
			nameOfConfigFile =charArray[i+1];
			i++;
			continue;
		}
		if(strcmp(charArray[i],"-o")==0)
		{
			nameOfOutputFile =charArray[i+1];
			i++;
			continue;
		}
		if(strcmp(charArray[i],"-i")==0)
		{
			nameOfInputDataFile =charArray[i+1];
			i++;
			continue;
		}
				if(strcmp(charArray[i],"-o")==0)
		{
			nameOfOutputFile =charArray[i+1];
			i++;
			continue;
		}
		if(strcmp(charArray[i],"-showInput")==0)
		{
			showInputData =true;
			continue;
		}
		if(strcmp(charArray[i],"-n")==0)
		{
			numberOfPartsLocal=atoi(charArray[i+1]);
			if(numberOfPartsLocal<1)
				numberOfPartsLocal=0;
			i++;
			continue;
		}

		if(strcmp(charArray[i],"-debug")==0)
		{
			debug=true;

		}
		if(strcmp(charArray[i],"-optim")==0)
		{
			optimalization=true;

		}
		if(strcmp(charArray[i],"-withSom")==0)
		{
			withSom=true;

		}
		if(strcmp(charArray[i],"-s")==0)
		{
			SizeOfChange=atoi(charArray[i+1]);
			if(SizeOfChange<1)
				SizeOfChange=1;
			i++;
			continue;
		}
		if(strcmp(charArray[i],"-SP")==0)
		{
			StandardParalelization=true;

		}

	}
	return false;
}
