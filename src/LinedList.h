#pragma once
#include <string.h>

struct data
{
	data* next;
	data* prev;
	int * value;
	int ID;
};
class LinedList
{
public:
	LinedList(void);
	~LinedList(void);
public:
	data *head;
	data *tail;
public:
	bool Add(int * val, int IDval);
	bool Delete(int * val);
	data *actual;
	int count;
	bool DeleteActual(void);
	void IncreaseByOne(int limit);
	void ChangeItem(int * old,int * newItem, int newID);

//	Neuron* FindBest(Neuron** neuronField);
	// iteration
	data* Reset();

	data* Next();
	int size();
};

