#pragma once
#include <vector>
#include <unordered_map>
#include <hash_map>
#include <stack> 
#include <mpi.h>
#include "Neuron.h"

using namespace std;
class GrowingNetwork
{
public:
	GrowingNetwork(void);
	GrowingNetwork(bool optimalization);
	~GrowingNetwork(void);
	bool InitNetwork();
    bool InitNetworkMass(vector<double*> vahy,vector<vector<int>> neighnor,vector<int> ID);
	void SetAtributes(double e_w, double e_n, double alpha, double beta, int gama, int a_max, int maxNeuron, int MaxIteration);
	bool LoadConfig(string fileName);
	bool LoadTrainDate(string fileName);
	bool LoadNeuronMap(string configName,int IDRanku);
	void CheckVaule(int Source);
	void CreateSimilar(int Source, int Dest);
	void FindTwoNeuron(int counterDataSet, int *firstID,  int *secondID);
	void AddNeuron();
	void DecreaseErrorBeta();
	 bool StartTrainingMPI();
	 bool StartTraining(bool onlyOneCore);
	 bool SaveToGraphViz(string fileName);
	 bool SaveToGnuplot(string fileName);
	 bool SaveToGraphVizTestovani(string fileName);
	 double ErrorNetwork();
	 double ErrorNetworkCompareSOM();
	 void GenerateAnnulus();
	 void SaveForAgainTest(string path);
	 void ConvertToRank();
	 bool LoadWithoutSOMOnlySparsData(string nameOfTreningSet, int numberOfPartsLocal);
	 bool LoadWithoutSOM(string configName);
	 bool LoadForOneProces(string configName);
	 void SaveToGephiCSV(string path);
	  void SaveToGephiGdf(string path, int SizeChange);
	  void SaveInputDataToGephiGdf(string path,int  SizeChange);
	  void InitMPI();
	  void UpdateMPI();

	     int Rank;
		 int Size;
		 MPI_Comm comm;
		 		 bool debug;
private:
	        /// <summary>
        /// Faktor u?ení pro vít?zný neuron
        /// </summary>
          double e_w;
        /// <summary>
        /// Faktor u?ení pro sousedy vít?zného neuronu
        /// </summary>
          double e_n;
        /// <summary>
        /// Koeficient snížení akumulované chyby v neuronech p?i p?idávání nového neuronu
        /// </summary>
          double alpha;
        /// <summary>
        /// Koeficient snižování akumulované chyby v neuronech
        /// </summary>
          double beta;
        /// <summary>
        /// Interval vstupních vzor? pro p?idání nového neuronu
        /// </summary>
          int gama;
        /// <summary>
        /// Maximální stá?í hrany
        /// </summary>
         int a_max;
         int NeuronMaxCount;
        /// <summary>
        /// Po?et iterací p?i u?ení
        /// </summary>
         int MaxIteration;
    //     map<int, Neuron*> neuronField;
		  //unordered_map<int, Neuron*> neuronField;
		 Neuron ** neuronField;
		 int maxOfActualused;
		 stack<int> mystackOfFree;
		 int realNumberOfNeurons;

		 double** trainingSet;
		 int** trainingSetNotNull;

		 int maxDimension;


		 bool optimalization;
		 double minkovNumber;
		 int numberOfRecords;
		 int *numberOfDataInRows;
		 double *SumOfInput;
		 int numberofAdd;
		 int numberofRemove;
};

