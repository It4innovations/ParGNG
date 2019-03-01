#pragma once
#include <unordered_map>
#include <hash_map>
#include<vector>
#include "LinedList.h"
class Neuron
{
public:
	Neuron(void);
    Neuron(int countAtributes, int a_max, double e_w, double e_n, double alpha, double beta, int ID,int rand,bool optimalization,double Mink, int maxNumberOfNeurons);
	~Neuron(void);

	double Distance();
	void UpgradeError();
	double Euklid2(void);
	double Minkowskeho(void);
	void SetWeight(double* first, double* second);
	void DestroyEdge(int ID);
	 void ChangeID(int oldID, int newID);
public:
	void UpdateWeights(double e);
	void UpdateWeightsOpt(double e);
	void UpdateWeightsWiners();
	void UpdateWeightsNeighbor();
	void IncEdgeToNeighbor();
	bool IncEdgeToWinner(int IDWinner);
	  void CreateOrSetEdgeFS(int ID);
	    void DecBetaError();
	      void DecAplhaError();	

		  double r2(void);

public: double *input;
         int *inputNonZero;
         int ID ;
         double error;        
      //  std::hash_map<int, int> neighbor;
		 int *neighbor;
		 LinedList listOfNeighbor;
		 int **tempArray;

         double* vahy;
		double VahyX;
		int numberOfdimension;
			double VstupY;
	int  numberOfseznamNenul;
	double MinkowskehoNumber;
private: double resultDistance;   

          double e_w;
          double e_n;
          double alpha;
          double beta;        
          int a_max;
		  int seed;
        bool  isWinner;
        //Optim
         bool optimalization;
         bool* seznamNenulAll;
         int* nonZeroVahy;
         int nonZeroCount;
         double ZERO;
};

