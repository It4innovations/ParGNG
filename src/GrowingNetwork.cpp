#include "GrowingNetwork.h"
#include <iostream>
#include <fstream>
#include <string>
#include<sstream>


typedef std::pair<int, Neuron*> MyPair;
using namespace std;


template <typename T>
string NumberToString ( T Number )
{
	stringstream ss;
	ss << Number;
	return ss.str();
}

template <typename T>
T StringToNumber ( const string &Text )//Text not by const reference so that the function can be used with a 
{                               //character array as argument
	stringstream ss(Text);
	T result;
	return ss >> result ? result : 0;
}

void split(string &s1, string &s2, char delim){
    size_t i;

    i=s1.find_first_of(delim);
    s2.append(s1, i+1, s1.size());
    s1.erase(i, s1.size());
  }

GrowingNetwork::GrowingNetwork(void)
{
	trainingSet=NULL;
				maxDimension=0;
				numberofAdd=0;
				numberofRemove=0;
}
  GrowingNetwork::GrowingNetwork(bool optimalization1)
   {          
            gama = 200; 
            e_w = 0.05;
            e_n = 0.006;
            alpha = 0.5;
            beta = 0.0005;
            a_max = 30;
            NeuronMaxCount = 500;
            MaxIteration = 200000;
			trainingSet=NULL;
            optimalization = optimalization1;
			debug = true;
			minkovNumber=2;
			trainingSet=NULL;
			maxDimension=0;
  }

GrowingNetwork::~GrowingNetwork(void)
{
}
 bool  GrowingNetwork::InitNetwork()
 {
            if(trainingSet==NULL)
                return false;

			maxOfActualused=2;
			realNumberOfNeurons=2;
			for (int i = NeuronMaxCount-1; i > 1; i--)
			{
				mystackOfFree.push(i);
			}
							if(debug)
				{
					cout<<"BEFORE INIT"<<endl;
					cout<<gama<<endl<<e_w<<endl<<e_n<<endl<<alpha<<endl<<beta<<endl<<a_max<<endl<<NeuronMaxCount<<endl<<MaxIteration<<endl<<maxDimension<<endl;
				}
			neuronField= new Neuron*[NeuronMaxCount];
			neuronField[0]=new Neuron(maxDimension, a_max, e_w, e_n, alpha, beta, 0, 0, optimalization,minkovNumber,NeuronMaxCount);
			neuronField[1]=new Neuron(maxDimension, a_max, e_w, e_n, alpha, beta, 1, 1, optimalization,minkovNumber,NeuronMaxCount);

		//	neuronField.insert(MyPair(0, new Neuron(maxDimension, a_max, e_w, e_n, alpha, beta, 0, 0, optimalization,minkovNumber)));
		//	neuronField.insert(MyPair(1, new Neuron(maxDimension, a_max, e_w, e_n, alpha, beta, 1, 1, optimalization,minkovNumber)));

            neuronField[0]->CreateOrSetEdgeFS(1);
            neuronField[1]->CreateOrSetEdgeFS(0);
            return true;
 }
  bool  GrowingNetwork::InitNetworkMass(vector<double*> vahy,vector<vector<int>> neighnor,vector<int> ID)
 {
            if (trainingSet == NULL)
                return false;

			int numberOfNeurons = neighnor.size();
		//	NeuronMaxCount=numberOfNeurons;
			realNumberOfNeurons=0;
		//	maxOfActualused=0;
			neuronField= new Neuron*[NeuronMaxCount];
			for (int i = 0; i < NeuronMaxCount; i++)
			{
				neuronField[i]=NULL;
			}
            for (int i = 0; i < numberOfNeurons; i++)
            {
                Neuron *helpr = new Neuron(maxDimension, a_max, e_w, e_n, alpha, beta, ID[i], ID[i], optimalization,minkovNumber,NeuronMaxCount);
                helpr->vahy=vahy[i];
                
				int temp=(neighnor[i]).size();
                for (int j = 0; j <temp ; j++)
                {
                    helpr->CreateOrSetEdgeFS(neighnor[i][j]);
                }

				neuronField[ID[i]]= helpr;
				realNumberOfNeurons++;
			//	maxOfActualused++;
            }

      return true;
  }
  void GrowingNetwork::InitMPI()
  {
	  comm=MPI_COMM_WORLD;
	 MPI_Comm_size(comm,&Size);
		MPI_Comm_rank(comm,&Rank);
  }
    void GrowingNetwork::UpdateMPI()
  {

	 MPI_Comm_size(comm,&Size);
		MPI_Comm_rank(comm,&Rank);
  }
 void GrowingNetwork::SetAtributes(double e_w1, double e_n1, double alpha1, double beta1, int gama1, int a_max1, int maxNeuron1, int MaxIteration1)
   {
            gama = gama1;
            e_w = e_w1;
            e_n = e_n1;
            alpha = alpha1;
            beta = beta1;
            a_max = a_max1;
            NeuronMaxCount = maxNeuron1;
            MaxIteration = MaxIteration1;
    }

 bool GrowingNetwork::LoadConfig(string fileName)
        {
           	string line;
			  ifstream myfile (fileName);
			  if (myfile.is_open())
			  {
				getline (myfile,line);
				gama=atoi(line.c_str());

				getline (myfile,line);
				e_w=atof(line.c_str());

				getline (myfile,line);
				e_n=atof(line.c_str());

				getline (myfile,line);
				alpha=atof(line.c_str());

				getline (myfile,line);
				beta=atof(line.c_str());

				getline (myfile,line);
				a_max=atoi(line.c_str());

				getline (myfile,line);
				NeuronMaxCount=atoi(line.c_str());

				getline (myfile,line);
				MaxIteration=atoi(line.c_str());

				myfile.close();
				if(debug)
				{
					cout<<gama<<endl<<e_w<<endl<<e_n<<endl<<alpha<<endl<<beta<<endl<<a_max<<endl<<NeuronMaxCount<<endl<<MaxIteration<<endl<<maxDimension<<endl;
				}
			  }
				else
				{
					  cout << "Cannot find file"; 
					  return false;
			  }
			  return true;
        }

 bool GrowingNetwork::LoadTrainDate(string fileName)
        {
            numberOfRecords = 0;
            int numberOfDimensions = 0;
			string line;
			  ifstream myfile (fileName);
			  if (myfile.is_open())
			  {
				  getline (myfile,line);
				numberOfRecords=atoi(line.c_str());

				getline (myfile,line);
				numberOfDimensions=atoi(line.c_str());

				trainingSet=new double*[numberOfRecords];
				   for (int i = 0; i < numberOfRecords; i++)
				   {
					    trainingSet[i] = new double[numberOfDimensions];
					istringstream iss(line);
					do
					{
						string sub;
						string secondPart;
						iss >> sub;
						if(sub=="")
							break;
						split(sub,secondPart,':');

						trainingSet[i][StringToNumber<int>(sub)] = StringToNumber<double>(secondPart);
			
					} while (iss);
				   }
				myfile.close();
			  }
				else
				{
					  cout << "Cannot find: training date"; 
					  return false;
			  }
			  return true;

        }

 bool GrowingNetwork::LoadNeuronMap(string configName,int IDRanku)
        {
			 ifstream myfile1 (configName);
			if (!myfile1.is_open())
            {
                cout<<"Missing:"<<configName<<endl;
                return false;
            }
            numberOfRecords = 0;
            int numberOfDimensions = 0;
            string nameOfTreningSet;
            string nameOfWiner;
            string nameOfFileNeuron;

                string line;

				getline (myfile1,line);
				numberOfRecords=atoi(line.c_str());


				getline (myfile1,line);
				numberOfDimensions=atoi(line.c_str());

             //   trainingSet = new double[numberOfRecords][];
              //  nameOfTreningSet = sr.ReadLine();
				getline (myfile1,nameOfTreningSet);
            //    nameOfWiner = sr.ReadLine();
				getline (myfile1,nameOfWiner);

				getline (myfile1,line);
				nameOfFileNeuron = line + NumberToString<int> (IDRanku) + ".1Bbin";            

				myfile1.close();


            //Nacteni dat

			 ifstream myfile2 (nameOfTreningSet);
			if (!myfile2.is_open())
			{
                cout<<"Missing:" << nameOfTreningSet<<endl;
                return false;
            }
     
	maxDimension = 0;
	vector< vector<double> > tempData;
	vector< vector<int> > tempPosition;
	if(myfile2.is_open())
	{
		    while ( myfile2.good() )
			{
				vector<double> tempOneLineData;
				vector<int> tempOneLinePosition;
			  getline (myfile2,line);

			  istringstream iss(line);
					do
					{
						string sub;
						string secondPart;
						iss >> sub;
						if(sub=="")
							break;
						split(sub,secondPart,':');
						tempOneLinePosition.push_back(StringToNumber<int>(sub));
						tempOneLineData.push_back(StringToNumber<double>(secondPart));

						 int ook = StringToNumber<int>(sub);

						if (ook > maxDimension)
                           maxDimension=ook;

					//	cout << "Substring: " << sub << endl;
					} while (iss);
					if(tempOneLineData.size()!=0)
					{
					tempData.push_back(tempOneLineData);
					tempPosition.push_back(tempOneLinePosition);
					}
			}
	}
				myfile2.close();
			maxDimension++;



			ifstream myfile (nameOfWiner);
			if (!myfile.is_open())
            {
                cout<<"Missing:" << nameOfWiner<<endl;
                return false;
            }
			myfile.close();

			vector<int> prirazeni;
			ifstream r(nameOfWiner, ios::binary);
			int tempx;
			for (int i = 0; i < numberOfRecords; i++)
                {
					r.read((char*)&tempx, sizeof(tempx));
					prirazeni.push_back(tempx);
                }
			r.close();

            //vyrazeni nevyuzitych trenovacich dat

			for (int i = prirazeni.size()-1; i >=0; i--)
            {
                if (prirazeni[i] != IDRanku)
                {
					tempData.erase(tempData.begin()+i);
                    tempPosition.erase(tempPosition.begin()+i);
                }
            }
			if (tempData.size() == 0)
                return false;

				numberOfRecords=tempData.size();
				trainingSet=new double*[numberOfRecords];
				trainingSetNotNull=new int*[numberOfRecords];
				numberOfDataInRows = new int[numberOfRecords];
				SumOfInput=new double[numberOfRecords];
				for (int i = 0; i < numberOfRecords; i++)
				{
					int tempNumber=tempData[i].size();
					trainingSet[i]=new double[tempNumber];
					trainingSetNotNull[i]=new int[tempNumber];
					numberOfDataInRows[i]=tempNumber;
					double sum=0;
					for (int j = 0; j < tempNumber; j++)
					{
						trainingSet[i][j]=tempData[i][j];
						trainingSetNotNull[i][j]=tempPosition[i][j];

						sum+=tempData[i][j]*tempData[i][j];
					}
					SumOfInput[i]=sum;
				}


            //nainicializuji sit
            InitNetwork();

         //   neuronField[0].vahy = vahy.ToArray();
            CheckVaule(0);
            CreateSimilar(0, 1);
            return true;
	
	}
	void GrowingNetwork::CheckVaule(int Source)
	{
            bool result = false;
			for (int i = 0; i < maxDimension; i++)
            {
                if (neuronField[Source]->vahy[i] < 0.1)
                    result = false;
            }
            if (!result)
            {               
                //    neuronField[Source].vahy = trainingSet[0];
                for (int i = 0; i < numberOfDataInRows[0]; i++)
                {
                    neuronField[Source]->vahy[trainingSetNotNull[0][i]] = trainingSet[0][i];
                }
            }            
	}

         void GrowingNetwork::CreateSimilar(int Source, int Dest)
        {
            srand(Dest);
            for (int i = 0; i < maxDimension; i++)
            {
                double help = (((double)rand() / (double)RAND_MAX)/100);
                neuronField[Dest]->vahy[i] = neuronField[Source]->vahy[i] + help;
            }
        }
 void GrowingNetwork::FindTwoNeuron(int counterDataSet,int *firstID, int* secondID)
        {

         //   ICollection<Neuron> icoll = neuronField.Values;
		//	std::map<char,Neuron *>::iterator it;

            int index = 0;
            double fiRS = 0;
            int pozfID = -1;
            double secRS = 0;
            int pozsID = -1;

            double result = 0;
      //      foreach (Neuron s in icoll)
			Neuron *s;
			//for (std::unordered_map<int,Neuron *>::iterator it=neuronField.begin(); it!=neuronField.end(); ++it)
			// for ( auto it = neuronField.begin(); it != neuronField.end(); ++it )
			for(int i=0;i<maxOfActualused;i++)
            {
				if(neuronField[i]==NULL)
					continue;
				s=neuronField[i];
                s->input = trainingSet[counterDataSet];
                s->inputNonZero = trainingSetNotNull[counterDataSet];
				s->numberOfseznamNenul=numberOfDataInRows[counterDataSet];
				s->VstupY=SumOfInput[counterDataSet];
                result = s->Distance();
                if (index == 0)
                {
                    fiRS = result;
                    pozfID = s->ID;
                }else
                    if (index == 1)
                    {
                        if (fiRS < result)
                        {
                            secRS = result;
                            pozsID = s->ID;
                        }
                        else
                        {
                            secRS = fiRS;
                            pozsID = pozfID;

                            fiRS = result;
                            pozfID = s->ID;
                        }
                    }
                    else
                    {
                        if ((fiRS < result) && (result < secRS))
                        {
                            secRS = result;
                            pozsID = s->ID;
                        }
                        if (result < fiRS)
                        {
                            secRS = fiRS;
                            pozsID = pozfID;

                            fiRS = result;
                            pozfID = s->ID;
                        }
                    }

                index++;
            }
            *firstID = pozfID;
            *secondID = pozsID;
        }

        void GrowingNetwork::AddNeuron()
        {
			numberofAdd++;
          //  ICollection<Neuron> icoll = neuronField.Values;
            double maxError = -1;
            int IDError = 0;
            int IDrecord = 0;
            //vyhledam nejvetsi chybu
			Neuron *s;
		//	for (std::map<int,Neuron *>::iterator it=neuronField.begin(); it!=neuronField.end(); ++it)
			for(int i=0;i<maxOfActualused;i++)
            {
				if(neuronField[i]==NULL)
					continue;
				s=neuronField[i];
                if (s->error > maxError)
                {
                    maxError = s->error;
                    IDError = s->ID;
                }
             //   if (s->ID > IDrecord)
            //        IDrecord = s->ID;
            }
         //   IDrecord++;
			IDrecord=mystackOfFree.top();
			mystackOfFree.pop();

            //vyhledam chybu u sousedu
            double maxErrorSecond = -1;
            int IDErrorSecond = 0;
          //  ICollection<int> listOfNeighbor = neuronField[IDError]->neighbor.Keys;

           // foreach (int s in listOfNeighbor)
		//	for (std::map<int,int>::iterator it=neuronField[IDError]->neighbor.begin(); it!=neuronField[IDError]->neighbor.end(); ++it)


			 //for ( auto it = neuronField[IDError]->neighbor.begin(); it != neuronField[IDError]->neighbor.end(); ++it )
    //        {
				//if (neuronField[it->first]->error > maxErrorSecond)
    //            {
    //                maxErrorSecond = neuronField[it->first]->error;
    //                IDErrorSecond = neuronField[it->first]->ID;
    //            }
    //        }

			LinedList*		tempLinked=&(neuronField[IDError]->listOfNeighbor);
				for (data *it= tempLinked->head;it!=NULL;it=it->next)
				{
				if (neuronField[it->ID]->error > maxErrorSecond)
                {
                    maxErrorSecond = neuronField[it->ID]->error;
                    IDErrorSecond = neuronField[it->ID]->ID;
                }
            }
			//Neuron *temp=neuronField[IDError]->listOfNeighbor.FindBest(neuronField);
			//maxErrorSecond=temp->error;
			//IDErrorSecond=temp->ID;

            if(debug)
				cout<<realNumberOfNeurons<<":"<<NeuronMaxCount<<".Vytvarim novy neuron:" << IDrecord << " s vaznabou na:" << IDError << " a " << IDErrorSecond<<endl;
            //vytvor uzel
       //     Random rand = new Random(50);

			Neuron *novy = new Neuron(maxDimension, a_max, e_w, e_n, alpha, beta, IDrecord, IDrecord*maxOfActualused, optimalization,minkovNumber,NeuronMaxCount);
            novy->SetWeight(neuronField[IDError]->vahy,neuronField[IDErrorSecond]->vahy);
			neuronField[IDrecord]=novy;

            //nastav hrany
            neuronField[IDError]->CreateOrSetEdgeFS(IDrecord);
            neuronField[IDrecord]->CreateOrSetEdgeFS(IDError);

            neuronField[IDErrorSecond]->CreateOrSetEdgeFS(IDrecord);
            neuronField[IDrecord]->CreateOrSetEdgeFS(IDErrorSecond);

            neuronField[IDErrorSecond]->DestroyEdge(IDError);
            neuronField[IDError]->DestroyEdge(IDErrorSecond);

            //Zmenseni chyby
            neuronField[IDError]->DecAplhaError();
            neuronField[IDErrorSecond]->DecAplhaError();
            neuronField[IDrecord]->error = neuronField[IDError]->error;

			realNumberOfNeurons++;
			if(IDrecord==maxOfActualused)
				maxOfActualused++;

        }

        void GrowingNetwork:: DecreaseErrorBeta()
        {
		//	for (std::unordered_map<int,Neuron *>::iterator it=neuronField.begin(); it!=neuronField.end(); ++it)
			 //for ( auto it = neuronField.begin(); it != neuronField.end(); ++it )
           	for(int i=0;i<maxOfActualused;i++)
            {
				if(neuronField[i]==NULL)
					continue;
					neuronField[i]->DecBetaError();
		    }
		}


        bool GrowingNetwork::StartTrainingMPI()
        {
            double localError = 0;
            if (trainingSet == NULL)
                return false;

            MaxIteration *= numberOfRecords;

            int counterDataSet = 0;
            int numberDataOfData = numberOfRecords;
            if (debug)
            {
                cout<<Rank << "\t number of record:" << numberDataOfData<<endl;
            }
            //   return false;
            if (numberDataOfData == 1)
            {
                for (int i = 0; i < numberOfDataInRows[0]; i++)
                {

                    neuronField[0]->vahy[trainingSetNotNull[0][i]] = trainingSet[0][i];
                    neuronField[1]->vahy[trainingSetNotNull[0][i]] = trainingSet[0][i];
                }
                //neuronField[0].vahy = trainingSet[0];
                //neuronField[1].vahy = trainingSet[0];
                neuronField[1]->vahy[0] += 0.01;
                if (debug)
                {
					cout<<Rank << "\tOnly one record, bye bye"<<endl;
                }
                
                return true;
            }
            int firstWin = 0, secondWin = 0;
			LinedList *tempLinked;
            for (int iterace = 0; iterace < MaxIteration; iterace++)
            {
                //Najit 2 uzly co nejmensi
                FindTwoNeuron(counterDataSet, & firstWin, &secondWin);

                //update chyby a vahy u viteze 
                neuronField[firstWin]->UpgradeError();
                neuronField[firstWin]->UpdateWeightsWiners();



              //  int[] icoll = neuronField[firstWin].neighbor.Keys.ToArray();

                //zvyseni veku hran od viteze
                neuronField[firstWin]->IncEdgeToNeighbor();

                //Vytvor spojeni mezi first a second
                neuronField[firstWin]->CreateOrSetEdgeFS(secondWin);
                neuronField[secondWin]->CreateOrSetEdgeFS(firstWin);

                //Aktualizuju vahy, zvys vek hran a odstran nezapojene uzly
				int kvp=0;
            //    for (std::map<int,int>::iterator it=neuronField[firstWin]->neighbor.begin(); it!=neuronField[firstWin]->neighbor.end(); ++it)
          //    for ( auto it = neuronField[firstWin]->neighbor.begin(); it != neuronField[firstWin]->neighbor.end(); ++it )
				tempLinked=&(neuronField[firstWin]->listOfNeighbor);
				for (data *it= tempLinked->head;it!=NULL;it=it->next)
				{
					kvp=it->ID;
					neuronField[kvp]->UpdateWeightsNeighbor();
                    if (kvp == secondWin)
                        continue;
                    if (!neuronField[kvp]->IncEdgeToWinner(firstWin))
                    {
						delete neuronField[kvp];
						neuronField[kvp]=NULL;
						mystackOfFree.push(kvp);
						realNumberOfNeurons--;
                        if (debug)
                            cout<<"Delete neuron:" << kvp<<endl;
                    }


                }


                //Novy uzel
			  if ((((iterace + 1) %gama) == 0) && (realNumberOfNeurons < NeuronMaxCount))
                    //   if ((iterace == gama) || (iterace == 2 * gama) || (iterace == 3 * gama) || (iterace == 4 * gama))
                        //    if (neuronField.Count <= (trainingSet.Length*10))
                        AddNeuron();



                //Snizeni chyby o beta
                DecreaseErrorBeta();



                counterDataSet++;
                if ((counterDataSet % numberDataOfData) == 0)
                {
                    counterDataSet = 0;
                }
                //if (debug)
                if (((iterace + 1) % 500) == 0)
                {
                    cout<<Rank << "\t" << ErrorNetwork() << "\t" << iterace<<endl;
                }
            }


            return true;
        }

    bool GrowingNetwork::StartTraining(bool onlyOneCore)
        {
            double localError = 0;

			            int firstWin=0, secondWin=0;
			int kvp=0;
			int counterEpoch=0;
			int counterDataSet = 0;
			int numberDataOfData=0;
            if (trainingSet == NULL)
                goto endFail;
                //return false;

            if (onlyOneCore)
            {
                MaxIteration = 10*numberOfRecords;
            }
            else
            {
                MaxIteration *= numberOfRecords;
            }
 
             numberDataOfData= numberOfRecords;
            if (debug)
            {
                cout<<Rank << "\t number of record:" << numberDataOfData<<endl;
            }
            if (numberDataOfData == 0)
                goto endFail;
             //   return false;
            if (numberDataOfData == 1)
            {
                for (int i = 0; i < numberOfDataInRows[0]; i++)
                {

                neuronField[0]->vahy[trainingSetNotNull[0][i]] = trainingSet[0][i];
                neuronField[1]->vahy[trainingSetNotNull[0][i]] = trainingSet[0][i];    
                }
                //neuronField[0].vahy = trainingSet[0];
                //neuronField[1].vahy = trainingSet[0];
                neuronField[1]->vahy[0] += 0.01;
                if (debug)
                {
					cout<<Rank<<"\tOnly one record, bye bye"<<endl;
                }
                goto end;
                //return true;
            }

            for (int iterace = 0; iterace < MaxIteration; iterace++)
            {
                //Najit 2 uzly co nejmensi
                FindTwoNeuron(counterDataSet,& firstWin,& secondWin);
                
                //update chyby a vahy u viteze 
                neuronField[firstWin]->UpgradeError();
                neuronField[firstWin]->UpdateWeightsWiners();

 

              //  int[] icoll = neuronField[firstWin].neighbor.Keys.ToArray();

				//Vytvor spojeni mezi first a second
                neuronField[firstWin]->CreateOrSetEdgeFS(secondWin);
                neuronField[secondWin]->CreateOrSetEdgeFS(firstWin);

                //Aktualizuju vahy, zvys vek hran a odstran nezapojene uzly

                			
           //     for (std::map<int,int>::iterator it=neuronField[firstWin]->neighbor.begin(); it!=neuronField[firstWin]->neighbor.end(); ++it)
           //   for ( auto it = neuronField[firstWin]->neighbor.begin(); it != neuronField[firstWin]->neighbor.end(); ++it )
				LinedList *	tempLinked=&(neuronField[firstWin]->listOfNeighbor);
				for (data *it= tempLinked->head;it!=NULL;it=it->next)
				{
					kvp=it->ID;
                
                    neuronField[kvp]->UpdateWeightsNeighbor();
                    if (kvp == secondWin)
                        continue;
                    if (!neuronField[kvp]->IncEdgeToWinner(firstWin))
                    {
						delete neuronField[kvp];
						neuronField[kvp]=NULL;
						mystackOfFree.push(kvp);
						realNumberOfNeurons--;
                            if (debug)
								cout<<"Delete neuron:" << kvp<<endl;
							numberofRemove++;
                    }


				}

				//zvyseni veku hran od viteze
                neuronField[firstWin]->IncEdgeToNeighbor();

                //Vytvor spojeni mezi first a second
                neuronField[firstWin]->CreateOrSetEdgeFS(secondWin);
                neuronField[secondWin]->CreateOrSetEdgeFS(firstWin);
			//	neuronField[secondWin]->UpdateWeightsNeighbor();
           

                //Novy uzel
				if (((iterace + 1) % gama == 0) && (realNumberOfNeurons < NeuronMaxCount))
                    //   if ((iterace == gama) || (iterace == 2 * gama) || (iterace == 3 * gama) || (iterace == 4 * gama))
                    if (onlyOneCore == false)
                    //    if (neuronField.Count <= (trainingSet.Length*10))
                            AddNeuron();



                //Snizeni chyby o beta
                DecreaseErrorBeta();



                counterDataSet++;
                if ((counterDataSet % numberDataOfData) == 0)
                {                    
                    counterDataSet = 0;
					counterEpoch++;
                }
                if (debug)
                if (((iterace+1) % 500) == 0)
                {
					if(counterEpoch==582)
						cout<<"OK"<<endl;
					cout<<Rank<<"\t"<<ErrorNetwork()<<"\t"<<counterEpoch<<endl;
					// cout<<Rank<<"\t"<<0<<"\t"<<iterace<<endl;
                }
            }
end:            if (debug)
            {                
               cout<<"Computing is done, rank:"<<Rank<<endl;                
            }
            localError = ErrorNetworkCompareSOM();
endFail: if (debug)
            {
                cout<<"Local error rank:" << Rank << " \t S with extension:"<<localError<<endl;
            }

		 double globalError=0;

		 MPI_Reduce(&localError,&globalError,1,MPI_DOUBLE,MPI_SUM ,0,comm);
         //   double globalError = comm.Reduce(localError, Operation<double>.Add, 0);



            if (Rank == 0)
				{
					cout<<"Global Error is " << globalError<<endl;
					cout<<"Number of add neruons:"<<numberofAdd<<". Number of remove:"<<numberofRemove<<endl;
			}
            return true;
        }

         bool GrowingNetwork::SaveToGraphViz(string fileName)
        {
			//using (StreamWriter sw = File.CreateText(fileName))
			// {
			ofstream myfile;
			myfile.open(fileName); 
			myfile<<"graph G {"<<endl;

			Neuron *s;
			int s1=0;
		//	for (std::map<int,Neuron *>::iterator it=neuronField.begin(); it!=neuronField.end(); ++it)
			for(int i=0;i<maxOfActualused;i++)
            {
				if(neuronField[i]==NULL)
					continue;
				s=neuronField[i];

				int actualID = s->ID;    			
				//for (std::map<int,int>::iterator it=s->neighbor.begin(); it!=s->neighbor.end(); ++it)
			//	 for ( auto it = s->neighbor.begin(); it != s->neighbor.end(); ++it )
				//{
					//s1=it->first;
								LinedList *	tempLinked=&(s->listOfNeighbor);
				for (data *it= tempLinked->head;it!=NULL;it=it->next)
				{
					s1=it->ID;
					if(actualID<s1)
						myfile<<actualID<<" -- "<<s1<<";"<<endl;
				}
				//sw.WriteLine("Dnes je: {0}", DateTime.Now);
			}
			myfile<<"}"<<endl;
			myfile.close();

			return true;
        }

         bool GrowingNetwork::SaveToGnuplot(string fileName)
        {
			MPI_Status stat;
            int prijem=0;
            if (Rank != 0)
                //comm.Receive(comm.Rank - 1, 10, out prijem);
				MPI_Recv(&prijem,1,MPI_INT,Rank - 1,10, comm,&stat);
            else
            {
				std::ofstream outfile (fileName);
				outfile.close();
            }

           // Console.WriteLine("Zapisuje:" + comm.Rank + " Nazev souboru:" + fileName + " Pocet Dat:" + neuronField.Values.Count);
			cout<<"Zapisuje:" <<Rank << " Nazev souboru:" << fileName << " Pocet Dat:" << realNumberOfNeurons<<endl;
           // StreamWriter sw = File.AppendText(fileName);
			ofstream myfile;
			myfile.open (fileName,ios::app);

			Neuron *s;
			//for (std::map<int,Neuron *>::iterator it=neuronField.begin(); it!=neuronField.end(); ++it)
			for(int j=0;j<maxOfActualused;j++)
            {
				if(neuronField[j]==NULL)
					continue;
				s=neuronField[j];
				int WeightCount = maxDimension;

                for (int i = 0; i < WeightCount; i++)
                {
					myfile<<s->vahy[i] << "\t";
                }
				myfile<<endl;
                //sw.WriteLine("Dnes je: {0}", DateTime.Now);
            }
            myfile.close();
            if (Rank != Size - 1)
			{
				int x=10;
				MPI_Send(&x,1,MPI_INT,Rank + 1,10,comm);
			}
			return true;
        }

         bool GrowingNetwork::SaveToGraphVizTestovani(string fileName)
        {


        			ofstream myfile;
			myfile.open(fileName); 
                
               myfile<<" digraph g {"<<endl;
               myfile<<"graph ["<<endl;
               myfile<<"rankdir = \"LR\""<<endl;
               myfile<<"];"<<endl;
               myfile<<"node ["<<endl;
               myfile<<"fontsize = \"16\""<<endl;
               myfile<<"shape = \"ellipse\""<<endl;
               myfile<<"];"<<endl;
               myfile<<"edge ["<<endl;
               myfile<<"];"<<endl;


			   int numberOfWight= maxDimension;
                int counter = 0;

			Neuron *s;
		//	for (std::map<int,Neuron *>::iterator it=neuronField.begin(); it!=neuronField.end(); ++it)
			for(int j=0;j<maxOfActualused;j++)
            {
				if(neuronField[j]==NULL)
					continue;
				s=neuronField[j];
                    myfile<<"\"node"<<s->ID<<"\" ["<<endl;
                     myfile<<"label = \" <f0> "<<s->ID<<"|";
                    for (int i = 0; i < numberOfWight; i++)
                    {
                         myfile<<s->vahy[i];
                        if (i == (numberOfWight - 1))
                            break;
                         myfile<<"|";
                    }
                    myfile<<"\""<<endl;
                   
                    myfile<<"shape = \"record\""<<endl;
                    myfile<<"];"<<endl;

					int s1=0;
                   int actualID = s->ID;    			
				//for (std::map<int,int>::iterator it=s->neighbor.begin(); it!=s->neighbor.end(); ++it)
			//for ( auto it = s->neighbor.begin(); it !=s->neighbor.end(); ++it )
			//	{
			//		s1=it->first;
				   	LinedList *	tempLinked=&(s->listOfNeighbor);
				for (data *it= tempLinked->head;it!=NULL;it=it->next)
				{
					s1=it->ID;

                        if (actualID < s1)
                        {
                            myfile<<"\"node"<<actualID<<"\":f0 -> \"node"<<s1<<"\":f0 ["<<endl;
                            myfile<<"id = "<<counter++<<endl;
                            myfile<<"];"<<endl;

                           // sw.WriteLine("{0} -- {1};", actualID, s1);
                        }
                    }
                    //sw.WriteLine("Dnes je: {0}", DateTime.Now);
                }
                myfile<<"}"<<endl;
            

            return true;
        }


         double GrowingNetwork::ErrorNetwork()
        {
            double result=0;
			Neuron *s;
		//	for (std::map<int,Neuron *>::iterator it=neuronField.begin(); it!=neuronField.end(); ++it)
			int counter=0;
			for(int j=0;j<maxOfActualused;j++)
            {
				if(neuronField[j]==NULL)
					continue;
				s=neuronField[j];
                result += s->error;
				counter++;
            }
			cout<<"Counter:"<<counter<<". RealNumberofNeurons:"<<realNumberOfNeurons<<".Results:"<<result<<endl;
			//return result / realNumberOfNeurons;
			return result/counter;

        }
         double  GrowingNetwork::ErrorNetworkCompareSOM()
        {
            double result = 0,min,pom;
          //  ICollection<Neuron> icoll = neuronField.Values;
			Neuron *s;
			bool first=true;
			for (int i = 0; i < numberOfRecords; i++)
            {
				min = 1000000;
				first=true;
		//	for (std::map<int,Neuron *>::iterator it=neuronField.begin(); it!=neuronField.end(); ++it)
			for(int j=0;j<maxOfActualused;j++)
            {
				if(neuronField[j]==NULL)
					continue;
				s=neuronField[j];
                    s->input = trainingSet[i];
                    s->inputNonZero = trainingSetNotNull[i];
					s->VstupY=SumOfInput[i];
					s->numberOfseznamNenul=numberOfDataInRows[i];
                    pom = s->Distance();

					if(first)
					{
						min=pom;
						first=false;
					}
					else
					{
						if (min > pom)
							min = pom;
					}
                    
                }
                result += min;
            }
            return result/numberOfRecords;
        }

         void GrowingNetwork::GenerateAnnulus()
        {
            const int NumberOfPatterns = 10000;
            const double R1 = 1.0;
            const double R2 = 0.7;
            double r;
            double coords [2];
           
          //  Random rnd = new Random();
            double** ipl = new double*[NumberOfPatterns*2];
            for (int i = 0; i < NumberOfPatterns*2; i=i+2)
            {
                ipl[i] = new double[2];
                ipl[i+1] = new double[2];
                do
                {
                    coords[0] = (double)rand() / (double)RAND_MAX;
                    coords[1] = (double)rand() / (double)RAND_MAX;

                    r = coords[0] * coords[0] + coords[1] * coords[1];
                }
                while (r < R2 * R2 || r > R1 * R1);
               
                ipl[i][0]=coords[0];
                ipl[i][1] = coords[1];
                coords[0] *= -1;
                coords[1] *= -1;

                ipl[i+1][0] = coords[0];
                ipl[i+1][1] = coords[1];
            }

            trainingSet= ipl;
        }

         void GrowingNetwork::SaveForAgainTest(string path)
        {
            int prijem = 0;
            path += ".help";
          //  ConvertToRank();
			MPI_Status mpiStat;
            if (Rank != 0)
				MPI_Recv(&prijem,1,MPI_INT,Rank - 1, 10,comm,&mpiStat);
            else
            {
				remove( path.c_str() );
               		std::ofstream outfile (path);
					outfile.close();
            }

          //  BinaryWriter bw1 = new BinaryWriter(File.Open(path,FileMode.Append));
			ofstream myFile (path, ios::out | ios::app | ios::binary);

      //      ICollection<Neuron> icoll = neuronField.Values;
            if (Rank == 0)
               myFile.write(reinterpret_cast<const char*>(&Size), sizeof Size);// bw1.Write(comm.Size);                                               //pocet procesoru
            
			int tempCount=realNumberOfNeurons;
           myFile.write(reinterpret_cast<const char*>(&tempCount), sizeof tempCount);// bw1.Write(icoll.Count);                                                 //pocet neuronu
			Neuron *s;
			int c=0;
			int base=NeuronMaxCount*Rank;
			int tempNumber=0;
		//	for (std::map<int,Neuron *>::iterator it=neuronField.begin(); it!=neuronField.end(); ++it)

			for(int j=0;j<maxOfActualused;j++)
            {
				if(neuronField[j]==NULL)
					continue;
				s=neuronField[j];
				tempNumber=s->ID+base;
				myFile.write(reinterpret_cast<const char*>(&tempNumber), sizeof tempNumber);//bw1.Write(s.ID);                                                    //ukladam ID
               // ICollection<int> icole = s.neighbor.Keys;
				tempCount=s->listOfNeighbor.size();
                myFile.write(reinterpret_cast<const char*>(&tempCount), sizeof tempCount);// bw1.Write(icole.Count);                                             //zapisu pocet sousedu
  					
                   int actualID = s->ID;    			
			//	for (std::map<int,int>::iterator it=s->neighbor.begin(); it!=s->neighbor.end(); ++it)
				//    for ( auto it = s->neighbor.begin(); it !=s->neighbor.end(); ++it )
				//{
				//	c=it->first;
				   	LinedList *	tempLinked=&(s->listOfNeighbor);
				for (data *it= tempLinked->head;it!=NULL;it=it->next)
				{
					c=it->ID;
					tempNumber=c+base;
                   myFile.write(reinterpret_cast<const char*>(&tempNumber), sizeof tempNumber);//  bw1.Write(c);                                                   //zapisu id souseda
                }
				tempCount=-1;
                myFile.write(reinterpret_cast<const char*>(&tempCount), sizeof tempCount);// bw1.Write(-1);                                                      //pro test zapisu odelovac -1
				for (int i = 0; i <maxDimension; i++)
                {
					myFile.write(reinterpret_cast<const char*>(&s->vahy[i]), sizeof s->vahy[i]);// bw1.Write(s.vahy[i]);                                            //ukladam vahy tohoto neuronu
                }
            }


            myFile.close();             

			int temp=10;
            if (Rank != Size - 1)
				MPI_Send(&temp,1,MPI_INT,Rank + 1, 10,comm);
        }

         void GrowingNetwork::ConvertToRank()
        {
        //    ICollection<Neuron> icoll = neuronField.Values;
            bool isnOK=true;
			Neuron *s;
            int referen = 1000 * (Rank+1);
            while (isnOK)
            {
                referen *= 10;
                isnOK=false;
			
			int c=0;
			for(int j=0;j<maxOfActualused;j++)
            {
				if(neuronField[j]==NULL)
					continue;
				s=neuronField[j];
                    if (s->ID >= referen)
                        isnOK = true;
                }
            }

            int forDestroy = 0;
            int start = referen;
znova:      
			for(int j=0;j<maxOfActualused;j++)
			{
				if(neuronField[j]==NULL)
					continue;
				s=neuronField[j];
				if (s->ID >= start)
					continue;
				forDestroy = s->ID;
				int si=0;   			
				LinedList *	tempLinked=&(s->listOfNeighbor);
				for (data *it= tempLinked->head;it!=NULL;it=it->next)
				{
					si=it->ID;
					neuronField[si]->ChangeID(s->ID, referen);
				}
				s->ID = referen;
				//  neuronField.add(referen, s);
				neuronField[referen]=s;
				//neuronField.erase(forDestroy);
				delete neuronField[forDestroy];

				neuronField[forDestroy]=NULL;
				mystackOfFree.push(forDestroy);

				referen++;
				goto znova;
			}
        }
         bool GrowingNetwork::LoadWithoutSOMOnlySparsData(string nameOfTreningSet, int numberOfPartsLocal)
        {

		ifstream myfile2 (nameOfTreningSet);
			if (!myfile2.is_open())
			{
                cout<<"Missing:" << nameOfTreningSet<<endl;
                return false;
            }

				vector< vector<double> > tempData;
	vector< vector<int> > tempPosition;
	string line;
 		    while ( myfile2.good() )
			{
				vector<double> tempOneLineData;
				vector<int> tempOneLinePosition;
			  getline (myfile2,line);

			  istringstream iss(line);
					do
					{
						string sub;
						string secondPart;
						iss >> sub;
						if(sub=="")
							break;
						split(sub,secondPart,':');
						tempOneLinePosition.push_back(StringToNumber<int>(sub));
						tempOneLineData.push_back(StringToNumber<double>(secondPart));

						 int ook = StringToNumber<int>(sub);

						if (ook > maxDimension)
                           maxDimension=ook;

					//	cout << "Substring: " << sub << endl;
					} while (iss);
					if(tempOneLineData.size()!=0)
					{
					tempData.push_back(tempOneLineData);
					tempPosition.push_back(tempOneLinePosition);
					}
			}
	
				myfile2.close();
				if(debug)
				{
					cout<<"AFTER LOAD"<<endl;
					cout<<gama<<endl<<e_w<<endl<<e_n<<endl<<alpha<<endl<<beta<<endl<<a_max<<endl<<NeuronMaxCount<<endl<<MaxIteration<<endl<<maxDimension<<endl;
				}

            maxDimension++;   //je zduvodu pocitani od 0
				numberOfRecords=tempData.size();
				trainingSet=new double*[numberOfRecords];
				trainingSetNotNull=new int*[numberOfRecords];

				numberOfDataInRows = new int[numberOfRecords];
				SumOfInput=new double[numberOfRecords];
				for (int i = 0; i < numberOfRecords; i++)
				{
					int tempNumber=tempData[i].size();
					trainingSet[i]=new double[tempNumber];
					trainingSetNotNull[i]=new int[tempNumber];
					numberOfDataInRows[i]=tempNumber;
					double sum=0;
					for (int j = 0; j < tempNumber; j++)
					{
						trainingSet[i][j]=tempData[i][j];
						trainingSetNotNull[i][j]=tempPosition[i][j];
					

						sum+=tempData[i][j]*tempData[i][j];
					}
					SumOfInput[i]=sum;
				}

				        NeuronMaxCount *= numberOfPartsLocal;


            InitNetwork();

    

            return true;
        }

        bool GrowingNetwork::LoadWithoutSOM(string configName)
        {
			ifstream myfile1 (configName);
			if (!myfile1.is_open())
			{
                cout<<"Missing:" << configName<<endl;
                return false;
            }
             numberOfRecords = 0;
            int numberOfDimensions = 0;
            string nameOfTreningSet;
            string nameOfWiner;
            string nameOfFileNeuron = configName + ".help";
            string temp="";
			string line;

            int numberOfParts=0;

				getline (myfile1,line);
				numberOfRecords=atoi(line.c_str());

				getline (myfile1,line);
				numberOfDimensions=atoi(line.c_str());

				getline (myfile1,nameOfTreningSet);
				getline (myfile1,nameOfWiner);
				getline (myfile1,temp);
				getline (myfile1,line);
				numberOfParts=atoi(line.c_str());

            //Nacteni dat

          ifstream myfile2 (nameOfTreningSet);
			if (!myfile2.is_open())
			{
                cout<<"Missing:" << nameOfTreningSet<<endl;
                return false;
            }

				vector< vector<double> > tempData;
			vector< vector<int> > tempPosition;
		
 		    while ( myfile2.good() )
			{
				vector<double> tempOneLineData;
				vector<int> tempOneLinePosition;
			  getline (myfile2,line);

			  istringstream iss(line);
					do
					{
						string sub;
						string secondPart;
						iss >> sub;
						if(sub=="")
							break;
						split(sub,secondPart,':');
						tempOneLinePosition.push_back(StringToNumber<int>(sub));
						tempOneLineData.push_back(StringToNumber<double>(secondPart));

						 int ook = StringToNumber<int>(sub);

						if (ook > maxDimension)
                           maxDimension=ook;

					//	cout << "Substring: " << sub << endl;
					} while (iss);
					if(tempOneLineData.size()!=0)
					{
					tempData.push_back(tempOneLineData);
					tempPosition.push_back(tempOneLinePosition);
					}
			}
	
				myfile2.close();


            maxDimension++;   //je zduvodu pocitani od 0
				numberOfRecords=tempData.size();
				trainingSet=new double*[numberOfRecords];
				trainingSetNotNull=new int*[numberOfRecords];
				numberOfDataInRows = new int[numberOfRecords];
				SumOfInput=new double[numberOfRecords];
				for (int i = 0; i < numberOfRecords; i++)
				{
					int tempNumber=tempData[i].size();
					trainingSet[i]=new double[tempNumber];
					trainingSetNotNull[i]=new int[tempNumber];
					numberOfDataInRows[i]=tempNumber;
					double sum=0;
					for (int j = 0; j < tempNumber; j++)
					{
						trainingSet[i][j]=tempData[i][j];
						trainingSetNotNull[i][j]=tempPosition[i][j];

						sum+=tempData[i][j]*tempData[i][j];
					}
					SumOfInput[i]=sum;
				}
			NeuronMaxCount *= numberOfParts;
            InitNetwork();


            return true;
        }

        bool GrowingNetwork::LoadForOneProces(string configName)
        {
			ifstream myfile1 (configName);
			if (!myfile1.is_open())
			{
                cout<<"Missing:" << configName<<endl;
                return false;
            }
             numberOfRecords = 0;
            int numberOfDimensions = 0;
            string nameOfTreningSet;
            string nameOfWiner;
            string nameOfFileNeuron = configName + ".help";
            int numberOfParts=0;
			string line;
				getline (myfile1,line);
				numberOfRecords=atoi(line.c_str());

				getline (myfile1,line);
				numberOfDimensions=atoi(line.c_str());

				getline (myfile1,nameOfTreningSet);
				getline (myfile1,nameOfWiner);

				myfile1.close();

       ifstream myfile2 (nameOfTreningSet);
			if (!myfile2.is_open())
			{
                cout<<"Missing:" << nameOfTreningSet<<endl;
                return false;
            }

				vector< vector<double> > tempData;
			vector< vector<int> > tempPosition;
			
 		    while ( myfile2.good() )
			{
				vector<double> tempOneLineData;
				vector<int> tempOneLinePosition;
			  getline (myfile2,line);

			  istringstream iss(line);
					do
					{
						string sub;
						string secondPart;
						iss >> sub;
						if(sub=="")
							break;
						split(sub,secondPart,':');
						tempOneLinePosition.push_back(StringToNumber<int>(sub));
						tempOneLineData.push_back(StringToNumber<double>(secondPart));

						 int ook = StringToNumber<int>(sub);

						if (ook > maxDimension)
                           maxDimension=ook;

					//	cout << "Substring: " << sub << endl;
					} while (iss);
					if(tempOneLineData.size()!=0)
					{
					tempData.push_back(tempOneLineData);
					tempPosition.push_back(tempOneLinePosition);
					}
			}
	
				myfile2.close();


            maxDimension++;   //je zduvodu pocitani od 0
				numberOfRecords=tempData.size();
				trainingSet=new double*[numberOfRecords];
				trainingSetNotNull=new int*[numberOfRecords];
				numberOfDataInRows = new int[numberOfRecords];
				SumOfInput=new double[numberOfRecords];
				for (int i = 0; i < numberOfRecords; i++)
				{
					int tempNumber=tempData[i].size();
					trainingSet[i]=new double[tempNumber];
					trainingSetNotNull[i]=new int[tempNumber];
					numberOfDataInRows[i]=tempNumber;
					double sum=0;
					for (int j = 0; j < tempNumber; j++)
					{
						trainingSet[i][j]=tempData[i][j];
						trainingSetNotNull[i][j]=tempPosition[i][j];

						sum+=tempData[i][j]*tempData[i][j];
					}
					SumOfInput[i]=sum;
				}


				ifstream myfile3 (nameOfFileNeuron,ios::binary);
			if (!myfile3.is_open())
			{
				cout<<"Missing:" << nameOfFileNeuron<<endl;
                return false;
            }
            int count = 0;
            int ID = 0, numberOfNeighbor = 0;
			int tempInt=0;
			double tempDouble=0;

			vector<double*> vahy;
            vector<vector<int>> poolNeighbor;
            vector<int> IDpool ;

                int commSize =0;// sr1.ReadInt32();                                                     //pocet procesoru
				myfile3.read( reinterpret_cast<char*>( &commSize ), sizeof commSize );
				int maxID=0;

                for (int k = 0; k < commSize; k++)
                {
                   // count = sr1.ReadInt32();                                                        //pocet neuronu
					myfile3.read( reinterpret_cast<char*>( &count ), sizeof count );

					double* vahy_help;//=new double[maxDimension];
                    vector<int> neigh_Help;
                    for (int j = 0; j < count; j++)
                    {
						vahy_help=new double[maxDimension];
                      //  ID = sr1.ReadInt32();                                                       //ID
						myfile3.read( reinterpret_cast<char*>( &ID ), sizeof ID );
						if(maxID<ID)
							maxID=ID;
                       // numberOfNeighbor = sr1.ReadInt32();                                         //pocet sousedu
					myfile3.read( reinterpret_cast<char*>( &numberOfNeighbor ), sizeof numberOfNeighbor );
					neigh_Help.clear();
                        for (int i = 0; i < numberOfNeighbor; i++)
                        {
                           // neigh_Help.Add(sr1.ReadInt32());                                        //sousedi ID
							myfile3.read( reinterpret_cast<char*>( &tempInt ), sizeof tempInt );
							 neigh_Help.push_back(tempInt);
                        }
                        poolNeighbor.push_back(neigh_Help);

						myfile3.read( reinterpret_cast<char*>( &tempInt ), sizeof tempInt );
                        if (tempInt!= -1)
                        {
                           cout<<"Chyba pri Intu"<<endl;
                            return false;
                        }
                        //vahy_help = new List<double>();
                        //for (int i = 0; i < numberOfDimensions; i++)
                        //{

                        //    vahy_help.Add(sr1.ReadDouble());                                        //Vahy
                        //}
					
					//	myfile3.read( reinterpret_cast<char*>( &vahy_help ), (sizeof tempDouble)*numberOfDimensions );
						for (int i = 0; i <maxDimension; i++)
						{
							myfile3.read( reinterpret_cast<char*>( &tempDouble ), (sizeof tempDouble) );
							vahy_help[i]=tempDouble;
						}

						vahy.push_back(vahy_help);
						IDpool.push_back(ID);
                    }

                }
            
			NeuronMaxCount=NeuronMaxCount*commSize;
			maxOfActualused=maxID+1;
            InitNetworkMass(vahy, poolNeighbor, IDpool);

            return true;

        }


        void SaveToHonza(string path)
        {

        }

        void GrowingNetwork::SaveToGephiCSV(string path)
        {
			ofstream myfile;
			myfile.open(path); 
//                 neuronField
                string firstLine = "";
                unordered_map<int, int> map1;
              //  ICollection<int>keys = neuronField.Keys;
                int j = 0;
			int s1=0;
			//for (std::map<int,Neuron *>::iterator it=neuronField.begin(); it!=neuronField.end(); ++it)
			for(int i=0;i<maxOfActualused;i++)
            {
				if(neuronField[i]==NULL)
					continue;
				s1=i;
                    firstLine += ";n" + s1;
					map1.insert(std::pair<int,int>(s1, j));
                    j++;
                }
                myfile<<firstLine<<endl;

				int tempSize=map1.size();
                double *pole = new double[tempSize];
				
				Neuron *s;
			//for (std::map<int,Neuron *>::iterator it=neuronField.begin(); it!=neuronField.end(); ++it)
			for(int i=0;i<maxOfActualused;i++)
            {
				if(neuronField[i]==NULL)
					continue;
				s=neuronField[i];
                    string otherLine = "n" + s->ID;
                    for (int k = 0; k < tempSize; k++)
                    {
                        pole[k] = 0;
                    }

					int s1=0;
			//	for (std::map<int,int>::iterator it=s->neighbor.begin(); it!=s->neighbor.end(); ++it)
				 //for ( auto it = s->neighbor.begin(); it != s->neighbor.end(); ++it )
					//{
					//s1=it->first;
						LinedList *	tempLinked=&(s->listOfNeighbor);
				for (data *it= tempLinked->head;it!=NULL;it=it->next)
				{
					s1=it->ID;
                        pole[map1[s1]] = 1.0;
                    }

                    for (int k = 0; k < tempSize; k++)
                    {
                        otherLine += ";" + NumberToString<double>( pole[k]);
                    }
                    myfile<<otherLine<<endl;
                }
			myfile.close();
            
        }
       void GrowingNetwork::SaveToGephiGdf(string path, int SizeChange)
      {
			ofstream myfile;
			myfile.open(path); 

              myfile<<"nodedef> name VARCHAR,label VARCHAR,width DOUBLE,height DOUBLE,x DOUBLE,y DOUBLE,color VARCHAR"<<endl;

              unordered_map<int, int> map1;
          //    ICollection<Neuron> icoll = neuronField.Values;
              int temp = 0;
				Neuron *s;
		//	for (std::map<int,Neuron *>::iterator it=neuronField.begin(); it!=neuronField.end(); ++it)
			for(int i=0;i<maxOfActualused;i++)
            {
				if(neuronField[i]==NULL)
					continue;
				s=neuronField[i];
                 map1.insert(std::pair<int,int>(s->ID, temp));
				 myfile<<s->ID<<",\"n"<<s->ID<<"\",10.0,10.0,"<<s->vahy[0] * SizeChange<<","<<s->vahy[1] * SizeChange<<",'153.153.153'"<<endl;
             //     sw.WriteLine(String.Format("{0},\"n{1}\",10.0,10.0,{2},{3},'153.153.153'", s.ID, s.ID, s.vahy[0] * SizeChange, s.vahy[1] * SizeChange));
                  temp++;
              }

               myfile<<"edgedef> node1,node2,weight DOUBLE,directed BOOLEAN"<<endl;
              vector<int> tempMap;
			  int s1=0;
	//		for (std::map<int,Neuron *>::iterator it=neuronField.begin(); it!=neuronField.end(); ++it)
		  	for(int i=0;i<maxOfActualused;i++)
            {
				if(neuronField[i]==NULL)
					continue;
				s=neuronField[i];
				int actualID = s->ID;    			
		//		for (std::hash_map<int,int>::iterator it=s->neighbor.begin(); it!=s->neighbor.end(); ++it)
		////	for(int j=0;j<maxOfActualused;j++)
  //          {
		//		s1=it->first;
					LinedList *	tempLinked=&(s->listOfNeighbor);
				for (data *it= tempLinked->head;it!=NULL;it=it->next)
				{
					s1=it->ID;
						if(std::find(tempMap.begin(), tempMap.end(), s1)!=tempMap.end())
							myfile<<actualID<<","<<s1<<",1.0,true"<<endl;
                       // sw.WriteLine(String.Format("{0},{1},1.0,true",actualID,s1));
                    }
                    tempMap.push_back(actualID);
              }
              
          myfile.close();
      }
      void GrowingNetwork::SaveInputDataToGephiGdf(string path,int  SizeChange)
      {
         // if (trainingSet[0].Length > 2)
          //    return;
			ofstream myfile;
			myfile.open(path); 
			       myfile<<"nodedef> name VARCHAR,label VARCHAR,width DOUBLE,height DOUBLE,x DOUBLE,y DOUBLE,color VARCHAR"<<endl;
          //    sw.WriteLine("nodedef> name VARCHAR,label VARCHAR,width DOUBLE,height DOUBLE,x DOUBLE,y DOUBLE,color VARCHAR");
				   for (int i = 0; i < numberOfRecords; i++)
              {
				   myfile<<i<<",\"n"<<i<<"\",10.0,10.0,"<<trainingSet[i][0]*SizeChange<<","<< trainingSet[i][1]*SizeChange<<",'153.153.153'"<<endl;
                 // sw.WriteLine(String.Format("{0},\"n{1}\",10.0,10.0,{2},{3},'153.153.153'", i, i, trainingSet[i][0]*SizeChange, trainingSet[i][1]*SizeChange));
              }
            //  sw.WriteLine("edgedef> node1,node2,weight DOUBLE,directed BOOLEAN");
			       myfile<<"edgedef> node1,node2,weight DOUBLE,directed BOOLEAN"<<endl;
           myfile.close();
      }
    
