#include "Neuron.h"
#include <vector>
#include <iostream>
#include <cmath>

using namespace std;
Neuron::Neuron(void)
{
}

Neuron::Neuron(int countAtributes, int a_max1, double e_w1, double e_n1, double alpha1, double beta1, int ID1,int rand1,bool optimalization1,double Mink, int maxNumberOfNeurons)
{
	srand(rand1);
	e_n = e_n1;
	e_w = e_w1;
	alpha = alpha1;
	beta = beta1;
	a_max = a_max1;
	vahy = new double[countAtributes];
	MinkowskehoNumber=Mink;
	numberOfdimension=countAtributes;
	//neighbor = new  Dictionary<int,int>();
	error = 0;
	ID = ID1;
	ZERO = 0.0000000001;
	neighbor= new int[maxNumberOfNeurons];
	tempArray= new int*[maxNumberOfNeurons];


	//optim
	optimalization = optimalization1;
	seznamNenulAll = new bool[countAtributes];
	nonZeroCount = countAtributes;
	nonZeroVahy = new int[countAtributes];
	VahyX=0;
	///

	for (int i = 0; i < maxNumberOfNeurons; i++)
	{
		neighbor[i]=-1;
	}

	for (int i = 0; i < countAtributes; i++)
	{
		  vahy[i]=r2();
		  VahyX += vahy[i] * vahy[i];
		  nonZeroVahy[i] = i;
		  seznamNenulAll[i]=false;
	}
	nonZeroCount = countAtributes;
}

Neuron::~Neuron(void)
{
		delete vahy;
	delete nonZeroVahy;
	delete seznamNenulAll;
}
 double Neuron::Distance()
        {
			return Euklid2();

            int coutner = 0;
            resultDistance = 0;
            int count = numberOfdimension;
            double pom = 0;
            for (int i = 0; i < count; i++)
            {
                if (i == inputNonZero[coutner])
                {
                    pom = vahy[i] - input[coutner];
                    resultDistance += pom * pom;
                    coutner++;
                    if (coutner == numberOfseznamNenul)
                        coutner = 0;
                }
                else
                {
                    resultDistance += vahy[i] * vahy[i];
                }
            }
			resultDistance=sqrt(resultDistance);
			
            return resultDistance;
        }
double Neuron::Minkowskeho(void)
{
	if(MinkowskehoNumber==2)
	{
		return Euklid2();
	}else
	{
		    double suma = 0;
            double pom = 0;
            int coutner = 0;

            for (int i = 0; i < numberOfdimension; i++)
            {
				if (i == inputNonZero[coutner])
				 {
					 pom = input[coutner] - vahy[i];
					suma += pow(abs(pom),MinkowskehoNumber);
					coutner++;
					if (coutner == numberOfseznamNenul)
                        coutner = 0;
                }
                else
                {
					 suma += pow(vahy[i],MinkowskehoNumber);
				}
            }
			return pow(suma,1/MinkowskehoNumber);// sqrt(suma);
	}
}

double Neuron::Euklid2(void)
{
	 double suma = 0;
	  resultDistance = 0;
            int pocet = numberOfseznamNenul;

                for (int i = 0; i < pocet; i++)
                {
					suma -= 2 * input[i] * vahy[inputNonZero[i]];
                }
				
            suma += VahyX;
            suma += VstupY;					
	
			double testX2=0;
			for (int i = 0; i < numberOfdimension; i++)
				{
					testX2+=vahy[i]*vahy[i];
				}
			//if(testX2!=VahyX)
			//	cout<<"Error:Improve"<<endl;
			if(suma<0) 
			{
				suma=0;
 			//	cout<<"error"<<endl;
				//double testX=0;
				//double testY=0;
				//double testALL=0;

				//for (int i = 0; i < pocet; i++)
				//{
				//	testY+=input[i]*input[i];

				//}

				//for (int i = 0; i < numberOfdimension; i++)
				//{
				//	testX+=vahy[i]*vahy[i];
				//}

				//int coutner = 0;
				//int count = numberOfdimension;
				//double pom = 0;
				//for (int i = 0; i < count; i++)
				//{
				//	if (i == inputNonZero[coutner])
				//	{
				//		pom = vahy[i] - input[coutner];
				//		testALL += pom * pom;
				//		coutner++;
				//		if (coutner == numberOfseznamNenul)
				//			coutner = 0;
				//	}
				//	else
				//	{
				//		testALL += vahy[i] * vahy[i];
				//	}
				//}


			}
					resultDistance=sqrt(suma);
            return resultDistance;
	
}
double Neuron::r2(void)
{
	return (double)0.6+(((double)rand() / (double)RAND_MAX)/100) ;

}

        void Neuron::UpgradeError()
        {
            error += resultDistance;
        }

   void Neuron::UpdateWeights(double e)
        {
			int count = numberOfdimension;
            int coutner = 0;
			VahyX = 0;
            for (int i = 0; i < count; i++)
            {
                if (i == inputNonZero[coutner])
                {
                    vahy[i] += e * (input[coutner] - vahy[i]);					
                    coutner++;
                    if (coutner ==numberOfseznamNenul)
                        coutner = 0;
                }
                else
                {
                    vahy[i] += e * (0 - vahy[i]);
                }
				VahyX += vahy[i] * vahy[i];
             }
        }

         void Neuron::UpdateWeightsOpt(double e)
        {
			VahyX = 0;
			for (int i = 0; i < numberOfseznamNenul; i++)
            {
                seznamNenulAll[inputNonZero[i]] = true;
            }

            int tmpNonZero = 0;
            for (int i = 0; i < nonZeroCount; i++)
            {
                int index = nonZeroVahy[i];

                if (seznamNenulAll[index]) continue;

                double vahyI = vahy[index];

                vahyI -= e * vahyI;
                if (vahyI < ZERO)
                {
                    vahy[index] = 0;
                }
                else
                {
					//VahyX += vahyI * vahyI;
                    nonZeroVahy[tmpNonZero] = index;
                    tmpNonZero++;
                    vahy[index] = vahyI;
					VahyX +=vahy[index]*vahy[index];
                }
            }

			for (int i = 0; i < numberOfseznamNenul; i++)
            {
                int index = inputNonZero[i];
                double vahyI = vahy[index];
                vahyI += e * (input[i] - vahyI);
			//	VahyX += vahyI * vahyI;
                vahy[index] = vahyI;
				VahyX +=  vahy[index]  *  vahy[index] ;
                seznamNenulAll[index] = false;
                nonZeroVahy[tmpNonZero] = index;
                tmpNonZero++;
            }
            nonZeroCount = tmpNonZero;
			//double testX=0;
			//				for (int i = 0; i < numberOfdimension; i++)
			//	{
			//		testX+=vahy[i]*vahy[i];
			//	}

        }


         void Neuron::UpdateWeightsWiners()
        {
      //      if(!optimalization)
        //        UpdateWeights(e_w);
     //       else
             UpdateWeightsOpt(e_w);
            isWinner = true;
        }

         void Neuron::UpdateWeightsNeighbor()
        {
       //     if (!optimalization)
         //       UpdateWeights(e_n);
      //      else
               UpdateWeightsOpt(e_n);
        }
         void Neuron::IncEdgeToNeighbor()
        {
		//	vector<int> icoll;
		////	for(map<int,int>::iterator it = neighbor.begin(); it != neighbor.end(); ++it) {
		//	 for ( auto it = neighbor.begin(); it != neighbor.end(); ++it )
		//	 {
		//	  icoll.push_back(it->first);
		//	}
  //    
		//	for (std::vector<int>::iterator it = icoll.begin(); it != icoll.end(); ++it)
		//	{
		//		 if(neighbor[*it]==(a_max-1))
  //              {
		//			neighbor.erase(*it);
  //               //   Console.WriteLine("Odstranuji hranu od " + ID + " do " + kvp);
  //                  continue;
  //              }
  //              neighbor[*it]++;
		//	}

			//int tempSize=listOfNeighbor.size();

			//for (int i = tempSize-1; i >=0; i--)
			//{
			//	if(*listOfNeighbor[i]==(a_max-1))
			//	{
			//		*listOfNeighbor[i]=-1;
			//		listOfNeighbor.erase(listOfNeighbor.begin()+i);
			//	}else
			//	{
			//		(*listOfNeighbor[i])++;
			//	}
			//}
			
			listOfNeighbor.IncreaseByOne(a_max-1);


            //for each( int kvp in icoll)           
            //{               
            //    if(neighbor[kvp]==(a_max-1))
            //    {
            //        neighbor.Remove(kvp);
            //     //   Console.WriteLine("Odstranuji hranu od " + ID + " do " + kvp);
            //        continue;
            //    }
            //    neighbor[kvp]++;
            //}     

        }
         bool Neuron::IncEdgeToWinner(int IDWinner)
        {
            if (neighbor[IDWinner]== (a_max - 1))
                {
					//int tempSize=listOfNeighbor.size();
					//	for (int i = tempSize-1; i >=0; i--)
					//	{
					//		if(listOfNeighbor[i]==&neighbor[IDWinner])
					//		{
					//			*listOfNeighbor[i]=-1;
					//			listOfNeighbor.erase(listOfNeighbor.begin()+i);
					//			break;
					//		}
					//	}
					listOfNeighbor.Delete(&neighbor[IDWinner]);
					
					//cout<<"Odstranuji hranu od " << ID <<" do " << IDWinner<<endl;
						if (listOfNeighbor.count!= 0)
                        return true;
                    else
                        return false;
                }
            neighbor[IDWinner] ++;                        
                return true;            
        }

         void Neuron::CreateOrSetEdgeFS(int ID)
        {
            //int value;
            //bool test = false;
            //if (neighbor.TryGetValue(ID, out value))
            //{
            //    neighbor[ID] = 0;
            //}
            //else
            //{
            //    neighbor.Add(ID, 0);
            //    test = true;
            //}
			if (neighbor[ID]>=0)
			{
                neighbor[ID] = 0;
            }
            else
            {
				neighbor[ID] = 0;
			//	listOfNeighbor.push_back(&neighbor[ID]);
				listOfNeighbor.Add(&neighbor[ID],ID);
				//neighbor.insert( std::pair<int,int>(ID, 0));
            }
            if (isWinner)
            {
                isWinner = false;
                //if (test)
              //      Console.WriteLine("new edge:{0} -> {1}",this.ID,ID);
            }
        }
         void Neuron::DecBetaError()
        {
            error -= beta * error;
        }

         void Neuron::DecAplhaError()
        {
            error = alpha* error;
        }

         void Neuron::SetWeight(double* first, double* second)
        {
			VahyX=0;
			 for (int i = 0; i < numberOfdimension ; i++)
            {
                vahy[i]=(first[i]+second[i])/2;
				VahyX += vahy[i] * vahy[i];
            }
        }

         void Neuron::DestroyEdge(int ID)
        {
			//neighbor.erase(ID);

			//int tempSize=listOfNeighbor.size();
			//			for (int i = tempSize-1; i >=0; i--)
			//			{
			//				if(listOfNeighbor[i]==&neighbor[ID])
			//				{
			//					*listOfNeighbor[i]=-1;
			//					listOfNeighbor.erase(listOfNeighbor.begin()+i);
			//					break;
			//				}
			//			}
						listOfNeighbor.Delete(&neighbor[ID]);
        }

         void Neuron::ChangeID(int oldID, int newID)
        {
        //    neighbor.insert( std::pair<int,int>(newID, neighbor[oldID]));
        //    neighbor.erase(oldID);

			neighbor[newID] = neighbor[oldID];
			neighbor[oldID]=-1;
			listOfNeighbor.ChangeItem(&neighbor[oldID],&neighbor[newID],newID);

			/*	listOfNeighbor.push_back(&neighbor[newID]);

					int tempSize=listOfNeighbor.size();
						for (int i = tempSize-1; i >=0; i--)
						{
							if(listOfNeighbor[i]==&neighbor[oldID])
							{
								*listOfNeighbor[i]=-1;
								listOfNeighbor.erase(listOfNeighbor.begin()+i);
								break;
							}
						}	*/
        }