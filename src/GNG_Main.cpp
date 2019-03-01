#include <mpi.h>
#include <stdio.h>
#include <iostream>
#include<sstream>
#include "ParameterParser.h"
#include "GrowingNetwork.h"

using namespace std;
//On GCC compile under -std=c++0x -pthread to get it to work.
template <typename T>
string NumberToString ( T Number )
{
	stringstream ss;
	ss << Number;
	return ss.str();
}

int main(int argc, char *argv[])
{
	  MPI_Init(&argc,&argv);

	ParameterParser parameters;
	parameters.Parse(argc,argv);

	GrowingNetwork gng;
	gng.InitMPI();
	gng.LoadConfig(parameters.nameOfConfigFile);
	gng.debug=parameters.debug;

	 if (gng.Size == 1)  // --debug
   // if(false) 
	{
	//	 gng.LoadConfig(parameters.nameOfConfigFile);

                        bool result = false;
						if (parameters.withSom)
							result = gng.LoadWithoutSOM(parameters.nameOfInputDataFile);
                        else
                            result = gng.LoadWithoutSOMOnlySparsData(parameters.nameOfInputDataFile, parameters.numberOfPartsLocal);

                        if (result)
                        {
							 double Start = MPI_Wtime();

                            gng.StartTraining(false);

                            double End = MPI_Wtime();
							cout<<"Total time without saving is:" << End-Start<<endl;

							string temp(parameters.nameOfOutputFile);
							gng.SaveToGraphViz( temp+ " ONLYONE " + NumberToString(gng.Rank));
                           //         gng.SaveToGraphVizTestovani(String.Format("test{0}.txt",comm.Rank));
                            gng.SaveToGnuplot(temp+ "ONLYONE.gn");
                            gng.SaveToGephiGdf(temp + "ONLYONE.gdf", parameters.SizeOfChange);
                            if (parameters.showInputData)
                                gng.SaveInputDataToGephiGdf(temp + "InputDataOnlyOne.gdf", parameters.SizeOfChange);
						}
	 }
	 else
	 {
		 if (gng.LoadNeuronMap(parameters.nameOfInputDataFile, gng.Rank))  // --debug
                           {
							   cout<<"Start rank:"<<gng.Rank<<endl;
                          //  gng.comm = (Intracommunicator)comm.Split(1, comm.Rank); // --debug
							    MPI_Comm_split(MPI_COMM_WORLD, 1, gng.Rank, &gng.comm);
   						   
								gng.UpdateMPI();
                            double Start = MPI_Wtime();
                            gng.StartTraining(false);

  
                           cout<<"Finish process:" << gng.Rank<<endl;
                            gng.SaveForAgainTest(parameters.nameOfInputDataFile);


							MPI_Barrier(gng.comm);
                           
                            if (gng.Rank == 0)
                            {
                                            double End = MPI_Wtime();
							cout<<"The time of first phase is:" << End-Start<<endl;

                                GrowingNetwork gng1 = new GrowingNetwork(parameters.optimalization);
                               gng1.LoadConfig(parameters.nameOfConfigFile);

                                if (gng1.LoadForOneProces(parameters.nameOfInputDataFile))
                                {
									gng1.InitMPI();
									gng1.comm=MPI_COMM_SELF;
									gng1.UpdateMPI();
                                    //gng1.comm = Communicator.self;
                                   double Start2 = MPI_Wtime();
                                    gng1.StartTraining(true);


                                     double End = MPI_Wtime();
										cout<<"The time of second phase is:" << End-Start2<<endl;
									 
									cout<<"Total time is:" << End-Start<<endl;



									string temp(parameters.nameOfOutputFile);
                                    gng1.SaveToGraphViz( temp+ " PARALLEL " + NumberToString(gng.Rank));
                                    //       gng.SaveToGraphVizTestovani(String.Format("test{0}.txt",comm.Rank));
                                    gng1.SaveToGnuplot(temp + ".gn");
                                    gng1.SaveToGephiGdf(temp + ".gdf",  parameters.SizeOfChange);
                                    if (parameters.showInputData)
                                        gng1.SaveInputDataToGephiGdf(temp + "InputData.gdf",  parameters.SizeOfChange);
                                }

							}
		 }
		  else
		  {
			  cout<<"The process without data:" << gng.Rank<<endl;
						 MPI_Comm_split(MPI_COMM_WORLD, 0, gng.Rank, &gng.comm);
	
		 }
	 }
		 MPI_Finalize();
	return 0;

}