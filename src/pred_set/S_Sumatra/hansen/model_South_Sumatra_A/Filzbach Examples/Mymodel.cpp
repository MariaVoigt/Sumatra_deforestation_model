//Predicts the potential magnitude and spatial pattern of Amazon deforestation
///Dynamic and spatially-explicit model of deforestation
//Differs from previous models in three ways: 
//(1) probabilistic rather than deterministic - quantify the uncertainty around our predictions
//(2) the overall deforestation rate emerges ‘bottom up’, as the sum of local-scale deforestation driven by local processes
//(3) deforestation is contagious, such that local deforestation rate increases through time if nearby locations are deforested.  
//REFERENCE: Rosa, IMD, Purves, D, Souza Jr, C and Ewers, RM "Predictive modelling of contagious deforestation in the Brazilian Amazon" 
//Acceptef for publication in PLOS One 2013
//Written by Isabel Rosa and Drew Purves 


// check upstairs RELEASE and x64 (thats faster)

//IMPORTANT: first need to install (freely available) Filzbach C++ library from http://research.microsoft.com/en-us/um/cambridge/groups/science/tools/filzbach/filzbach.htm

#include "preamble.h"
#ifdef Mymodel	//name of the model c++ file

//libraries to include
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "filzbach.h"

void pause() { PAUSE }

/************************************************************/
/* global variables                                         */
/************************************************************/

/* some things needed for book-keeping (see below for how each one is used) */
//confirm ascii file dimensions (number of columns and rows, change if appropriate)

// Maria - change dimensions of input data

//South Sumatra
#define GRIDNX 2458
#define GRIDNY 2157 // dimensions of grid we are using
#define XLLCORNER -2440978.004000000190
#define YLLCORNER  1090629.379999999888// coerner dimensions

#define CELLSIZE 180  //cellsize

//define all grids to be used here 
float forest[GRIDNX][GRIDNY];			// 0 = previous deforestation, 1= forest, -9999 = sea, non-forest vegetation. This is observed
float def[GRIDNX][GRIDNY];				// observed deforestation, -9999 = not deforested or not considered, 1 = forest converted to non-forest 
float slope[GRIDNX][GRIDNY];			//slope
float fire[GRIDNX][GRIDNY];			//fire
float access_hrs[GRIDNX][GRIDNY];			//access in hrs
float roads[GRIDNX][GRIDNY];			//roads distance in m
float rivers[GRIDNX][GRIDNY];			//rivers distance in m

float pp_1[GRIDNX][GRIDNY];			//pop pressure sigma 1

float subsistenceLH[GRIDNX][GRIDNY];			//subsistence_distance_non_forest	
float plantationLH[GRIDNX][GRIDNY];			//plantation_distance_non_forest
float non_agriLH[GRIDNX][GRIDNY];			//non_agri_distance_non_forest
float transmigrant[GRIDNX][GRIDNY];	    //LC_transmigrant_distance


float peat[GRIDNX][GRIDNY];			//peat, 2 levels, 0-1
float plantations[GRIDNX][GRIDNY];	    //plantations, 3 levels, 0-2
float mining[GRIDNX][GRIDNY];        //mining, 3 levels, 0-2
float landuse[GRIDNX][GRIDNY];        //landuse, 5 levels,0 - 4 
									 // 
float modelForest[GRIDNX][GRIDNY];		//This is going to be updated every year
float dataNeib[GRIDNX][GRIDNY];		//This is going to be updated every year
float bigDataNeib[GRIDNX][GRIDNY][12]; // final dimension is rho
float dataRoaded[GRIDNX][GRIDNY];		//1 or 0 for road, or not, in data
float modelRoaded[GRIDNX][GRIDNY];		//1 or 0 for road, or not, modelled (not used...)

float predprob[GRIDNX][GRIDNY];		//store predicted probabilities from model (used to get ROC and AUC value)

									//annual deforestation maps (1 - deforested, -9999 all the rest)
float modelDef[GRIDNX][GRIDNY];

int xval[GRIDNX * GRIDNY]; //vector to store sample rows ID
int yval[GRIDNX * GRIDNY]; //vector to store sample columns ID

						 //temporary copies of parameters to be included in the model
						 // Maria - change here for the number of variables you have
						 // IS THIS CORRECT THE 7 LEVELS OF GAMLU
double mybeta_0, mybeta_1, mybeta_2, mybeta_3, mybeta_4, mybeta_5, mybeta_6, mybeta_7,
mybeta_8, mybeta_9, mybeta_10, mybeta_11, gampeat[1], gamplant[2], gammine[2], gamlu[3], myrho;
int beta0_global = 1;

// Maria - change here for the total number of parameters
int tfree[15]; // for stepwise regression, flag for whether parameters are currently free (1), or fixed at 0.0 (0) 

			  //other global variables
double training_likelihood, current_likelihood;
int training = 2;
int resample = 1; // when this is one, likelihood selects a different set of cells at random
int ibest;

#define NBORMODE 1 // 0 means, flat nbhoorhood (all nbors count equally), 1 means proper nborhood (effect drops off with distance according to parameter rho)


/************************************************************/
/* function headers                                         */
/************************************************************/
void fake_data();
void read_data();
void setup_model();
void fit_model();
void final_output();
void readAsciiGrid(char filename[], float destinationArray[GRIDNX][GRIDNY]);
void writeAsciiGrid(char filename[], float sourceArray[GRIDNX][GRIDNY], int numx, int numy, double xcorner, double ycorner, double cellsize, double nodata);
void applyModel(int plookup);
double get_pdefor(int xx, int yy, int datamod);
double get_edist(int xx, int yy, int datamod);
double get_edist2(int xx, int yy, int datamod, int code);
double edist(int nx, int ny);
void get_params();
double get_neib(int xx, int yy, int datamod, double rho);
double get_neib2(int xx, int yy, int datamod, double rho);
void test_getneib(double rho);
void fill_dataNeib(double rho);
void fillBigDataNeib();
double nweight(int nx, int ny, double rho);


/************************************************************/
void readAsciiGrid(char filename[], float destinationArray[][GRIDNY])
{
	FILE* ifile;
	int ncols, nrows, xx, yy;
	char wst[100];
	double wsd;

	printf("\n readAsciiGrid: beginning %s \n", filename);
	//PAUSE;

	if (!(ifile = fopen(filename, "r"))) {
		printf("\n readAsciiGrid: cannot open file named %s \n", filename);
		system("pause");
	}
	else
	{
		printf("\n readAsciiGrid: file named %s opened OK \n", filename);


		// get ncols, nrows
		fscanf(ifile, "%s %d", &wst, &ncols);
		fscanf(ifile, "%s %d", &wst, &nrows);

		printf("\n reading ascii grid, ncols nrows is: %d %d \n", ncols, nrows);

		// ignore the rest
		fscanf(ifile, "%s %lf %s %lf %s %lf %s %lf", &wst, &wsd, &wst, &wsd, &wst, &wsd, &wst, &wsd);

		// now run over all cells and transfer to destination array as appropriate
		for (yy = nrows - 1; yy >= 0; yy--)
		{
			for (xx = 0; xx < ncols; xx++)
			{
				fscanf(ifile, "%lf", &wsd);
				if (xx < GRIDNX && yy < GRIDNY)
				{
					destinationArray[xx][yy] = wsd;
				}
			}
		}
	}
	printf("\n reading asciigrid (%s) finished \n", filename);
	fclose(ifile);

	return;

}

/************************************************************/
void writeAsciiGrid(char filename[], float sourceArray[GRIDNX][GRIDNY], int numx, int numy, double xcorner, double ycorner, double cellsize, double nodata)
{
	FILE* ofile;
	int ncols, nrows, xx, yy;
	char wst[100];
	double wsd;

	if (!(ofile = fopen(filename, "w"))) {
		printf("\n writeAsciiGrid: cannot open file named %s \n", filename);
		system("pause");
	}
	else
	{
		// write header info
		fprintf(ofile, "ncols %d\n", numx);
		fprintf(ofile, "nrows %d\n", numy);
		fprintf(ofile, "xllcorner %lf\n", xcorner);
		fprintf(ofile, "yllcorner %lf\n", ycorner);
		fprintf(ofile, "cellsize %lf\n", cellsize);
		fprintf(ofile, "NODATA_value %lf\n", nodata);

		// now write data
		for (yy = numy - 1; yy >= 0; yy--)
		{
			for (xx = 0; xx < numx; xx++)
			{
				fprintf(ofile, "%lf ", sourceArray[xx][yy]);
			}
			fprintf(ofile, "\n");
		}
	}

	fclose(ofile);
}

/************************************************************/
double get_pdefor(int xx, int yy, int datamod)
{

	// Maria - change here  for your model 

	//function to get deforestation probability based on several ascii grids

	double pdefor, k, mygamma_peat, mygamma_plant, mygamma_mine, mygamma_lu, mygamma_b0;
	//xx, yy are the coordinates of the grid-cell
	//datamod==0 means, get current pattern of forest, from our data grid
	//datamod==1 means, get current pattern of forest, from our model grid
	//in both cases, refer to several other grids (rivers, soil fert etc.)


	if (peat[xx][yy] == 0)
	{
		mygamma_peat = 0.0;
	}
	else
	{
		// otherwise, get relevant parameter
		mygamma_peat = gampeat[(int)peat[xx][yy] - 1];
	}

	if (plantations[xx][yy] == 0)
	{
		mygamma_plant = 0.0;
	}
	else
	{
		// otherwise, get relevant parameter
		mygamma_plant = gamplant[(int)plantations[xx][yy] - 1];
	}


	if (landuse[xx][yy] == 0)
	{
		mygamma_lu = 0.0;
	}
	else
	{
		// otherwise, get relevant parameter
		mygamma_lu = gamlu[(int)landuse[xx][yy] - 1];
	}


	if (mining[xx][yy] == 0)
	{
		mygamma_mine = 0.0;
	}
	else
	{
		// otherwise, get relevant parameter
		mygamma_mine = gammine[(int)mining[xx][yy] - 1];
	}

	//get correct gammine depending on level of protection
	//int p;



	//model
	k = mybeta_0 + mybeta_1 * get_neib(xx, yy, datamod, myrho) +
		mybeta_2 * slope[xx][yy] +
		mybeta_3 * fire[xx][yy] +
		mybeta_4 * access_hrs[xx][yy] +
		mybeta_5 * roads[xx][yy] +
		mybeta_6 * rivers[xx][yy] +
		mybeta_7 * pp_1[xx][yy] +
		mybeta_8 * subsistenceLH[xx][yy] +
		mybeta_9 * plantationLH[xx][yy] +
		mybeta_10 * non_agriLH[xx][yy] +
		mybeta_11 * transmigrant[xx][yy] +

		mygamma_peat * peat[xx][yy] +
		mygamma_plant * plantations[xx][yy] +
		mygamma_mine * mining[xx][yy] +
		mygamma_lu * landuse[xx][yy];

	pdefor = logistic(k);

	return pdefor;
}

/************************************************************/
void get_params() //small trick to make the model faster
{
	// Maria - change here for the number of variables 

	//continuous variables
	mybeta_0 = cv("beta_0");
	mybeta_1 = cv("beta_1");
	mybeta_2 = cv("beta_2");
	mybeta_3 = cv("beta_3");
	mybeta_4 = cv("beta_4");
	mybeta_5 = cv("beta_5");
	mybeta_6 = cv("beta_6");
	mybeta_7 = cv("beta_7");
	mybeta_8 = cv("beta_8");
	mybeta_9 = cv("beta_9");
	mybeta_10 = cv("beta_10");
	mybeta_11 = cv("beta_11");

	//categorical variable
	gampeat[0] = cv("gamma_peat", 0);

	gamplant[0] = cv("gamma_plant", 0);
	gamplant[1] = cv("gamma_plant", 1);

	gammine[0] = cv("gamma_mine", 0);
	gammine[1] = cv("gamma_mine", 1);

	gamlu[0] = cv("gamma_lu", 0);
	gamlu[1] = cv("gamma_lu", 1);
	gamlu[2] = cv("gamma_lu", 2);


	// scale for nborhood
	myrho = cv("rho");

	return;
}

/************************************************************/
double get_edist(int xx, int yy, int datamod)
{
	//function to get euclidean distance
	int nx, ny, fx, fy;
	double tedist, mindist;

	// initially set mindist to largest allowable value
	mindist = CELLSIZE * edist(10, 10);

	// return euc distance to nearest deforested cell
	// either for real pattern of forest at time 0 (datamod==0)
	// or from the model (datamod==1)
	for (nx = -10; nx <= 10; nx++)						//a moving window of 5 cells on each size
	{
		fx = xx + nx;
		if (fx >= 0 && fx < GRIDNX)
		{
			for (ny = -10; ny <= 10; ny++)
			{
				fy = yy + ny;
				if (fy >= 0 && fy < GRIDNY)
				{
					if (!(nx == 0 && ny == 0))			//To make sure it doesn't look at the corners of the grid
					{
						if ((datamod == 0 && forest[fx][fy] > (-900.0) && forest[fx][fy] < 1.0) || (datamod == 1 && modelForest[fx][fy] > (-900.0) && modelForest[fx][fy] < 1.0))
						{
							//edist = sqrt((double)nx*(double)nx+(double)ny*(double)ny);
							tedist = CELLSIZE * edist(nx, ny);
							if (tedist < mindist)
								mindist = tedist;
						}
					}
				}
			}
		}
	}

	return mindist;

}

/************************************************************/


/************************************************************/
double edist(int nx, int ny)
{
	double euc;
	static double euctable[20][20];
	static int first = 1;

	// first time the function runs, populate the table
	if (first == 1)
	{
		first = 0;
		for (int nx = 0; nx < 20; nx++)
			for (int ny = 0; ny < 20; ny++)
				euctable[nx][ny] = sqrt((double)nx * (double)nx + (double)ny * (double)ny);
	}

	if (abs(nx) < 20 && abs(ny) < 20)
		return euctable[abs(nx)][abs(ny)];
	else
	{
		return sqrt((double)nx * (double)nx + (double)ny * (double)ny);
	}

}

/************************************************************/
double nweight(int nx, int ny, double rho)
{
	double euc, weight;
	static double wtable[20][20];
	static int first = 1;
	static double prevrho = 0.0;

	// first time the function runs, populate the table
	if (first == 1 || (rho - prevrho) > 0.10 || (prevrho - rho) > 0.10)
	{
		printf("\n recalc nweight array \n");
		first = 0;
		for (int nx = 0; nx < 20; nx++)
			for (int ny = 0; ny < 20; ny++)
				wtable[nx][ny] = exp((-1.0) * pow((CELLSIZE * edist(nx, ny)) / rho, 2.0));
	}

	if (abs(nx) < 20 && abs(ny) < 20)
		weight = wtable[abs(nx)][abs(ny)];
	else
	{
		weight = exp((-1.0) * pow((CELLSIZE * edist(nx, ny)) / rho, 2.0));
	}

	prevrho = rho;

	return weight;

}

/************************************************************/
void fill_dataNeib(double rho)
{
	static int first = 1;
	static double prevrho;

	if (first > 0 || (prevrho - rho) > 0.10 || (rho - prevrho) > 0.10)
	{
		first = 0;
		printf("\n filling dataneib array again \n");
		for (int xx = 0; xx < GRIDNX; xx++)
		{
			for (int yy = 0; yy < GRIDNY; yy++)
			{
				if (forest[xx][yy] > (-900.0))
					dataNeib[xx][yy] = get_neib2(xx, yy, 0, rho);
				else
					dataNeib[xx][yy] = 0.0;
			}
		}
		printf("\n done \n");

	}

	prevrho = rho;

	return;
}

/************************************************************/
double get_neib2(int xx, int yy, int datamod, double rho)
{
	int nx, ny, fx, fy;
	double sneib1 = 0.0, sneib2 = 0.0;

	for (nx = -10; nx <= 10; nx++)						//a moving window of 5 cells on each size
	{
		fx = xx + nx;
		if (fx >= 0 && fx < GRIDNX)						//To make sure is inside the GRID x limit
		{
			for (ny = -10; ny <= 10; ny++)
			{
				fy = yy + ny;
				if (fy >= 0 && fy < GRIDNY)				//To make sure is inside the GRID y limit
				{
					if (!(nx == 0 && ny == 0))			//To make sure it doesn't look at the corners of the grid
					{
						// is this neighbouring cell suitable for counting?
						if ((datamod == 0 && forest[fx][fy] > (-900.0)) || (datamod == 1 && modelForest[fx][fy] > (-900.0)))
						{
							if (NBORMODE == 0)
								sneib2 += 1.0;      //sum2 wherever forest is 0 or 1
							else
								sneib2 += nweight(nx, ny, rho);
						}
						// is this suitable cell, deforested?
						if ((datamod == 0 && forest[fx][fy] > (-900.0) && forest[fx][fy] < 1.0) || (datamod == 1 && modelForest[fx][fy] > (-900.0) && modelForest[fx][fy] < 1.0))
						{
							if (NBORMODE == 0)
								sneib1 += 1.0;		//sum 1 whenever forest == 0.0
							else
								sneib1 += nweight(nx, ny, rho);
						}
					}
				}
			}
		}
	}

	if (sneib2 > 0.0)
		sneib1 /= sneib2; //fraction of neighbours that are deforested
	else
		sneib1 = 0.0;

	//if(random(0.0,1.0)>0.99)
	//printf(" :) %lf",sneib1);

	return sneib1;

}

/************************************************************/
double get_neib(int xx, int yy, int datamod, double rho)
{
	static int first = 1;
	double v1, v2, l1, l2;
	//return the sum of deforestation events within a 5x5 window of the deforested cell
	//either for real pattern of forest at time 0 (datamod==0)
	//or from the model (datamod==1)

	if (datamod == 0 && first == 1)
	{
		first = 0;
		fillBigDataNeib();
	}

	if (datamod == 0)
	{
		//fill_dataNeib(rho); // fill array if neccessary
		//return dataNeib[xx][yy]; // from data
		int kk = (rho - CELLSIZE) / (CELLSIZE * 2.0);
		if (kk < 0)
			kk = 0;
		if (kk > 10)
			kk = 10;
		v1 = bigDataNeib[xx][yy][kk];
		l1 = CELLSIZE + (CELLSIZE * 2.0) * ((double)kk);
		kk++;
		if (kk < 0)
			kk = 0;
		if (kk > 10)
			kk = 10;
		v2 = bigDataNeib[xx][yy][kk];
		l2 = CELLSIZE + (CELLSIZE * 2.0) * ((double)kk);

		return v1 + (v2 - v1) * ((rho - l1) / (l2 - l1));
	}
	else
	{
		return get_neib2(xx, yy, 1, rho); // from model grid
	}

}

/************************************************************/
void fillBigDataNeib()
{
	int kk;
	double rho;

	for (kk = 0; kk <= 10; kk++)
	{
		printf("\n filling big data neib, going to next rho value... \n");
		// set value of rho
		rho = CELLSIZE + (CELLSIZE * 2.0) * ((double)kk);

		// run over grid and fill in big data table
		for (int xx = 0; xx < GRIDNX; xx++)
		{
			for (int yy = 0; yy < GRIDNY; yy++)
			{
				if (forest[xx][yy] > (-900.0))
					bigDataNeib[xx][yy][kk] = get_neib2(xx, yy, 0, rho);
				else
					bigDataNeib[xx][yy][kk] = 0.0;
			}
		}
	}

	return;
}

/************************************************************/
void applyModel(int plookup)
{
	//function for deforestation simulations
	//auxiliary variables
	static int first = 1;
	static int pp = 0; //counter for output
	char fname[200];

	//get copy of parameters
	if (plookup == 1)
		get_params();

	//create table for output of various rates - holds forest data through time
	if (first == 1)
	{
		first = 0;
		table_create("forestdata");
		table_addcolumn("forestdata", "Iteration");
		table_addcolumn("forestdata", "Year");
		table_addcolumn("forestdata", "ForestLeft");
		table_addcolumn("forestdata", "TotalDef");
		table_addcolumn("forestdata", "RateDef");

		//create table to output parameters used in each iteration - uncertainty in parameters
		table_create("parameters");
		table_addcolumn("parameters", "Iteration");
		table_addcolumn("parameters", "beta_0");
		table_addcolumn("parameters", "beta_1");
		table_addcolumn("parameters", "beta_2");
		table_addcolumn("parameters", "beta_3");
		table_addcolumn("parameters", "beta_4");
		table_addcolumn("parameters", "beta_5");
		table_addcolumn("parameters", "beta_6");
		table_addcolumn("parameters", "beta_7");
		table_addcolumn("parameters", "beta_8");
		table_addcolumn("parameters", "beta_9");
		table_addcolumn("parameters", "beta_10");
		table_addcolumn("parameters", "beta_11");

		table_addcolumn("parameters", "gampeat_0");

		table_addcolumn("parameters", "gammine_0");
		table_addcolumn("parameters", "gammine_1");

		table_addcolumn("parameters", "gamplant_0");
		table_addcolumn("parameters", "gamplant_1");

		table_addcolumn("parameters", "gamlu_0");
		table_addcolumn("parameters", "gamlu_1");
		table_addcolumn("parameters", "gamlu_2");
	}

	//other auxiliary variables
	int j, n, xx, yy;
	double ypred, pdefor;


	//define the number of years we want the model to run and the number of iterations (error propagation loop)
	const int numsteps = 8, iterations = 100;


	//Iterations loop starts here (number of times the model will run)
	for (j = 0; j < iterations; j++)
	{
		printf("\n deforest model, started iteration %d", j);

		//draw a new set of input parameter values - uncertainty in the parameters
		params_draw_random_vector();

		//initialize model first - forest cover
		for (xx = 0; xx < GRIDNX; xx++)
			for (yy = 0; yy < GRIDNY; yy++)
				modelForest[xx][yy] = forest[xx][yy];

		//at this point, model forest is same as observed forest at time 0
		//but now, we are going to start altering the model forest grid according to the model

		//this is just to set our modelDef grig to be zero at the beginnig of a new iteration loop as well as the predicted probability from the model
		//modelDef will store 1 if the pixel was deforested in a given year or -9999.0 for everything else
		//predprob will store the predicted probability for each pixel where forest = 1 in a given year or -9999.0 for everything else
		for (xx = 0; xx < GRIDNX; xx++)
			for (yy = 0; yy < GRIDNY; yy++)
				if (modelForest[xx][yy] < (-99.0))
				{
					modelDef[xx][yy] = -9999.0;
					predprob[xx][yy] = -9999.0;
				}


		//Time steps loop starts here (number of years the model will run in each iteration)
		for (n = 0; n < numsteps; n++)
		{
			table_writevalue("forestdata", "Iteration", pp, j);
			table_writevalue("forestdata", "Year", pp, n);

			//set everything to 0 - updated every year
			double totfor = 0.0, totdef = 0.0, ratedef = 0.0;
			for (xx = 0; xx < GRIDNX; xx++)
			{
				for (yy = 0; yy < GRIDNY; yy++)
				{
					//to get the annual deforestation, need to reset this data to zero every year
					modelDef[xx][yy] = 0.0;

					if (modelForest[xx][yy] > 0.0)
					{
						//sum forest data to calculate remaining forest	
						totfor = totfor + modelForest[xx][yy];
					}

				}
			}
			//store forest cover
			table_writevalue("forestdata", "ForestLeft", pp, totfor);

			//get deforestation rate
			int forn, forn1;
			if (n > 0)
			{
				forn = table_getvalue("forestdata", "ForestLeft", pp);				//gets forest in time step n
				forn1 = table_getvalue("forestdata", "ForestLeft", pp - 1);				//gets forest in time step n-1
				totdef = forn1 - forn;												//total deforestation = forest in n-1 - forest in n
				table_writevalue("forestdata", "TotalDef", pp, totdef);
				ratedef = ((double)totdef / (double)totfor) * 100;						//rate of defor = (tot def n / tot for n )*100
				table_writevalue("forestdata", "RateDef", pp, ratedef);
			}

			//xx, yy is our 'focal' cell that we may or may not deforest
			for (xx = 0; xx < GRIDNX; xx++)
			{
				for (yy = 0; yy < GRIDNY; yy++)
				{
					if (modelForest[xx][yy] > 0.0)			// only need to look at forested cells
					{
						//calls pdefor function to get probability of deforestation
						pdefor = get_pdefor(xx, yy, 1);

						//store the probability values in a matrix to export as ascii and build ROC curve and get AUC value
						predprob[xx][yy] = pdefor;

						//Pixel to be deforested - Flip weighted coin selection
						if (random(0.0, 1.0) < pdefor)
						{
							modelForest[xx][yy] = 0.0;
							modelDef[xx][yy] = 1.0;
						}

					}

				} //loop over yy
			} //loop over xx


			  //Hardcode this to S
			  //EDIT

			  //Exporting predicted probabilites to get AUC -- MARIA change xll yll corners
			printf("\n Exporting predicted probabilities \n");
			sprintf(fname, "E:/Sumatra_model_August22/tmf/S_Sumatra/model_S_Sumatra_A/predprob_yr%d_i%d.asc", n, j);
			writeAsciiGrid(fname, predprob, GRIDNX, GRIDNY, XLLCORNER, YLLCORNER, CELLSIZE, -9999.0);

			printf("\n deforest model, year %d", n);

			//move output counter on 1
			pp++;


			//export the annual forest cover map as an ASCII file
			printf("\n Exporting new forest cover map \n");
			sprintf(fname, "E:/Sumatra_model_August22/tmf/S_Sumatra/model_S_Sumatra_A/predfor_i%d_%dyr.asc", j, n);
			writeAsciiGrid(fname, modelForest, GRIDNX, GRIDNY, XLLCORNER, YLLCORNER, CELLSIZE, -9999.0);

			//export the annual deforestation map as an ASCII file
			printf("\n Exporting deforestation map for time %d \n", n);
			sprintf(fname, "E:/Sumatra_model_August22/tmf/S_Sumatra/model_S_Sumatra_A/preddef_i%d_%dyr.asc", j, n);
			writeAsciiGrid(fname, modelDef, GRIDNX, GRIDNY, XLLCORNER, YLLCORNER, CELLSIZE, -9999.0);



		} //End time steps loop 


		  //write to table the values of the parameters used in this iteration
		table_writevalue("parameters", "beta_0", j, cv("beta_0"));
		table_writevalue("parameters", "beta_1", j, cv("beta_1"));
		table_writevalue("parameters", "beta_2", j, cv("beta_2"));
		table_writevalue("parameters", "beta_3", j, cv("beta_3"));
		table_writevalue("parameters", "beta_4", j, cv("beta_4"));
		table_writevalue("parameters", "beta_5", j, cv("beta_5"));
		table_writevalue("parameters", "beta_6", j, cv("beta_6"));
		table_writevalue("parameters", "beta_7", j, cv("beta_7"));
		table_writevalue("parameters", "beta_8", j, cv("beta_8"));
		table_writevalue("parameters", "beta_9", j, cv("beta_9"));
		table_writevalue("parameters", "beta_10", j, cv("beta_10"));
		table_writevalue("parameters", "beta_11", j, cv("beta_11"));

		table_writevalue("parameters", "gampeat_0", j, cv("gamma_peat", 0));

		table_writevalue("parameters", "gamplant_0", j, cv("gamma_plant", 0));
		table_writevalue("parameters", "gamplant_1", j, cv("gamma_plant", 1));

		table_writevalue("parameters", "gammine_0", j, cv("gamma_mine", 0));
		table_writevalue("parameters", "gammine_1", j, cv("gamma_mine", 1));


		table_writevalue("parameters", "gamlu_0", j, cv("gamma_lu", 0));
		table_writevalue("parameters", "gamlu_1", j, cv("gamma_lu", 1));
		table_writevalue("parameters", "gamlu_2", j, cv("gamma_lu", 2));

		//Export forestdata table
		sprintf(fname, "./workspace/Model_35yrs_temp.txt");
		table_output("forestdata", fname);

	} //End error propagation loop 

	  //EDIT
	  //Export forestdata table
	sprintf(fname, "./workspace/Model_35yrs.txt");
	table_output("forestdata", fname);

	//Export parameters table 
	table_output("parameters", "./workspace/Parameters_Values.txt");


	return;
} //applymodel() function ends

  /************************************************************/
void read_data()
{

	// Maria - reads in the data 

	//read various ascii grids with input data
	//EDIT
	readAsciiGrid("./workspace/forest_2017_21_180m_repro_res_tmf_S_Sumatra.ascii", forest);		//
	readAsciiGrid("./workspace/deforestation_2017_21_180m_repro_res_tmf_S_Sumatra.ascii", def);		//
	readAsciiGrid("./workspace/slope_180m_repro_res_tmf_S_Sumatra.ascii", slope);		//
	readAsciiGrid("./workspace/fire_yearly_average_180m_repro_res_tmf_S_Sumatra.ascii", fire);		//
	readAsciiGrid("./workspace/IDN_TTCSM_hrs_180m_repro_res_tmf_S_Sumatra.ascii", access_hrs);		//
	readAsciiGrid("./workspace/road_distance_180m_repro_res_tmf_S_Sumatra.ascii", roads);		//
	readAsciiGrid("./workspace/river_distance_180m_repro_res_tmf_S_Sumatra.ascii", rivers);		//
	readAsciiGrid("./workspace/pressurelog10_sigma1_180m_repro_res_tmf_S_Sumatra.ascii", pp_1);		//
	readAsciiGrid("./workspace/subsistence_LH_distance_180m_repro_res_tmf_S_Sumatra.ascii", subsistenceLH);		//
	readAsciiGrid("./workspace/plantation_LH_distance_180m_repro_res_tmf_S_Sumatra.ascii", plantationLH);		//
	readAsciiGrid("./workspace/non_agri_LH_distance_180m_repro_res_tmf_S_Sumatra.ascii", non_agriLH);		//
	readAsciiGrid("./workspace/transmigrant_distance_180m_repro_res_tmf_S_Sumatra.ascii", transmigrant);		//

	readAsciiGrid("./workspace/peat_180m_repro_res_tmf_S_Sumatra.ascii", peat);		//peat
	readAsciiGrid("./workspace/plantations_180m_repro_res_tmf_S_Sumatra.ascii", plantations);		//
	readAsciiGrid("./workspace/mining_180m_repro_res_tmf_S_Sumatra.ascii", mining);		//mining
	readAsciiGrid("./workspace/lu_180m_repro_res_tmf_S_Sumatra.ascii", landuse);		//landuse
	return;

}


/************************************************************/
void setup_model()
{

	// Maria - add parameter intervals here

	/* each line defines one new parameter for the model */
	if (beta0_global == 1)
		parameter_create("beta_0", -6.0, 6.0, 0.0, 0, 0, 1);		//intercept
	else
		parameter_create_vector("gamma_b0", -6.0, 6.0, 0.0, 0, 0, 1, 9);	//vector to store different b0 for each state (not used...)

																			// continuous variables
	parameter_create("beta_1", -6.0, 8.0, 0.0, 0, 0, 1);		//proportion of deforested neighbours (updated every year, dynamic variable)
	parameter_create("beta_2", -1.0000, 0.0001, 0.0, 0, 0, 1);		//slope
	parameter_create("beta_3", -2.0, 4.0, 0.0, 0, 0, 1);		//fire
	parameter_create("beta_4", -2.0, 2.0, 0.0, 0, 0, 1);		//access_hrs
	parameter_create("beta_5", -1.0, 1.0, 0.0, 0, 0, 1);		//roads
	parameter_create("beta_6", -1.0, 1.0, 0.0, 0, 0, 1);		//rivers
	parameter_create("beta_7", -2.0, 2.0, 0.0, 0, 0, 1);		//population pressure sigma 1
	parameter_create("beta_8", -1.0, 1.0, 0.0, 0, 0, 1);		//distance to subsistence
	parameter_create("beta_9", -1.0, 1.0, 0.0, 0, 0, 1);		//distance to plantation
	parameter_create("beta_10", -1.0, 1.0, 0.0, 0, 0, 1);		//distance to non_agri
	parameter_create("beta_11", -2.0, 2.0, 0.0, 0, 0, 1);		//transmigrant

	parameter_create_vector("gamma_peat", -2.0, 2.0, 0.0, 0, 0, 1, 1);// // last digit is 1+ the last element in get_params, to allow a parameter value for each class (excluding zero!)

	parameter_create_vector("gamma_plant", -2.0, 2.0, 0.0, 0, 0, 1, 2);// // last digit is 1+ the last element in get_params, to allow a parameter value for each class (excluding zero!)

	parameter_create_vector("gamma_mine", -2.0, 2.0, 0.0, 0, 0, 1, 2);// // last digit is 1+ the last element in get_params, to allow a parameter value for each class (excluding zero!)

	parameter_create_vector("gamma_lu", -2.0, 6.0, 0.0, 0, 0, 1, 3);// EDIT
	parameter_create("rho", CELLSIZE, (CELLSIZE * 20.0), (CELLSIZE * 2.0), 0, 0, 1); // distance, in km, at which the effect of a deforested nbor drops to 1/e of the max value
															   // For stepwise regression
	if (tfree[0] == 0)
		parameter_fix("beta_1", 0.0); // name, and value to fix at
	if (tfree[1] == 0)
		parameter_fix("beta_2", 0.0); // name, and value to fix at
	if (tfree[2] == 0)
		parameter_fix("beta_3", 0.0); // name, and value to fix at
	if (tfree[3] == 0)
		parameter_fix("beta_4", 0.0); // name, and value to fix at
	if (tfree[4] == 0)
		parameter_fix("beta_5", 0.0); // name, and value to fix at
	if (tfree[5] == 0)
		parameter_fix("beta_6", 0.0); // name, and value to fix at
	if (tfree[6] == 0)
		parameter_fix("beta_7", 0.0); // name, and value to fix at
	if (tfree[7] == 0)
		parameter_fix("beta_8", 0.0); // name, and value to fix at
	if (tfree[8] == 0)
		parameter_fix("beta_9", 0.0); // name, and value to fix at
	if (tfree[9] == 0)
		parameter_fix("beta_10", 0.0); // name, and value to fix at
	if (tfree[10] == 0)
		parameter_fix("beta_11", 0.0); // name, and value to fix at

	if (tfree[11] == 0)
		parameter_fix("gamma_peat", 0.0); // name, and value to fix at
	if (tfree[12] == 0)
		parameter_fix("gamma_plant", 0.0); // name, and value to fix at
	if (tfree[13] == 0)
		parameter_fix("gamma_mine", 0.0); // name, and value to fix at
	if (tfree[14] == 0)
		parameter_fix("gamma_lu", 0.0); // name, and value to fix at
	parameter_showall();

	return;
}


/************************************************************/
void likelihood()
{
	/* writes the log-likelihood given current parameter values */

	double prob1, ypred, pdefor;
	double ltot = 0.0, ssize = 0.0;

	// make sure we get temp copies of params
	get_params();

	/* set sum over log-likelihood to zero */
	set_metr_ltotnew(0.0);
	set_metr_number_ok(0.0);

	/* loop over the data (e.g. taking a random sample from the many grid-cells) */
	const int samplesize = 10000;
	int z, xx, yy, numrow, numcol;

	int hits = 0;
	static int osamplesize;

	if (training == 1 && (resample == 1 || random(0.0, 1.0) > 0.90))
	{
		// if in training phase, take random sample from even y values, any x values


		while (hits < samplesize)
		{
			xx = random_integer(0, GRIDNX - 1);
			yy = random_integer(0, GRIDNY - 1);
			if (yy % 2 == 0 && forest[xx][yy] > 0.0)
			{
				xval[hits] = xx;
				yval[hits] = yy;
				hits++;
			}
		}
		osamplesize = samplesize;

		resample = 0;
		params_set_to_old();
		get_params();
		force_accept();
	}

	if (training == 0 || training == 2)
	{
		// next time we call training mode, we're going to need to resample
		//resample=1;
		int bb;
		if (training == 0)
			bb = 1;
		else
			bb = 0;
		// if in testing mode, use ALL suitable cells from odd y values
		for (xx = 0; xx < GRIDNX; xx++)
		{
			for (yy = 0; yy < GRIDNY; yy++)
			{
				if ((yy + bb) % 2 == 0 && forest[xx][yy] > 0.0)
				{
					xval[hits] = xx;
					yval[hits] = yy;
					hits++;
				}
			}
		}
		// set osamplesize
		osamplesize = hits;
	}

	for (z = 0; z < osamplesize; z++)
	{
		numrow = xval[z];
		numcol = yval[z];

		/* for each sample, calculate prob defor and compare to obs defor */
		if (forest[numrow][numcol] > 0.0) // forest is the initial state, at 2007: 1 means forest, 0 means deforestation, -9999 means everything else
		{
			pdefor = get_pdefor(numrow, numcol, 0); // 0 means, do this for the observed data, rather than the model grid
		}

		/* compare to real data */
		if (def[numrow][numcol] > 0.0) // def=1 means this cell was deforested between 2001 and 2002
			prob1 = pdefor;
		else
			prob1 = 1.0 - pdefor;

		if (prob1 < 0.000000010)
			prob1 = 0.000000010;

		if (prob1 > 0.999999999)
			prob1 = 0.999999999;

		ltot += log(prob1);
		ssize += 1.0;
	}

	set_metr_ltotnew(ltot);
	set_metr_number_ok(ssize);

	// also set our global 'current likelihood' metric
	current_likelihood = ltot;

	return;
}


/* ************************************************* */
/* The control function that calls everything else.  */
/* ************************************************* */
int main()
{
	atexit(pause);

	// set the likelihood function pointer
	pfn_likelihood = &likelihood;

	/***********************************/
	/* Stepwise regression starts here */
	/***********************************/

	int turn; // one turn, is one run over all params trying each one as free rather than fixed
	int tp; // parameter we are currently testing as free, rather than fixed
	int p;
	int bfree; // parameter which, in this turn, is giving us the best increase in fit when we make it free
	double dfit; // the difference in fit when we make tp free
	double bfit = 0.0; // baseline fit
	double bestdfit; // the best dfit we have found in this turn

					 // Maria - change to total number of parameters
	int pfree[15]; // this is 1 if the parameter is set to be ALWAYS free
	int npara = 15; // number of params we are considering as fixed or not

	static int rr = 0; //counter for output model likelihood table

	char aname[200];


	//if 1 - runs forward stepwise regression
#if 1
	//Create "likelihood" table
	//add 8 columns to store if the parameters if fixed (0) or free (1) in each run of the stepwise regression
	table_create("model_likelihood");
	table_addcolumn("model_likelihood", "Run");				//Model run
	table_addcolumn("model_likelihood", "beta_0");
	table_addcolumn("model_likelihood", "beta_1");
	table_addcolumn("model_likelihood", "beta_2");
	table_addcolumn("model_likelihood", "beta_3");
	table_addcolumn("model_likelihood", "beta_4");
	table_addcolumn("model_likelihood", "beta_5");
	table_addcolumn("model_likelihood", "beta_6");
	table_addcolumn("model_likelihood", "beta_7");
	table_addcolumn("model_likelihood", "beta_8");
	table_addcolumn("model_likelihood", "beta_9");
	table_addcolumn("model_likelihood", "beta_10");
	table_addcolumn("model_likelihood", "beta_11");

	table_addcolumn("model_likelihood", "gamma_peat");
	table_addcolumn("model_likelihood", "gamma_plant");
	table_addcolumn("model_likelihood", "gamma_mine");
	table_addcolumn("model_likelihood", "gamma_lu");
	table_addcolumn("model_likelihood", "traininglikelihood");
	table_addcolumn("model_likelihood", "testlikelihood");

	//add 8 columns to store the posterior mean
	table_addcolumn("model_likelihood", "beta_0pm");
	table_addcolumn("model_likelihood", "beta_1pm");
	table_addcolumn("model_likelihood", "beta_2pm");
	table_addcolumn("model_likelihood", "beta_3pm");
	table_addcolumn("model_likelihood", "beta_4pm");
	table_addcolumn("model_likelihood", "beta_5pm");
	table_addcolumn("model_likelihood", "beta_6pm");
	table_addcolumn("model_likelihood", "beta_7pm");
	table_addcolumn("model_likelihood", "beta_8pm");
	table_addcolumn("model_likelihood", "beta_9pm");
	table_addcolumn("model_likelihood", "beta_10pm");
	table_addcolumn("model_likelihood", "beta_11pm");

	table_addcolumn("model_likelihood", "gamma_0peat");

	table_addcolumn("model_likelihood", "gamma_0plant");
	table_addcolumn("model_likelihood", "gamma_1plant");

	table_addcolumn("model_likelihood", "gamma_0mine");
	table_addcolumn("model_likelihood", "gamma_1mine");

	table_addcolumn("model_likelihood", "gamma_0lu");
	table_addcolumn("model_likelihood", "gamma_1lu");
	table_addcolumn("model_likelihood", "gamma_2lu");

	//Read input data from files:
	read_data();

	// now do our turns. -1 turn is just for set up
	for (turn = 0; turn <= npara; turn++)
	{
		if (turn == 0)
		{
			// initialize by setting all permanent frees to zero
			for (p = 0; p < npara; p++)
				pfree[p] = 0;
		}

		// loop over params and try freeing one at a time
		int hits = 0;
		for (tp = -1; tp < npara; tp++)
		{
			// set any non-permanent params to fixed
			for (p = 0; p < npara; p++)
			{
				if (pfree[p] == 0)
					tfree[p] = 0;
				else
					tfree[p] = 1;
			}


			// now free up our focal param and fit the resultant model, if not fixed already, unless tp negative
			if (tp >= 0 && pfree[tp] == 1)
			{
				// this param already permanently fixed at free -- do nothing
			}
			else
			{
				// if tp negative, don't touch anything
				if (tp >= 0)
					tfree[tp] = 1;

				// for info, print current pattern of tfree,pfree to screen
				printf("\n");
				for (p = 0; p < npara; p++)
					printf(" %d:%d ", pfree[p], tfree[p]);


				//start filzbach and setup model
				initialize_filzbach();
				beta0_global = 1;
				setup_model();
				sprintf(aname, "´S_Sumatra_1_%d", rr);
				name_analysis(aname);					//This gives a name to the analysis. All the outputs will have this extension.
				set_chains(1);

				// tell likelihood to use training data
				training = 2;
				// run mcmc
				runmcmc(10000, 10000, 1000, 1000);

				training = 2;
				resample = 0;
				likelihood();
				double training_likelihood = current_likelihood;

				// test model
				training = 0;
				likelihood();


				// out row of info to table here
				// which params free, or fixed, and training likelihood, and test likelihood
				table_writevalue("model_likelihood", "beta_0", rr, 1);
				table_writevalue("model_likelihood", "beta_1", rr, tfree[0]);
				table_writevalue("model_likelihood", "beta_2", rr, tfree[1]);
				table_writevalue("model_likelihood", "beta_3", rr, tfree[2]);
				table_writevalue("model_likelihood", "beta_4", rr, tfree[3]);
				table_writevalue("model_likelihood", "beta_5", rr, tfree[4]);
				table_writevalue("model_likelihood", "beta_6", rr, tfree[5]);
				table_writevalue("model_likelihood", "beta_7", rr, tfree[6]);
				table_writevalue("model_likelihood", "beta_8", rr, tfree[7]);
				table_writevalue("model_likelihood", "beta_9", rr, tfree[8]);
				table_writevalue("model_likelihood", "beta_10", rr, tfree[9]);
				table_writevalue("model_likelihood", "beta_11", rr, tfree[10]);

				table_writevalue("model_likelihood", "gamma_peat", rr, tfree[11]);
				table_writevalue("model_likelihood", "gamma_plant", rr, tfree[12]);
				table_writevalue("model_likelihood", "gamma_mine", rr, tfree[13]);
				table_writevalue("model_likelihood", "gamma_lu", rr, tfree[14]);



				table_writevalue("model_likelihood", "Run", rr, rr);
				table_writevalue("model_likelihood", "traininglikelihood", rr, training_likelihood);
				table_writevalue("model_likelihood", "testlikelihood", rr, current_likelihood);


				//get posterior mean for each parameter and store it in the model likelihood table
				table_writevalue("model_likelihood", "beta_0pm", rr, cv("beta_0"));
				table_writevalue("model_likelihood", "beta_1pm", rr, cv("beta_1"));
				table_writevalue("model_likelihood", "beta_2pm", rr, cv("beta_2"));
				table_writevalue("model_likelihood", "beta_3pm", rr, cv("beta_3"));
				table_writevalue("model_likelihood", "beta_4pm", rr, cv("beta_4"));
				table_writevalue("model_likelihood", "beta_5pm", rr, cv("beta_5"));
				table_writevalue("model_likelihood", "beta_6pm", rr, cv("beta_6"));
				table_writevalue("model_likelihood", "beta_7pm", rr, cv("beta_7"));
				table_writevalue("model_likelihood", "beta_8pm", rr, cv("beta_8"));
				table_writevalue("model_likelihood", "beta_9pm", rr, cv("beta_9"));
				table_writevalue("model_likelihood", "beta_10pm", rr, cv("beta_10"));
				table_writevalue("model_likelihood", "beta_11pm", rr, cv("beta_11"));

				table_writevalue("model_likelihood", "gamma_0peat", rr, cv("gamma_peat", 0));

				table_writevalue("model_likelihood", "gamma_0plant", rr, cv("gamma_plant", 0));
				table_writevalue("model_likelihood", "gamma_1plant", rr, cv("gamma_plant", 1));

				table_writevalue("model_likelihood", "gamma_0mine", rr, cv("gamma_mine", 0));
				table_writevalue("model_likelihood", "gamma_1mine", rr, cv("gamma_mine", 1));

				table_writevalue("model_likelihood", "gamma_0lu", rr, cv("gamma_lu", 0));
				table_writevalue("model_likelihood", "gamma_1lu", rr, cv("gamma_lu", 1));
				table_writevalue("model_likelihood", "gamma_2lu", rr, cv("gamma_lu", 2));

				//output likelihood table 
				table_output("model_likelihood", "./workspace/TempLhoodDdef_1.txt");

				//move counter
				rr++;

				if (tp < 0) {
					bfit = current_likelihood;
					dfit = 0.0;
					bestdfit = 0.0;
					//printf("\n doing baseline for this turn \n");
				}
				else
				{
					dfit = current_likelihood - bfit;
					if (hits == 0 || dfit > bestdfit)
					{
						printf("\n new bestdfit found: %lf for param %d \n", dfit, tp);
						bestdfit = dfit;
						bfree = tp;
						hits++;
					}
				}
			}
		} // end of tp loop

		  // now we have the info to set one param at permanently free
		printf("\n permanently freeing param: %d ", bfree);
		pfree[bfree] = 1;
	}

	//EDIT

	//output likelihood table 
	table_output("model_likelihood", "./workspace/Models_Stepwise_1.txt");

#endif

	//if 1 runs simulations using defined model
#if 0
	//read data
	read_data();

	//start filzbach and setup best model
	initialize_filzbach();
	beta0_global = 1;	//global (1) or "regional"(0)? regional option not being used
						//if tfree[]=1 parameter is included in the model
	tfree[0] = 1; //Deforestation
	tfree[1] = 1; //slope
	tfree[2] = 1; //fire
	tfree[3] = 0; //access_hrs
	tfree[3] = 0; //roads
	tfree[3] = 0; //rivers
	tfree[4] = 0; //pp1
	tfree[5] = 1; //dist sustainability
	tfree[6] = 1; //dist plantation
	tfree[7] = 1; //dist non-agri
	tfree[8] = 1; //transmigrant

	tfree[9] = 1; //peat
	tfree[10] = 1; //plantations
	tfree[11] = 1; //mine
	tfree[12] = 1; //landuse

	setup_model();
	name_analysis("S_Sumatra_model");
	set_chains(1);

	//tell likelihood to use training data
	training = 2;
	//run mcmc
	// run here with 10000, 10000, outputting every thousand iterations and you can put 5000- > then a little bit faster
	runmcmc(20000, 20000, 5000, 5000);
	training = 2;
	resample = 0;
	likelihood();
	double training_likelihood = current_likelihood;
	// test model
	training = 0;
	likelihood();

	//run simulations
	applyModel(1);	//look-up values from Filzbach

#endif

}

#endif

