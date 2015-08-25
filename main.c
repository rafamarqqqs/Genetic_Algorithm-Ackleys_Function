#include <stdio.h>
#include <limits.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#define CROSSOVER_TYPE ("Average beetween parent")
#define PARENTS_SELECTION ("Roullete method")
#define INITIALIZATION_TYPE ("Aleatory")
#define STOP_CONDITION (" fitness calls OR optimum achived")

#define DIM 8
#define C1 20.0
#define C2 0.20
#define C3 2.0*M_PI
#define CALLS_LIMIT 50000
#define POPULATION_SIZE 200
#define CROSSOVER_RATE 100		//max 100
#define MUTATION_SIZE 0.01
#define MUTATION_RATE 10		//max 100
#define RAND_DOUBLE(interval) ((double) rand() * interval / ((double) RAND_MAX))

//struct used to store data from a gene
 typedef struct chrom {
	double genes[DIM];
	double fitness;
} Chrom;

//struct used to create the roullete
typedef struct roulleteStruct {
	double value;
	Chrom chrom;
} RoulleteStruct;

//function to specify the type of parent selection
typedef Chrom* (*ParentSelectionFunction)(Chrom *, Chrom *);

//function to specify the type of crossover
typedef void (*CrossoverFunction)(Chrom *, Chrom *);

//function to specify the type of mutation 
typedef void (*MutationFunction)(Chrom *);

void printGeneration(Chrom *chroms){
	int i, j;

	for(i = 0; i < POPULATION_SIZE; i++){
		for(j = 0; j < DIM; j++)
			printf("%lf ", chroms[i].genes[j]);
		printf("\n");
	}
}

//calculates the ackley's function
double ackleysFunction(double *x, double c1, double c2, double c3){
	int i;
	double s1 = 0;
	double s2 = 0;
	
	for(i = 0; i < DIM; i++){
		s1 += x[i]*x[i];
		s2 += cos(c3*x[i]);
	}

	return (-1.0)*c1 * exp((-1.0) * c2 * sqrt((1.0/DIM) * s1)) - exp((1.0/DIM) * s2) + M_E + 20;
}

//roullete type classification
Chrom* roulleteSelection(Chrom *parent, Chrom *chroms){
	int i, j;
	double last = 0.0;
	double r;
	double sum = 0;
	double bigger = 0;
	RoulleteStruct *roullete = NULL;

 	for(i = 0; i < POPULATION_SIZE; i++){
		sum += fabs(chroms[i].fitness);
		bigger = fabs(chroms[i].fitness) > bigger ? fabs(chroms[i].fitness) : bigger;
	}
	
	//creating the roullete
	roullete = (RoulleteStruct *) malloc(sizeof(RoulleteStruct) * (POPULATION_SIZE + 1));

	roullete[0].value = 0.0;

	//the lowest the fitness, the closest to zero that gene got, then the roullete space will be biger for him
	for(i = 0; i < POPULATION_SIZE; i++){
		roullete[i + 1].value = last + bigger/sum - (fabs(chroms[i].fitness)/sum);
		roullete[i].chrom = chroms[i];
		last += bigger/sum - (fabs(chroms[i].fitness)/sum);
	}

	for(i = 0; i < POPULATION_SIZE; i++){
		r = RAND_DOUBLE(last);

		//CAN BE OVERWRITED BY BINARY SEARCH
		//selecting the parent
		for(j = 0; j < POPULATION_SIZE; j++){
			if(r >= roullete[j].value && r < roullete[j + 1].value){
				memcpy(&(parent[i]), &(roullete[j].chrom), sizeof(Chrom));
				break;
			}
		}
	}
	
	free(roullete);

	return parent;
}

void averageBetweenParentsCrossover(Chrom *parent, Chrom *chroms){
	int i, k;
	
	for(i = 0; i < POPULATION_SIZE; i++){
		if(RAND_DOUBLE(100) <= CROSSOVER_RATE){	
			//assign new value to the offspring
			for(k = 0; k < DIM; k++)
				chroms[i].genes[k] = (parent[i].genes[k] + parent[(i + 1)%POPULATION_SIZE].genes[k])/2.0;
		}
		//crossover didnt happened, both chromossomes pass to the next generation without modifications
		else{
			for(k = 0; k < DIM; k++){
				chroms[i].genes[k] = parent[i].genes[k];
			}
		}
	}
}

//aplies mutation to the genes
void sumMutation(Chrom *chroms){
	int i, k;

	for(i = 0; i < POPULATION_SIZE; i++){
		for(k = 0; k < DIM; k++){
			if(RAND_DOUBLE(100) <= MUTATION_RATE)
				chroms[i].genes[k] += rand()%2 == 0 ? MUTATION_SIZE : (-1.0)*MUTATION_SIZE;
		}
	}
}

//creates a new generation
Chrom* createGeneration(Chrom *chroms, ParentSelectionFunction f, CrossoverFunction c, MutationFunction m){
	Chrom *parent = NULL;

	parent = (Chrom *) malloc(sizeof(Chrom) * POPULATION_SIZE);

	//selecting parent
	parent = f(parent, chroms);

	//crossover
	c(parent, chroms);

	//mutating offspring
	m(chroms);
	
	free(parent);

	return chroms;
}

//randomly creates a initial population
Chrom *generateInitialPopulation(Chrom *chroms){
	int i, j;
	
	chroms = (Chrom *) malloc(sizeof(Chrom) * POPULATION_SIZE);

	srand(time(NULL));
	
	for(i = 0; i < POPULATION_SIZE; i++){
		for(j = 0; j < DIM; j++){
			chroms[i].genes[j] = RAND_DOUBLE(5);
		
			if(rand() % 2 == 0)
				chroms[i].genes[j] *= (-1.0);
		}
	}

	return chroms;
}

void printInformations(Chrom *chroms){
	printf("\nCrossover: %s.\n", CROSSOVER_TYPE);
	printf("Crossover rate: %d%s.\n", CROSSOVER_RATE, "%");
	printf("Mutation: +/- %lf.\n", MUTATION_SIZE);
	printf("Mutation rate: %d%s.\n", MUTATION_RATE, "%");
	printf("Parents selection: %s.\n", PARENTS_SELECTION);
	printf("Population size: %d.\n", POPULATION_SIZE);
	printf("Initialization: %s.\n", INITIALIZATION_TYPE);
	printf("Stop condition: %d %s.\n", CALLS_LIMIT, STOP_CONDITION);
	printf("Dimensions: %d.\n", DIM);
	printf("\n");
}

int main(int argc, char *argv[]){
	int counter = 0;
	double best = INT_MAX;
	double globalBest = INT_MAX;
	int i, j;
	Chrom globalChromBest;
	Chrom chromBest;
	Chrom *chroms = NULL;

	chroms = generateInitialPopulation(chroms);
	
	printInformations(chroms);

	while(++counter <= CALLS_LIMIT){	
		
//		printf("Generation %d\n", counter);

//		printGeneration(chroms);

		for(i = 0; i < POPULATION_SIZE; i++){
			chroms[i].fitness = ackleysFunction(chroms[i].genes, C1, C2, C3);

			if(fabs(chroms[i].fitness) < fabs(best)){
				best = chroms[i].fitness;
				chromBest = chroms[i];
			}

			//global minimum found
			if(fabs(chroms[i].fitness) == 0){
				printf("Chromossome genes: ");
				for(j = 0; j < DIM; j++)
					printf("%lf ", chroms[i].genes[j]);
				printf("- Iterations: %d - Ackely's Function = %.8lf\n", counter, chroms[i].fitness);
				free(chroms);
				return 0;
			}
		}

		if(fabs(best) < fabs(globalBest)){
			globalBest = best;
			globalChromBest = chromBest;
		}

		if(counter%1000 == 0){	
			printf("Best result at generation %d: %.8lf\n", counter, best);
			best = INT_MAX;
		}
		
		//starts a new generation to try to minimize the function
		chroms = createGeneration(chroms, roulleteSelection, averageBetweenParentsCrossover, sumMutation);
	}

	printf("Best overall result: %.8lf\n**Genes**\n" , globalBest);
	
	for(i = 0; i < DIM; i++)
	   	printf("%lf\n", globalChromBest.genes[i]);
	
	free(chroms);

	return 0;
}
