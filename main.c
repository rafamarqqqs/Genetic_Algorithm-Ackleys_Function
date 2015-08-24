#include <stdio.h>
#include <limits.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#define CROSSOVER_TYPE ("Average beetween parents")
#define MUTATION_RATE_STR ("10%")
#define PARENTS_SELECTION ("Roullete method")
#define INITIALIZATION_TYPE ("Aleatory")
#define STOP_CONDITION ("50.000 fitness calls OR optimum achived")
#define MUTATION ("+/- 0.01")

#define DIM 2
#define C1 20.0
#define C2 0.20
#define C3 2.0*M_PI
#define CALLS_LIMIT 50000
#define POPULATION_SIZE 200
#define CROSSOVER_RATE 70		//max 100
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
typedef void (*ParentSelectionFunction)(Chrom *, Chrom *);

//function to specify the type of crossover
typedef void (*CrossoverFunction)(Chrom *, Chrom *);

//function to specify the type of crossover
typedef void (*MutationFunction)(Chrom *);

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

int compareChrom(const void *a, const void *b){
	Chrom genA = *((Chrom *) a);
	Chrom genB = *((Chrom *) b);
	return fabs(genA.fitness) - fabs(genB.fitness); 
}

//roullete type classification
void roulleteSelection(Chrom *parents, Chrom *chroms){
	int i, j;
	double last = 0.0;
	double r;
	RoulleteStruct *roullete = NULL;
	
	//creating the roullete
	roullete = (RoulleteStruct *) malloc(sizeof(RoulleteStruct) * (POPULATION_SIZE + 1));

	roullete[0].value = 0.0;
	
	//5 - fitness -> the lowest the fitness, the closest to zero that gene got, then the roullete value will be biger for him
	for(i = 0; i < POPULATION_SIZE; i++){
		roullete[i + 1].value = last + (5 - fabs(chroms[i].fitness));
		roullete[i].chrom = chroms[i];
		last += 5 - fabs(chroms[i].fitness);
	}
	 
	for(i = 0; i < 2*POPULATION_SIZE; i++){
		r = RAND_DOUBLE(last);

		//CAN BE OVERWRITED BY BINARY SEARCH
		//selecting the parents
		for(j = 0; j < POPULATION_SIZE; j++){
			if(r > roullete[j].value && r <= roullete[j + 1].value){
				parents[i] = roullete[j].chrom;
				break;
				//r = RAND_DOUBLE(last);
			}
		}
	}
	
	free(roullete);
}

void averageBetweenParentsCrossover(Chrom *parent, Chrom *chroms){
	int i, k;
	
	for(i = 0; i < POPULATION_SIZE; i++){
		if(RAND_DOUBLE(100) <= CROSSOVER_RATE){	
			//assign new value to the offspring
			for(k = 0; k < DIM; k++)
				chroms[i].genes[k] = (parent[i].genes[k] + parent[i + POPULATION_SIZE - 1].genes[k])/2.0;
			}
		//crossover didnt happened, both chromossomes pass to the next generation without modifications
		else
			chroms[i] = parent[i];
	
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

void createGeneration(Chrom *chroms, ParentSelectionFunction f, CrossoverFunction c, MutationFunction m){
	Chrom *parents = NULL;

	parents = (Chrom *) malloc(sizeof(Chrom) * 2 * POPULATION_SIZE);
	
	//selecting parents
	f(parents, chroms);

	//crossover
	c(parents, chroms);

	//mutating offspring
	m(chroms);
	
	free(parents);
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

void printGeneration(Chrom *chroms){
	int i, j;

	for(i = 0; i < POPULATION_SIZE; i++){
		for(j = 0; j < DIM; j++)
			printf("%lf ", chroms[i].genes[j]);
		printf("\n");
	}
}

void printInformations(Chrom *chroms){
	printf("\nCrossOver: %s.\n", CROSSOVER_TYPE);
	printf("Mutation: %s.\n", MUTATION);
	printf("Mutation rate: %s.\n", MUTATION_RATE_STR);
	printf("Parents selection: %s.\n", PARENTS_SELECTION);
	printf("Population size: %d.\n", POPULATION_SIZE);
	printf("Initialization: %s.\n", INITIALIZATION_TYPE);
	printf("Stop condition: %s.\n", STOP_CONDITION);
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

		for(i = 0; i < POPULATION_SIZE; i++){
			chroms[i].fitness = ackleysFunction(chroms[i].genes, C1, C2, C3);
		
			//printf("%lf %lf %lf\n", chroms[i].genes[0], chroms[i].genes[1], chroms[i].fitness);

			if(fabs(chroms[i].fitness) < fabs(best)){
				best = chroms[i].fitness;
				chromBest = chroms[i];
			}

			//global minimum found
			if(fabs(chroms[i].fitness) == 0){
				printf("Chromossome genes: ");
				for(j = 0; j < DIM; j++)
					printf("%lf ", chroms[i].genes[j]);
				printf("- Iterations: %d - Ackely's Function = %lf\n", counter, chroms[i].fitness);
				free(chroms);
				return 0;
			}
		}

		if(fabs(best) < fabs(globalBest)){
			globalBest = best;
			globalChromBest = chromBest;
		}

		if(counter%1000 == 0){	
			printf("Best result at generation %d: %lf\n", counter, best);
			best = INT_MAX;
		}
		
		//starts a new generation to try to minimize the function
		createGeneration(chroms, roulleteSelection, averageBetweenParentsCrossover, sumMutation);
	}

	printf("Best overall result: %lf\n **Genes**\n " , globalBest);
	
	for(i = 0; i < DIM; i++)
	   	printf("%lf\n", globalChromBest.genes[i]);

	printf(".\n");

	free(chroms);

	return 0;
}
