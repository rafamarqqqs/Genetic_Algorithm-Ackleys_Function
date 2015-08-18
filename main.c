#include <stdio.h>
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

#define C1 20.0
#define C2 0.20
#define C3 2.0*M_PI
#define CALLS_LIMIT 50000
#define POPULATION_SIZE 200
#define CROSSOVER_RATE 80		//max 100
#define MUTATION_SIZE 0.01
#define MUTATION_RATE 10		//max 100
#define RAND_DOUBLE(interval) ((double) rand() * interval / ((double) RAND_MAX))

//struct used to store data from a gene
typedef struct gene {
	double value;
	double rank;
} Gene;

//struct used to create the roullete
typedef struct roulleteStruct {
	double value;
	double gene;
} RoulleteStruct;

//function to specify the type of parent selection
typedef void (*ParentSelectionFunction)(Gene *, Gene *);

//calculates the ackley's function
double ackleysFunction(double x1, double x2, double c1, double c2, double c3){
	int i;
	double s1 = 0;
	double s2 = 0;
	int n = 2;

	for(i = 0; i < n; i++){
		s1 += (x1*x1) + (x2*x2);
		s2 += cos(c3*x1) + cos(c3*x2);
	}

	return (-1)*c1 * exp( (-1.0) * c2 * sqrt( (1.0/n) * s1)) - exp( (1.0/n) * s2) + c1 + 1;
}

int compareGene(const void *a, const void *b){
	Gene genA = *((Gene *) a);
	Gene genB = *((Gene *) b);
	return fabs(genA.rank) - fabs(genB.rank); 
}

//the best 20 parents will be always one of the parents
void bestParentsFitness(Gene *parents, Gene *genes){
	int i;

	qsort(parents, POPULATION_SIZE, sizeof(Gene), compareGene);

	//uses the 20 best parents to produce the new values, calculating the average with all genes values two by two
	for(i = 0; i < POPULATION_SIZE; i++){

		if(RAND_DOUBLE(100) <= CROSSOVER_RATE)
			genes[i].value = (parents[i % 20].value + parents[i].value)/2;
		else{
			genes[i].value = parents[i%20].value;
			genes[i + 1].value = parents[i].value;

			if(RAND_DOUBLE(100) <= MUTATION_RATE)
				genes[i].value += rand()%2 == 0 ? MUTATION_SIZE : (-1.0)*MUTATION_SIZE;
			if(RAND_DOUBLE(100) <= MUTATION_RATE)
				genes[i + 1].value += rand()%2 == 0 ? MUTATION_SIZE : (-1.0)*MUTATION_SIZE;
			
			i += 2;
		}
	
	}
}

//roullete type classification
void roulleteFitness(Gene *parents, Gene *genes){
	int i, j, k;
	double last = 0.0;
	double n;
	double parent[2];
	RoulleteStruct *roullete;

	//creating the roullete
	roullete = (RoulleteStruct *) malloc(sizeof(RoulleteStruct) * (POPULATION_SIZE + 1));

	roullete[0].value = 0.0;
	
	//5 - rank -> the lowest the rank, the closest to zero that gene got, then the roullete value will be biger for him
	for(i = 1; i <= POPULATION_SIZE; i++){
		roullete[i].value = last + (5 - fabs(parents[i].rank));
		roullete[i - 1].gene = parents[i].value;
		last += 5 - fabs(parents[i].rank);
	}
	 
	for(i = 0; i < POPULATION_SIZE; i++){
		n = RAND_DOUBLE(last);
		
		//CAN BE OVERWRITED BY BINARY SEARCH
		//selecting the parents
		for(j = 0, k = 0; j < POPULATION_SIZE && k < 2; j++){
			if(n >= roullete[j].value && n <= roullete[j + 1].value){
				parent[k] = roullete[j].gene;
				j = 0;
				k++;
				n = RAND_DOUBLE(last);
			}
		}
	
		if(RAND_DOUBLE(100) <= CROSSOVER_RATE){	
			//assign new value to the offspring
			genes[i].value = (parent[0] + parent[1])/2.0;
	
		}
		else{
			genes[i].value = parent[0];
			genes[i + 1].value = parent[1];

			if(RAND_DOUBLE(100) <= MUTATION_RATE)
				genes[i].value += rand()%2 == 0 ? MUTATION_SIZE : (-1.0)*MUTATION_SIZE;
			if(RAND_DOUBLE(100) <= MUTATION_RATE)
				genes[i + 1].value += rand()%2 == 0 ? MUTATION_SIZE : (-1.0)*MUTATION_SIZE;

			i += 2;
		}

	}

	free(roullete);
}

//fitness function
void fitness(Gene *genes, ParentSelectionFunction f){
	Gene *parents;

	parents = (Gene *) malloc(sizeof(Gene) * POPULATION_SIZE);
	
	//creating parents	
	memcpy(parents, genes, POPULATION_SIZE*sizeof(Gene));
	
	//selecting parents/creating offspring
	f(parents, genes);

	free(parents);
}

//randomly creates a initial population
void generateInitialPopulation(Gene *genes){
	int i;
	
	srand(time(NULL));
	
	for(i = 0; i < POPULATION_SIZE; i++){
		genes[i].value = RAND_DOUBLE(5);
		
		if(rand() % 2 == 0)
			genes[i].value *= (-1.0);
	}
}

void printInformations(Gene *genes){
	int i;	

	printf("CrossOver: %s.\n", CROSSOVER_TYPE);
	printf("Mutation: %s.\n", MUTATION);
	printf("Mutation rate: %s.\n", MUTATION_RATE_STR);
	printf("Parents selection: %s.\n", PARENTS_SELECTION);
	printf("Population size: %d.\n", POPULATION_SIZE);
	printf("Initialization: %s.\n", INITIALIZATION_TYPE);
	printf("Stop condition: %s.\n", STOP_CONDITION);

	printf("\n---- First Offspring ----\n");

	for(i = 0; i < POPULATION_SIZE; i++){
		if(i % 5 == 0)
			printf("\n");
		printf("%lf\t", genes[i].value);
	}	

	printf("\n");
}

int main(int argc, char *argv[]){
	int counter = 0;
	double best;
	double globalBest = 1000;
	double globalGeneBest;
	double geneBest;
	int i;
	Gene *genes;

	genes = (Gene *) malloc(sizeof(Gene) * POPULATION_SIZE);

	generateInitialPopulation(genes);
	
	printInformations(genes);

	while(++counter < CALLS_LIMIT){	
		best = 1000;

		for(i = 0; i < POPULATION_SIZE; i++){
			genes[i].rank = ackleysFunction(genes[i].value, genes[i].value, C1, C2, C3);
			
			if(fabs(genes[i].rank) < fabs(best)){
				best = genes[i].rank;
				geneBest = genes[i].value;
			}

			//global minima found
			if(genes[i].rank == 0.0){
				printf("Result: %lf - Iterations: %d\nAckely's Function = %lf\n", genes[i].value, counter, ackleysFunction(genes[i].value, genes[i].value, C1, C2, C3));
				free(genes);
				return 0;
			}
		}
	

		if(counter%1000 == 0)	
			printf("Best result at generation %d: %lf\n", counter, best);

		if(fabs(best) < fabs(globalBest)){
			globalBest = best;
			globalGeneBest = geneBest;
		}

		//fitness(genes, bestParentsFitness);
		fitness(genes, roulleteFitness);
	}

	printf("Best overall result: %lf - Gene: %lf.\n", best, globalGeneBest);

	free(genes);

	return 0;
}
