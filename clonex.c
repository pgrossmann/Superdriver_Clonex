/*
  clonex.c

  Version : 0.0.02
  Author  : Niko Beerenwinkel
            Moritz Gerstung

  (c) Niko Beerenwinkel, Moritz Gerstung 2009
*/


#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <time.h>
#include <math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include <getopt.h>
#include <assert.h>


/* site for new mutation parameter */
#define d 220   // dimension (number of loci)
#define d0 100 // passengers
#define ds 10 // super drivers
#define dm 10 // mutator mutations 
// order of array holding types of mutations: (lowest index) drivers, superdrivers, mutators, mutators, passengers (highest index)
//#define M 3 // M = 3 or more/less?
#define MAX_GENOTYPES 12000000  // max. 1.2 x 10^7
#define MAX_K 220  // max. no. of mutations per genotype

/* CHECK:
 * 703
 * 489
 */
typedef int (*compfunc) (const void*, const void*);  // type for comparison function in qsort


struct Genotype
{
  int mutation[MAX_K+1];
  int k;  // number of mutations
  int k_s; //p number of superdriver mutations
  int k_m; //p number of mutator mutations
  int k_d; //p number of driver mutations
  int k_p; //p number of passenger mutations
//   double fit;
  long int count; //h
};


// struct Genotype *genotype;
// struct Genotype *genotype = malloc(sizeof(struct Genotype) * MAX_GENOTYPES);
// struct Genotype genotype[MAX_GENOTYPES];
// double geno_prob[MAX_GENOTYPES];
// unsigned int geno_count[MAX_GENOTYPES];
struct Genotype *genotype;
double *geno_prob;
unsigned int *geno_count;
long int N_g = 0;  // current number of genotypes (make sure that N_g <= MAX_GENOTYPES) //H
// double fitness[MAX_K+1];
// double *fitness = NULL;

// int sum_obs[d];
// int k_abs_freq[MAX_K+1];
// double k_rel_freq[MAX_K+1];
int *sum_obs;
int *k_abs_freq;
double *k_rel_freq;

int M = 3; //p

int outputIntermediate = 1;
//p note: index from 1, ... , !!!!
int lociEnd_d_d = d - d0 - dm - ds;
int lociEnd_d_s = d - d0 - dm;
int lociEnd_d_m = d - d0;

gsl_rng *RNG = NULL;  // random number generator


int hamming_distance(i1, i2)
{
  int j, l, g1[d], g2[d];

  for (j=0; j<d; j++)
    g1[j] = g2[j] = 0;

  for (l=0; l<genotype[i1].k; l++)
    g1[genotype[i1].mutation[l]] = 1;
  for (l=0; l<genotype[i2].k; l++)
    g1[genotype[i2].mutation[l]] = 1;

  int dist = 0;
  for (j=0; j<d; j++)
    dist += (g1[j] != g2[j]);

  return dist;
}


double expected_hamming_distance()
{
  long int i, i1, i2; //H

//   double rel_freq[MAX_GENOTYPES];
  double* rel_freq = (double*) malloc(sizeof(double)*MAX_GENOTYPES);
  assert(rel_freq);
//   if (!rel_freq)
//   {
//     printf
//   }
  int NN = 0;
  //printf("hellou4 %lu\n", N_g);
  for (i=0; i<N_g; i++) {
    NN += genotype[i].count;
//     printf("hellou3 %lu \n", i);
//     printf("hellou2 %d \n", genotype[i].count);
  }
  for (i=0; i<N_g; i++) {
//     printf("bevor hier: %g\n", 0/2);
//     printf("blbla");
//     printf("%d\n",2);
//     printf("hier1: %lu \n", N_g);
//     printf("hier2: %d \n", genotype[i].count);
//     printf("hier3: %lu \n", i);
//     printf("relfreg1 %g\n",rel_freq[i]);
    rel_freq[i] = (double) genotype[i].count / (double) NN;
//     printf("relfreg2 %g\n",rel_freq[i]);
  }
    
  
  double dist = 0.0;
  for (i1=0; i1<N_g; i1++)
    for (i2=i1+1; i2<N_g; i2++)
      dist += rel_freq[i1] * rel_freq[i2] * (double) hamming_distance(i1, i2);
  
    free(rel_freq);
  return dist;
}



double detect_mutations(int *k_min, int *k_max, double *k_mean, int *k_median, int *k_obs, double *k_plus, double *k_minus, double *prob_add, double *homo, double *div, double delta)
{
  int j, l;
  int k, count;
  long int i; //H

  for (j=0; j<d; j++)
    sum_obs[j] = 0;
  

  for (k=0; k<=MAX_K; k++)
    {
      k_abs_freq[k] = 0;
      k_rel_freq[k] = 0.0;
    }
  
  //#pragma omp parallel for private(k,count)
  for (i=0; i<N_g; i++)
    {
      k = genotype[i].k;
      count = genotype[i].count;

      k_abs_freq[k] += count;
      for (l=0; l<k; l++) {
	assert(i < MAX_GENOTYPES && i >= 0);
	assert(l < MAX_K+1 && l >= 0);
	assert(genotype[i].mutation[l] <= d && genotype[i].mutation[l] >= 0);
	sum_obs[genotype[i].mutation[l]] += count;
      }
	
    }

  int N = 0;
  int max_freq;
  for (k=0; k<=MAX_K; k++)
    {
      N += k_abs_freq[k];
      if ((k == 0) || (k_abs_freq[k] > max_freq))
	{
	  max_freq = k_abs_freq[k];
	  *k_median = k;
	}
    }
  

  int k_median_count = 0;
  int dominant_count = 0;
  int dominant_type;
  for (i=0; i<N_g; i++)
    {
      if (genotype[i].k == *k_median)
	{
	  k_median_count += genotype[i].count;
	  if ((dominant_count == 0) || (genotype[i].count > dominant_count))
	    {
	      dominant_type = i;
	      dominant_count = genotype[i].count;
	    }
	}
    }
  *homo = (double) dominant_count / (double) k_median_count;
 
  *div = -1.0; //expected_hamming_distance();

  for (k=0; k<=MAX_K; k++)
      k_rel_freq[k] = (double) k_abs_freq[k] / (double) N;

  *k_min = 0;
  while (k_abs_freq[*k_min] == 0)
    (*k_min)++;

  *k_max = MAX_K;
  while (k_abs_freq[*k_max] == 0)
    (*k_max)--;
  
  *k_mean = 0.0;
  for (k=0; k<=MAX_K; k++)
    *k_mean += k_rel_freq[k] * (double) k;

  int g_obs[d], g[d];
  for (j=0; j<d; j++)
    g_obs[j] = 0;
  
  *k_obs = 0;
  for (j=0; j<d; j++)
    if (sum_obs[j] > delta * (double) N)
      {
	(*k_obs)++;
	g_obs[j] = 1;
      }

  double E_dist = 0.0;
  *k_plus = *k_minus = 0.0;
  *prob_add = 0.0;
  
  double k_pp, k_mm, prob_aa;
  k_pp = k_mm = prob_aa = 0.0;
  
#pragma omp parallel for reduction(+: E_dist, k_pp, k_mm, prob_aa)
  for (i=0; i<N_g; i++)
    {
      for (j=0; j<d; j++)
	g[j] = 0;

      for (l=0; l<genotype[i].k; l++)
	g[genotype[i].mutation[l]] = 1;

      int dist = 0;
      int plus = 0;
      int minus = 0;
      for (j=0; j<d; j++)
	{
	  dist += (g[j] != g_obs[j]);
	  plus += (g[j] > g_obs[j]);  // mutation in cell i that is not observable in population
	  minus += (g[j] < g_obs[j]);  // mutation observable, but not present in cell i
	}

      double frac = (double) genotype[i].count / (double) N;
      E_dist += frac * (double) dist;
      k_pp += frac * (double) plus;
      k_mm += frac * (double) minus;
      if (plus > 0)
	prob_aa += frac;
    }
  
  *k_plus = k_pp;
  *k_minus = k_mm;
  *prob_add = prob_aa;
  

  /*
    double p = 0.0;
    for (k=(*k_obs)+1; k<=*k_max; k++)
      p += k_rel_freq[k];
  */

  return E_dist;
}



void summary()
{
  int k;
  long int i; //H

  printf("%s\t%s\t%s\t%s\n",
	 "genotypeIndex",
	 "number.mutations",
	 "genotype.count",
	 "indices.mutations");
  
  for (i=0; i<N_g; i++)
    {
      printf("%10lu ", i+1); 
      printf("%3d ", genotype[i].k);
      printf("%10lu -> ", genotype[i].count);
      for (k=0; k<genotype[i].k; k++)
	//for (k=0; k<MAX_K; k++)
	printf("%2d ", genotype[i].mutation[k]);
      printf("\n");
    }
  printf("\n");
  printf("%s\t%s\t%s\t%s\n",
	 "genotypeIndex",
	 "number.mutations",
	 "genotype.count",
	 "indices.mutations");
}



int count_cmp (const struct Genotype *g, const struct Genotype *h)
{
  return (h->count - g->count);  // decreasing
}


int no_of_mut_cmp (const struct Genotype *g, const struct Genotype *h)
{
  return (h->k - g->k);  // decreasing
}


int pos_zero_cmp (const struct Genotype *g, const struct Genotype *h)
{
  if ((h->count > 0) && (g->count == 0))
    return 1;
  if ((h->count == 0) && (g->count > 0))
    return -1;
  
  return 0;
}


void remove_zeros(int total_sort)
{

  if (!total_sort) //normally (0) only sort out zeros
    qsort(genotype, N_g, sizeof(struct Genotype), (compfunc) pos_zero_cmp);
  else
    qsort(genotype, N_g, sizeof(struct Genotype), (compfunc) count_cmp);
  
  while (genotype[N_g-1].count == 0)
    N_g--;

}



int int_cmp (const int *a, const int *b)
{
  return (*a - *b);  // increasing
}


int mutation_cmp (const struct Genotype *g, const struct Genotype *h)
{
  int j;

  if (h->k != g->k)
    return (h->k - g->k);  // decreasing
  else
    {
      for (j=0; j<g->k; j++)
	{
	  if (g->mutation[j] < h->mutation[j])
	    return -1;
	  if (g->mutation[j] > h->mutation[j])
	    return 1;
	}
      return 0;
    }


}


void remove_duplicates()
{
  long int i; //H

  for (i=0; i<N_g; i++)
    qsort(genotype[i].mutation, genotype[i].k, sizeof(int), (compfunc) int_cmp);  // sort mutation lists

  qsort(genotype, N_g, sizeof(struct Genotype), (compfunc) mutation_cmp);  // sort population

  for (i=0; i<N_g-1; i++)
    {
      if (mutation_cmp(&genotype[i], &genotype[i+1]) == 0)  // collect counts
	{
	  genotype[i+1].count += genotype[i].count;
	  genotype[i].count = 0;
	}
    }

  remove_zeros(0);

  
}

// to edit patrick: how adapt for mutator phenotype?
void print_pop(FILE* const DB1, int N, double s, double s_super, double s_mutator, int ggg, int rrr, char* filestem)
{
  int j, k, locus[d];
  long int i; //H
  double fit;
//   int nr_superdrivers, nr_drivers, nr_passengers, nr_mutators;
  //int nr_drivers;

//   char filename[255];
//   FILE *DB;
  
//   sprintf(filename, "%s/t%03d-g%05d.pop", filestem, rrr, ggg);
//     if ((DB = fopen(filename, "w")) == NULL)
//       {
// 	  fprintf(stderr, "In run %d: Cannot open output file -- %s\n", rrr, filename);
// 	  exit(1);
//        }
  
//   printf("A %p\n", DB);

  fprintf(DB1, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", 
	  "count",
	 "mean.count",
	 "number.mutations",
	 "fitness",
	 "genotype.driver_superdr_mutator_passenger",
	 "number.drivers",
	 "number.superdrivers",
	  "number.mutators",
	 "number.passengers"//,
// 	  "debug.number.drivers",
// 	 "debug.number.superdrivers",
// 	  "debug.number.mutators",
// 	 "debug.number.passengers"
 	);
  
//   printf("B %p\n", DB);
  
  for (i=0; i<N_g; i++)
    {
      for (j=0; j<d; j++)
	locus[j] = 0;
      for (k=0; k< genotype[i].k; k++)
	locus[genotype[i].mutation[k] - 1] = 1;

      // edit patrick: why are number of mutations, fitness, etc. re-counted and not stored in simulate() ?
      fit = 1.0;
      
      /*
      nr_superdrivers = 0;
      nr_drivers = 0;	
      nr_mutators = 0;
      nr_passengers = 0;
      */
      
      /*
      for (j=0;j<genotype[i].k;j++)
	{
	  if (genotype[i].mutation[j] > 0 && genotype[i].mutation[j] < d - d0) //p if not passenger mutation, //p check: if < or really <=
	    {  

	      if (genotype[i].mutation[j] < ds) // to edit patrick: "20" to "ds" and > to <=, //p check: if < or really <=
		fit *= 1 + s;

	      else

		fit *= 1 + s_super;

	    }

	} */
      
	for (j=0;j<genotype[i].k;j++)
	    { 
	      int currentLoci = genotype[i].mutation[j]; //p loci j of genotype i
	      //p /* site for new mutation parameter */
	      if (currentLoci > 0 && currentLoci < d - d0) //p if not passenger mutation, //p check: if < or really <=
		{  // edit patrick: where are passenger, still last 10?
		  if (currentLoci < d - d0 - dm - ds) //p check: if < or really <=
		    fit *= 1 + s;
		  else if (currentLoci >= d - d0 - dm - ds && currentLoci < d - d0 - dm)
		    fit *= 1 + s_super;
		  else
		    fit *= 1 + s_mutator;
		}
	      //      printf("%f\t%i\n", fit, genotype[i].k);
	    }  
	    
// 	    printf("ff %f, %f, %d \n", fit, genotype[i].fit, i);
// 	    assert(fit == genotype[i].fit);

      /*
      int inc_d = 0;
	  int inc_s = 0;
	  int inc_m = 0;
	  int inc_p = 0;
	  int kk;

// 	  printf("%d\n",genotype[i].k);
	  for (kk=0; kk < genotype[i].k; kk++) {
	    int tt = genotype[i].mutation[kk];
	    if (tt <= lociEnd_dhttp://www.youtube.com/watch?v=2i99AVsLznA_d)
	      inc_d++;
	    else if (tt > lociEnd_d_d && tt <= lociEnd_d_s) 
	      inc_s++;

	    else if (tt > lociEnd_d_s && tt <= lociEnd_d_m) 
	      inc_m++;
	    else
	      inc_p++;
// 	    printf("%d\n", inc_d);
	    genotype[i].k_d = inc_d;
	    genotype[i].k_s = inc_s;
	    genotype[i].k_m = inc_m;
	    genotype[i].k_p = inc_p;

	  }*/
      

	  
// 	  hellou 3125, 3001, 745, 62
// 0x179f85290
// 0

// 0x79f85290 
// A 0x79f85290
// B 0x79f85290
// hellou 2692, 201, 0, 117665468
// 0x179f85290
//       printf("hellou %lu, %d, %lu, %lu\n", N_g, ggg, i, genotype[i].count);
//       printf("0 %p\n", DB);
//       printf("1 %p\n", DB1);
//       printf("%p\n", DB);
//       fprintf(DB, "%d\n", 1);
//       printf("hellou %lu, %d, %lu, %lu\n", N_g, ggg, i, genotype[i].count);
      fprintf(DB1, "%lu\t", genotype[i].count);
      fprintf(DB1, "%f\t", (double)genotype[i].count/(double)N);
      fprintf(DB1, "%d\t", genotype[i].k);
      fprintf(DB1, "%f\t", fit);
      
      
      for (j=0; j<d; j++)
	{
	  if(j==d-d0-dm-ds || j==d-d0-dm || j==d-d0) 
	    fprintf(DB1,"|");
	  fprintf(DB1, "%d", locus[j]);
	}
	/*
      for (j=0;j<d-dm-d0-ds;j++) //p driver loci
	if (locus[j] == 1)
		nr_drivers = nr_drivers +1;
      for (j=d-dm-d0-ds; j<d-dm-d0; j++) //p superdriver loci
	if (locus[j] == 1)
		nr_superdrivers = nr_superdrivers + 1;
      for (j=d-dm-d0; j<d-d0; j++) //p mutator loci
	if (locus[j] == 1)
		nr_mutators = nr_mutators + 1;
      for (j=d-d0; j<d; j++) //p passenger loci
	if (locus[j] == 1)
		nr_passengers = nr_passengers + 1;
	
      

      printf("dd %d %d \n", nr_drivers, genotype[i].k_d);
      assert(nr_drivers == genotype[i].k_d);
      printf("ss %d %d \n", nr_superdrivers, genotype[i].k_s);
      assert(nr_superdrivers == genotype[i].k_s);
      printf("mm %d %d \n", nr_mutators, genotype[i].k_m);
      assert(nr_mutators == genotype[i].k_m);
      printf("pp %d %d \n", nr_passengers, genotype[i].k_p);
      assert(nr_passengers == genotype[i].k_p);

      fprintf(DB1, "\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", nr_drivers, nr_superdrivers, nr_mutators, nr_passengers,
	genotype[i].k_d, genotype[i].k_s, genotype[i].k_m, genotype[i].k_p);
	
      */

      fprintf(DB1, "\t%d\t%d\t%d\t%d\n", genotype[i].k_d, genotype[i].k_s, genotype[i].k_m, genotype[i].k_p);

    }
//     fclose(DB);
}

void addMutationCounts(int i) {
 
// 	  printf("%d",genotype[i].k_s);
	  int inc_d = 0;
	  int inc_s = 0;
	  int inc_m = 0;
	  int inc_p = 0;
	  int kk;
// 	  if (gen == 1800)
// 	  printf("%d\n",genotype[i].k);
	  assert(i < MAX_GENOTYPES && i >= 0);
	  for (kk=0; kk < genotype[i].k; kk++) {
	    assert(kk < MAX_K+1 && kk >= 0);
	    assert(i < MAX_GENOTYPES && i >= 0);
	    int tt = genotype[i].mutation[kk];
	    if (tt <= lociEnd_d_d)
	      inc_d++;
	    else if (tt > lociEnd_d_d && tt <= lociEnd_d_s) 
	      inc_s++;
	    else if (tt > lociEnd_d_s && tt <= lociEnd_d_m) 
	      inc_m++;
	    else
	      inc_p++;
// 	    if (gen == 1800)
// 	    printf("%d\n", inc_d);
	    assert(i < MAX_GENOTYPES && i >= 0);
	    genotype[i].k_d = inc_d;
	    genotype[i].k_s = inc_s;
	    genotype[i].k_m = inc_m;
	    genotype[i].k_p = inc_p;
	  } 
  
}

int simulate(FILE* DB, int N_init, int N_fin, int gen_max, double u, double u_mutator, 
	     double s, double s_super, double s_mutator, int run, int verbose, char* filestem)
{ // to edit patrick: parameter double u_m

  /*unsigned long*/ int gen;
  long int N; //H
  int j, c; //k;
  long int i; //H
  int mutation; //p, index_mut; //, mutation_number, free_slots;
  long int index_mut; //H
  int mut[MAX_K + 1], mut_k;//, mut_kd, mut_ks, mut_km, mut_kp;
  
  double p, N_exp_growth, p0, ps, pm;
  double u0 = u;
  // to edit patrick: double u_mutator = u_m; double u_orig = u;
  // double u_mutator = u_m;
  double u_orig = u;
  double k_mean, p_more, k_plus, k_minus, prob_add, homo, div;
  int k_min, k_max, k_median;
  int k_obs;
  char filename[255];
  FILE *DB2;
//   fitness = (double *) malloc(sizeof(double) * (MAX_K + 1));


  double a = exp( ( log(N_fin) - log(N_init) ) / gen_max );  // exponential growth rate
  /* ensures growth in gen_max generations from N_init to N_fin */
  
  //a = exp( log(2) / 60 );  // doubling time is 60 generations
  //a = exp( log(2) / 30 );  // doubling time is 30 generations

  printf("a = %f\n", a);
  printf("doubling time = %f generations\n", log(2) / log(a));


  /*
  for (j=0; j<=MAX_K; j++) {
    assert(j < MAX_K+1 && j >= 0);
    fitness[j] = pow(1.0 + s, (double) j); // edit patrick: fitness array what for?
  }
  */
    

  
  // initialize:
  N = N_init;
  N_exp_growth = (double) N;

  for (j=0; j<=MAX_K; j++) {
    assert(j < MAX_K+1 && j >= 0);
    genotype[0].mutation[j] = 0;
  }
    
  genotype[0].k = 0;
  //p initialize count for each distinct class of mutations
  genotype[0].k_m = 0; 
  genotype[0].k_s = 0;
  genotype[0].k_d = 0;
  genotype[0].k_p = 0;
  genotype[0].count = N;
//   genotype[0].fit = 1.0;
  N_g = 1;
  
  // grow population:
  gen = 0;
  do
    {
      gen++;
      assert(gen <= gen_max && gen >= 0);
//       printf("gen %d \n", gen);
      
      /* population growth */
      
      N_exp_growth *= a;
      N = (int) (N_exp_growth + 0.5);
      
      if (N > 2000000000)
	{
	  printf("Polyp has grown too large!\n");
	  exit(1);
	}
      
//       long int N_g_local = N_g; //H
      
      /* selection */
      double fit = 1.0;

      // compute probabilities:
      for (i=0; i<N_g; i++)
	{
	  assert(i < MAX_GENOTYPES && i >= 0);
	  for (j=0;j<genotype[i].k;j++)
	    { 
	      assert(j < MAX_K && j >= 0);
	      int currentLoci = genotype[i].mutation[j]; //p loci j of genotype i
	      //p /* site for new mutation parameter */
	      if (currentLoci > 0 && currentLoci < d - d0) //p if not passenger mutation, //p check: if < or really <=
		{  // edit patrick: where are passenger, still last 10?
		  if (currentLoci < d - d0 - dm - ds) //p check: if < or really <=
		    fit *= 1 + s;
		  else if (currentLoci >= d - d0 - dm - ds && currentLoci < d - d0 - dm)
		    fit *= 1 + s_super;
		  else
		    fit *= 1 + s_mutator;
		}
	      //      printf("%f\t%i\n", fit, genotype[i].k);
	    }
// 	  printf("ffs %f, %f \n", fit, genotype[i].fit);
// 	  assert(fit == genotype[i].fit);
	  //geno_prob[i] = fitness[genotype[i].k] * genotype[i].count;  // no need to normalize for gsl function
	  geno_prob[i] = fit * genotype[i].count;  // no need to normalize for gsl function
	  fit = 1.0; // reset fitness to 1.0 for next genotype.
	}
      
      
      for (i=0; i<N_g; i++) {
	assert(i < MAX_GENOTYPES && i >= 0);
	geno_count[i] = genotype[i].count;
      }
	
      
      gsl_ran_multinomial(RNG, N_g, N, geno_prob, geno_count);
      
      for (i=0; i<N_g; i++) {
	assert(i < MAX_GENOTYPES && i >= 0);
	genotype[i].count = geno_count[i];
      }
	
      
      remove_zeros(0);  // because low-frequency mutants are likely not to be sampled
      // and we need to consider less genotypes for mutation
      
      
      /* mutation */
      
      for (i=0; i<N_g; i++)
	{ // edit patrick: super drivers also here?
	  // to edit patrick: if (genotype[i].k >= M) u = u_mutator; and afterwars u = u_init
	  
	  addMutationCounts(i);
	  
	  assert(i < MAX_GENOTYPES && i >= 0);
	  if (genotype[i].k_m >= M) //p /* site for mutator phenotype */
	    u = u_mutator;
	    
	  p = 1.0 - gsl_ran_binomial_pdf(0, u, d - d0 - ds - dm);  // prob. of >= 1 mutation //p prob. for driver mut.
	  assert(i < MAX_GENOTYPES && i >= 0);
	  int N_mut_cells = gsl_ran_binomial(RNG, p, genotype[i].count);  // number of mutated cells //p with current genotype and new driver mut.
 	  p0 = 1.0 - gsl_ran_binomial_pdf(0, u0, d0);  // prob. of >= 1 pass. mutation
	  assert(i < MAX_GENOTYPES && i >= 0);
	  int N_mut_cells0 = gsl_ran_binomial(RNG, p0, genotype[i].count);  // number of mutated cells
	  ps = 1.0 - gsl_ran_binomial_pdf(0, u, ds);  //p prob. of >= 1 superdriver mutation
	  assert(i < MAX_GENOTYPES && i >= 0);
	  int N_mut_cells_s = gsl_ran_binomial(RNG, ps, genotype[i].count);  //p number of mutated cells with new superdriver mut.
	  pm = 1.0 - gsl_ran_binomial_pdf(0, u, dm);  //p prob. of >= 1 mutator mutation
	  assert(i < MAX_GENOTYPES && i >= 0);
	  int N_mut_cells_m = gsl_ran_binomial(RNG, pm, genotype[i].count);  //p number of mutated cells with new mutator mut.
	  
// 	  genotype[i].fit = pow((double) 1+s, N_mut_cells) * pow((double) 1+s_super, N_mut_cells_s) * pow((double) 1+s_mutator, N_mut_cells_m);
	  int n_new_mutations = N_mut_cells + N_mut_cells0 + N_mut_cells_s + N_mut_cells_m; //p comp once
	  //p each new mutation will increase number of distinct genotypes by 1
	  if (N_g + n_new_mutations > MAX_GENOTYPES)
	    {
	      printf("Too many genotypes: out of memory\n");
	      printf("MAX_GENOTYPES = %d\n", MAX_GENOTYPES);
	      printf("N_g = %lu\n", N_g);
	      printf("n_new_mutations = %d\n", n_new_mutations);
	      printf("gen = %d\n", gen);
	      summary();
	      
// 	      if(gen%10==1 && outputIntermediate)
// 		{
		  sprintf(filename, "%s/t%03d-g%05d_failedRun.pop", filestem, run, gen);
		   if ((DB2 = fopen(filename, "w")) == NULL)
		    {
	  		fprintf(stderr, "In run %d: Cannot open output file -- %s\n", run, filename);
	  		exit(1);
		    }
		  print_pop(DB2,N, s, s_super, s_mutator, gen, run, filestem);
		  fclose(DB2);
// 		} 
	
	      exit(1);
	    }
	  
	  //printf("%i\n", N_g);
	  
	  
	  for (c = 0; c < n_new_mutations; c++) //p for all new mutations
	    {
// 	      printf("blub %lu \n", i);
	      // make a copy of genotype i for mutation, decrease count: // edit patrick: why?
	      assert(i < MAX_GENOTYPES && i >= 0);
	      (genotype[i].count)--;
	      assert(i < MAX_GENOTYPES && i >= 0);
	      if (genotype[i].count == 0)
		index_mut = i;  // overwrite
	      else
		{
		  index_mut = N_g;
		  assert(i < MAX_GENOTYPES && i >= 0);
		  assert(index_mut < MAX_GENOTYPES && index_mut >= 0);
		  genotype[index_mut] = genotype[i];  // copy
		  N_g++;
		}
	      assert(index_mut < MAX_GENOTYPES && index_mut >= 0);
	      genotype[index_mut].count = 1;
	      /*
	      120259084314
	      219043332144
	      42949672967 bei .count--
	      */
// 	      printf("blub2 %lu, %lu \n", i, index_mut);
	      
	      assert(index_mut < MAX_GENOTYPES && index_mut >= 0);
	      if (genotype[index_mut].k - 2 > MAX_K) //p why - 2? refers to index?
		{
		  printf("Warning: More than %d mutations generated!\n", MAX_K);
		  assert(index_mut < MAX_GENOTYPES && index_mut >= 0);
		  genotype[index_mut].k = MAX_K - 2;  // reset
		  /*
		    printf("Too many mutations: out of memory\n");
		    printf("MAX_K = %d\n", MAX_K);
		    exit(1);
		  */
		}
// 	      printf("blub3 %lu, %lu \n", i, index_mut);
	      // edit patrick: what does add first and second mutation mean?
	      // add first mutation:
	      double a, b;
	      // edit patrick: super drivers also here? what is this block for
	      
	      /* outsource c counter booleans for re-usage */
		
	      int driver_bool = c < N_mut_cells;
	      int passenger_bool = c >= N_mut_cells && c < n_new_mutations - N_mut_cells_s - N_mut_cells_m;
	      int superdriver_bool = c >= N_mut_cells + N_mut_cells0 && c < n_new_mutations - N_mut_cells_m;
	      
// 	      printf("blub4 %lu, %lu \n", i, index_mut);
	      
	      //p Note: order of mutation types: drivers, superdrivers, passengers
	      if (driver_bool) //p for each driver mutations
		{
		  a = 1.0;
		  b = (double)d - (double)dm - (double)d0 - (double)ds + 1.0;
// 		  (genotype[index_mut].k_d)++;
		}
	      else if (passenger_bool) //p for each passenger mutations
		{
		  a = (double)d - (double)d0 + 1.0;
		  b = (double)d + 1.0;
// 		  (genotype[index_mut].k_p)++;
		}
		//else if (c >= n_new_mutations - N_mut_cells_s - N_mut_cells_m && c < n_new_mutations - N_mut_cells_m)//p for each superdriver mutation
		else if (superdriver_bool)
		{
		  a = (double)d - (double)dm - (double)d0 - (double)ds + 1.0;
		  b = (double)d - (double)dm - (double)d0 + 1.0;
// 		  (genotype[index_mut].k_s)++;
		} 
		else //p fo reach mutator mutation, should be c >= N_mut_cells + N_mut_cells0 + N_mut_cells_s 
		{
		  a = (double)d - (double)dm - (double)d0 + 1.0;
		  b = (double)d - (double)d0 + 1.0;
// 		  (genotype[index_mut].k_m)++;
		}
	      mutation = (int) floor(gsl_ran_flat(RNG, a, b));
	      

// 	      printf("blub5 %lu, %lu \n", i, index_mut);
	      
	      //printf("%i\n", mutation);
	      
	      assert(index_mut < MAX_GENOTYPES && index_mut >= 0);
	      assert(genotype[index_mut].k < MAX_K+1 && genotype[index_mut].k >= 0);
	      
	      genotype[index_mut].mutation[genotype[index_mut].k] = mutation;
	      
	      assert(index_mut < MAX_GENOTYPES && index_mut >= 0);
	      (genotype[index_mut].k)++;
	      
// 	      printf("blub6 %lu, %lu \n", i, index_mut);
	      
	      // maybe add a second mutation:
	      if (driver_bool) //p driver
		p = gsl_ran_binomial_pdf(0, u, d - d0 - ds - dm);
	      else if (passenger_bool) //p passenger mutations
		p = gsl_ran_binomial_pdf(0, u0, d0);
	      else if (superdriver_bool) //p superdriver
		p = gsl_ran_binomial_pdf(0, u, ds);	
	      else //p mutator
		p = gsl_ran_binomial_pdf(0, u, dm);	
				     
	     // prob of zero additional mutations
	      if (gsl_ran_flat(RNG, 0.0, 1.0) > p)  // i.e., additional mutations, assume exactly 1
		{
		  mutation = (int) floor(gsl_ran_flat(RNG, a, b));
		  //printf("%i\n", mutation);
		  assert(index_mut < MAX_GENOTYPES && index_mut >= 0);
		  assert(genotype[index_mut].k < MAX_K+1 && genotype[index_mut].k >= 0);
		  
		  genotype[index_mut].mutation[genotype[index_mut].k] = mutation;
		  
		  assert(index_mut < MAX_GENOTYPES && index_mut >= 0);
		  (genotype[index_mut].k)++;
		  /*
		  if (driver_bool)
		    (genotype[index_mut].k_d)++;
		  else if (passenger_bool)
		    (genotype[index_mut].k_p)++;
		  else if (superdriver_bool)
		    (genotype[index_mut].k_s)++;
		  else
		    (genotype[index_mut].k_m)++; */
		  
		}		      
	      
// 	      printf("blub7 %lu, %lu \n", i, index_mut);
	      // edit patrick: why mutations duplicated with this script?
	      // remove duplicate mutations (mutation[] is just a list!):
	      assert(index_mut < MAX_GENOTYPES && index_mut >= 0);
	      if (genotype[index_mut].k > 1)
		{
// 		  printf("blub771 %lu, %lu \n", i, index_mut);
		  assert(index_mut < MAX_GENOTYPES && index_mut >= 0);
		  qsort(genotype[index_mut].mutation, genotype[index_mut].k, sizeof(int), (compfunc) int_cmp);
		  mut[0] = genotype[index_mut].mutation[0];
		  mut_k = 1;
		  /*
		  mut_kd = 1;
		  mut_ks = 1;
		  mut_km = 1;
		  mut_kp = 1;
		  */
		  assert(index_mut < MAX_GENOTYPES && index_mut >= 0);
// 		  printf("blub77 %lu, %lu, %lu, %d, %d, %d\n", i, N_g, index_mut, genotype[index_mut].k, j, gen); // da
		  for (j=1; j<genotype[index_mut].k; j++)
		    {	
// 		      printf("blub777 %lu, %lu \n", i, index_mut);
		      assert(index_mut < MAX_GENOTYPES && index_mut >= 0);
		      assert(j < MAX_K+1 && j > 0);
		      if (genotype[index_mut].mutation[j] != genotype[index_mut].mutation[j-1])
			{
// 			  printf("blub111 %lu, %lu \n", i, index_mut);
			  assert(mut_k<MAX_K+1 && mut_k>0);
			  mut[mut_k] = genotype[index_mut].mutation[j];
			  //int tmp_mut = mut[mut_k];
			  mut_k++;
			  /*
			  if (tmp_mut <= lociEnd_d_d)
			    mut_kd++;
			  else if (tmp_mut > lociEnd_d_d && tmp_mut <= lociEnd_d_s)
			    mut_ks++;
			  else if (tmp_mut > lociEnd_d_s && tmp_mut <= lociEnd_d_m)
			    mut_km++;
			  else
			    mut_kp++; */
// 			  printf("blub222 %lu, %lu \n", i, index_mut);
			}
// 			printf("blub7777 %lu, %lu \n", i, index_mut);
			assert(index_mut < MAX_GENOTYPES && index_mut >= 0);		
		    }
// 		    printf("blub8 %lu, %lu \n", i, index_mut);
		  // overwrite:
		  for (j=0; j<mut_k; j++) {
		    assert(index_mut < MAX_GENOTYPES && index_mut >= 0);
		    assert(j < MAX_K+1 && j >= 0);
		    
		    genotype[index_mut].mutation[j] = mut[j];
		    /*
		    int tmp_mut = mut[j];
		    if (tmp_mut <= lociEnd_d_d)
			    mut_kd++;
			  else if (tmp_mut > lociEnd_d_d && tmp_mut <= lociEnd_d_s)
			    mut_ks++;
			  else if (tmp_mut > lociEnd_d_s && tmp_mut <= lociEnd_d_m)
			    mut_km++;
			  else
			    mut_kp++; */
		    
		  }
		  assert(index_mut < MAX_GENOTYPES && index_mut >= 0);
		  genotype[index_mut].k = mut_k;
		  /*
		  genotype[index_mut].k_d = mut_kd;
		  genotype[index_mut].k_s = mut_ks;
		  genotype[index_mut].k_m = mut_km;
		  genotype[index_mut].k_p = mut_kp; */

		}

	      
	    }
	    
	    u = u_orig;
	    
	    
// 	    printf("%lu, %lu, %lu, %d, %d\n", i, N_g_local, N_g, genotype[i].k, gen); // da
// 	    if (genotype[i].k - 2 > MAX_K) //p why - 2? refers to index?
// 		{
// 		  printf("Warning: More than %d mutations generated!\n", MAX_K);
// 		  genotype[i].k = MAX_K - 2;  // reset
// 		/*
// 		    printf("Too many mutations: out of memory\n");
// 		    printf("MAX_K = %d\n", MAX_K);
// 		    exit(1);
// 		  */
// 		}
	    
	    addMutationCounts(i);
	    
	}
      
      //p_more = detect_mutations(&k_min, &k_max, &k_mean, &k_obs, 0.5);
      
      remove_zeros(0);
      
      //fprintf(DB, "%d\t%d\t%d\t%f\t%d\t%d\t%f\n", gen, N, k_min, k_mean, k_max, k_obs, p_more);

      
      if (gen % 10 == 0)
	{
	  remove_duplicates();
	  //summary();
	}
      
      if (gen == gen_max)
	{
// 	  printf("before remove duplicates: %d\n", genotype[1].count);
	  remove_duplicates();
	  remove_zeros(1);
	  summary();
	  
	  //for (j=0; j<d; j++)
	  //printf("%d\n", sum_obs[j]);
// 	  printf("before detect mutations: %lu\n", N_g);
// 	  printf("before detect mutations2: %d\n", genotype[1].count);
	  p_more = detect_mutations(&k_min, &k_max, &k_mean, &k_median, &k_obs, &k_plus, &k_minus, &prob_add, &homo, &div, 0.5);
// 	  printf("nach detect mutations: %lu\n", N_g);
	  fprintf(DB, "%d\t%d\t%d\t%lu\t%g\t%g\t%d\t%g\t%g\t%d\t%g\t%d\t%d\t%d\t%g\t%g\t%g\t%g\t%g\t%g\n", run, gen, N_init, N, u, u_mutator, M, s, s_super, k_min, k_mean, k_median, k_max, k_obs, p_more, k_plus, k_minus, prob_add, homo, div);
// 	  printf("nach detect mutations2: %lu\n", N_g);
	  
	}
	
// 	printf("hellou1 %d \n",gen);
     if(gen%10==1 && outputIntermediate)
	{
	  remove_duplicates();
	  sprintf(filename, "%s/t%03d-g%05d.pop", filestem, run, gen);
// 	     sprintf(filename, "%s/test_%d_%d.pop", filestem, run, gen);
	  if ((DB2 = fopen(filename, "w")) == NULL)
	    {
	      fprintf(stderr, "In run %d: Cannot open output file -- %s\n", run, filename);
	      exit(1);
	    }
			
// 	  printf("hellou1 %d \n",gen); // da
// 	  printf("%p \n",DB2);
	  print_pop(DB2,N, s, s_super, s_mutator, gen, run, filestem);
	  fclose(DB2);
	}
	
    }
  while(gen < gen_max);
  
      
  
  return 0;
}



int main(int argc, char **argv)
{
    
    // defaults:
  int     N = 1000000000;  // population size
  int N_ini = 1;
//   int n = N_ini;
  //  int    n_d = 100;  // drivers
  //  int    n_p = 100;  // passengers
  double u = 1e-8;  // mutation rate
  // to edit patrick: double u_mutator = 1e-5 // what rate?
  double u_mutator = 1e-5; //p or in dependency of u
  double s = 1e-2;  // Selective advantage
  double superdriverFactor = 1.5;
  double s_super; // = superdriverFactor*s; // Selective advantage of super drivers // edit patrick: should we really make s_super dependent on s?
  // edit patrick: Datta 2012 defines maximums of types of mutations, we too?
  double s_mutator = 0;

  int    g = 1800;  // number of Generations
  int    R = 1;  // number of simulations Runs
  char *filestem;
  unsigned int seed = (unsigned) time(NULL);  // r, random seed // edit patrick: reproducibility?
  unsigned int userSeed = 0;
  int verbose = 0;
  
  int error_flag = 0;
  int f_flag = 0;
  
  int c = 0;
  while((c = getopt(argc, argv, "N:n:u:w:M:s:t:c:l:g:R:f:r:vh")) != EOF )
    {
      switch(c)  // to edit patrick: add u_mutator
        {
	case 'N':
	  if (atoi(optarg) > 0)
	    N = atoi(optarg);
	  else
	    error_flag++;
	  break;
	  
	case 'n':
	  if (atoi(optarg) > 0)
	    N_ini = atoi(optarg);
	  else
	    error_flag++;
	  break;  
	  
	case 'u':
	  if (atof(optarg) > 0)
	    u = atof(optarg);
	  else
	    error_flag++;
	  break;
	  
	  case 'w':
	  if (atof(optarg) > 0)
	    u_mutator = atof(optarg);
	  else
	    error_flag++;
	  break;
	  
	  case 'M':
	  if (atof(optarg) > 0)
	    M = atof(optarg);
	  else
	    error_flag++;
	  break;
	  
	case 's':
	  if (atof(optarg) >= 0)
	    s = atof(optarg);
	  else
	    error_flag++;
	  break;

// 	case 't':
// 	  if (atof(optarg) >= 0)
// 	    s_super = atof(optarg);
// 	  else
// 	    error_flag++;
// 	  break;
	  
	case 'c':
	  if (atof(optarg) >= 0)
	    superdriverFactor = atof(optarg);
	  else
	    error_flag++;
	  break;
	 
	case 'l':
	  if (atof(optarg) >= 0)
	    s_mutator = atof(optarg);
	  else
	    error_flag++;
	  break;
	  
	case 'g':
	  if (atoi(optarg) > 0)
	    g = atoi(optarg);
	  else
	    error_flag++;
	  break;
	  
	case 'R':
	  if (atoi(optarg) >= 0)
	    R = atoi(optarg);
	  else
	    error_flag++;
	  break;
	  
	case 'f':
	  filestem = optarg;
	  f_flag++;
	  break;
	  
	case 'r':
	  seed = atoi(optarg);
	  break;
	  
	case 'v':
	  verbose = 1;
	  break;
	  
	case 'h':
	  printf("usage: clonex [-NnuwMslcgRrh] -f <directory> \n");
	  printf("  N - Final population size (default = %d)\n", N);
	  printf("  n - Initial population size (default = %d)\n", N_ini);
	  printf("  u - Mutation rate (default = %g)\n", u); // to edit patrick: add u_mutator
	  printf("  w - Mutation rate for 'mutator phenotype' (default = %g)\n", u_mutator);
	  printf("  M - Number of mutations needed to evolve to 'mutator phenotype' (default = %d)\n", M);
	  printf("  s - Selective advantage (default = %g)\n", s);
// 	  printf("  t - >0 if selective advantage of 'super drivers' should be turned on (default = %g)\n", s_super);
	  printf("  l - Selective advantage of 'mutator mutations' (default = %g)\n", s_mutator);
	  printf("  c - Factor multiplied to s to obtain final t (default = %g)\n", superdriverFactor);
	  printf("  g - Number of generations (default = %d)\n", g);
	  printf("  R - Replicates (default = %d)\n", R);
	  printf("  r - Random seed (default = time)\n");
	  printf("  f - File directory (Required! Make sure that the directory exists!)\n");
	  printf("  h - This help\n");
	  exit(0);
	  
	default :
	  exit(1);
        }
    }
    
//     if (s_super > 0)
      s_super = superdriverFactor*s;
  
  if (error_flag || (! f_flag))
    {
      fprintf(stderr, "Error!\n");
      exit(1);
    }
  
  // random number generator:
  RNG = gsl_rng_alloc (gsl_rng_taus);  // global variable
  printf("rng %p \n", RNG);
//   gsl_rng_set(RNG, seed);  // seed rng
  if (gsl_rng_max(RNG) - gsl_rng_min(RNG) < N)
    {
      printf("Population size N = %d too large for random number generator!\n", N);
      printf("RNG range = [%lu, %lu]\n", gsl_rng_min(RNG), gsl_rng_max(RNG));
      exit(1);
    }
  
  
  char summary_filename[255], filename[255];
  sprintf(summary_filename, "%s/summary.tab", filestem);
  FILE *DB;
  if ((DB = fopen(summary_filename, "w")) == NULL)
    {
      fprintf(stderr, "Cannot open output file -- %s\n", summary_filename);
      exit(1);
    }
  fprintf(DB, "Run\tGen\tNinit\tNfin\tu\tu_m\tM\ts\ts_s\tMin\tMean\tMedian\tMax\tObs\tDist\tPlus\tMinus\tProbAdd\tHomo\tDiv\n");
  fclose(DB);
  
  
  int r;
  //#pragma omp parallel for private(DB, genotype, geno_prob, geno_count, N_g, fitness, sum_obs, k_abs_freq, k_rel_freq)
//   printf("%d\n", R);
  for (r=0; r<R; r++)
    {
      if (!userSeed)
	seed = (unsigned) time(NULL);
      else
	seed = userSeed;
      genotype = malloc(sizeof(struct Genotype) * MAX_GENOTYPES);
      assert(genotype);
      geno_prob = malloc(sizeof(double) * MAX_GENOTYPES);
      assert(geno_prob);
      geno_count = malloc(sizeof(unsigned int) * MAX_GENOTYPES);
      assert(geno_count);
      sum_obs = malloc(sizeof(int) * (d+1));
      assert(sum_obs);
      k_abs_freq = malloc(sizeof(int) * (MAX_K+1));
      assert(k_abs_freq);
      k_rel_freq = malloc(sizeof(double) * (MAX_K+1));
      assert(k_rel_freq);
      gsl_rng_set(RNG, seed);
      if ((DB = fopen(summary_filename, "a")) == NULL)
	{
	  fprintf(stderr, "In run %d: Cannot open output file -- %s\n", r+1, summary_filename);
	  exit(1);
	}
// 	printf("%g\n", s_super);
      simulate(DB, N_ini, N, g, u, u_mutator, s, s_super, s_mutator, r+1, verbose, filestem);
      // edit patrick: do we start with 1 or 10^6 cells?
      fclose(DB);
      
      //printf("Succeed.\n");
      sprintf(filename, "%s/r%03d.pop", filestem, r+1);
      if ((DB = fopen(filename, "w")) == NULL)
	{
	  fprintf(stderr, "In run %d: Cannot open output file -- %s\n", r+1, filename);
	  exit(1);
	}
      print_pop(DB,N, s, s_super, s_mutator, g, r, filestem);
      fclose(DB);
      free(genotype);
      free(geno_prob);
      free(geno_count);
      free(sum_obs);
      free(k_abs_freq);
      free(k_rel_freq);
    }
  

  return 0;
}



