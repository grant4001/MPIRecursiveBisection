#include "mpi.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define NUM_POINTS 524288
#define NUM_QUADRANT_MAX 524288

unsigned int X_axis[NUM_POINTS];
unsigned int Y_axis[NUM_POINTS];
unsigned int Quadrants[NUM_QUADRANT_MAX][4];
double Cost; // Please place the total cost into the Cost variable on rank 0.

unsigned int Quadrants_temp[NUM_QUADRANT_MAX][4];

typedef struct 
{
  unsigned long X;
  unsigned long Y;
} Coord;

Coord Coords[NUM_POINTS];

int cmpfuncX(const void * a, const void * b)
{
  const Coord A = *(const Coord*) a;
  const Coord B = *(const Coord*) b;
  unsigned long a_x = (unsigned long)(A.X);
  unsigned long b_x = (unsigned long)(B.X);
  if (a_x < b_x)
  {
    return -1;
  }
  else if (a_x > b_x)
  {
    return 1;
  }
  else
  {
    return 0;
  }
  
}

int cmpfuncY(const void * a, const void * b)
{
  const Coord A = *(const Coord*) a;
  const Coord B = *(const Coord*) b;
  unsigned long a_y = (unsigned long)(A.Y);
  unsigned long b_y = (unsigned long)(B.Y);
  if (a_y < b_y) 
  {
    return -1;
  }
  else if (a_y > b_y) 
  {
    return 1;
  }
  else
  {
    return 0;
  }
}

// DO NOT MODIFY PROTOTYPE
void find_quadrants (num_quadrants)
     int num_quadrants;
{
  // Get MPI Info
  int numprocs, myid;
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);

  unsigned int X_min, X_max;
  unsigned int Quadrant_X_length, Quadrant_Y_length;
  double quad_start, quad_stop, part_start, part_stop;
  unsigned int i, j, k, q;

  //// Sorting  
  quad_start = MPI_Wtime();

  for (i = 0; i < NUM_POINTS; i++)
  {
    Coords[i].X = X_axis[i];
    Coords[i].Y = Y_axis[i];
  }

  // Find 1st quadrant
  qsort(&Coords[0], NUM_POINTS, sizeof(Coord), cmpfuncX);
  X_min = Coords[0].X; 
  X_max = Coords[NUM_POINTS-1].X;
  qsort(&Coords[0], NUM_POINTS, sizeof(Coord), cmpfuncY);
  Quadrants[0][0] = X_min;
  Quadrants[0][1] = Coords[0].Y;
  Quadrants[0][2] = X_max + 1;
  Quadrants[0][3] = Coords[NUM_POINTS-1].Y + 1;

  // Find subsequent quadrants
  for (q = 1; q < num_quadrants; q *= 2) 
  {
    for (i = 0; i < q; i++)
    {
      for (j = 0; j < 4; j++)
      {
        Quadrants_temp[i][j] = Quadrants[i][j];
      }
    } 

    for (i = 0; i < q; i++)
    {
      
      Quadrant_X_length = abs((long)Quadrants_temp[i][2] - (long)Quadrants_temp[i][0]);
      Quadrant_Y_length = abs((long)Quadrants_temp[i][3] - (long)Quadrants_temp[i][1]);
      

      /*
      Quadrant_X_length = (unsigned int)rand();
      Quadrant_Y_length = (unsigned int)rand();
      */
      if (Quadrant_X_length > Quadrant_Y_length)
      {
        qsort(&Coords[i*(NUM_POINTS/q)], (NUM_POINTS/q), sizeof(Coord), cmpfuncX);
        Quadrants[i*2][0] = Quadrants_temp[i][0];
        Quadrants[i*2][1] = Quadrants_temp[i][1];
        Quadrants[i*2][2] = Coords[i*(NUM_POINTS/q) + ((NUM_POINTS/q)/2) - 1].X + 1;
        Quadrants[i*2][3] = Quadrants_temp[i][3];
        Quadrants[i*2+1][0] = Coords[i*(NUM_POINTS/q) + ((NUM_POINTS/q)/2) - 1].X + 1;
        Quadrants[i*2+1][1] = Quadrants_temp[i][1];
        Quadrants[i*2+1][2] = Quadrants_temp[i][2];
        Quadrants[i*2+1][3] = Quadrants_temp[i][3];   
      }
      else 
      {
        qsort(&Coords[i*(NUM_POINTS/q)], (NUM_POINTS/q), sizeof(Coord), cmpfuncY);
        Quadrants[i*2][0] = Quadrants_temp[i][0];
        Quadrants[i*2][1] = Quadrants_temp[i][1];
        Quadrants[i*2][2] = Quadrants_temp[i][2];
        Quadrants[i*2][3] = Coords[i*(NUM_POINTS/q) + ((NUM_POINTS/q)/2) - 1].Y + 1;
        Quadrants[i*2+1][0] = Quadrants_temp[i][0];
        Quadrants[i*2+1][1] = Coords[i*(NUM_POINTS/q) + ((NUM_POINTS/q)/2) - 1].Y + 1;
        Quadrants[i*2+1][2] = Quadrants_temp[i][2];
        Quadrants[i*2+1][3] = Quadrants_temp[i][3];
      }
    }
  }
  quad_stop = MPI_Wtime();

  // Cost Calculation
  part_start = MPI_Wtime();
  int chunk = NUM_POINTS / num_quadrants;
  double X_dist, Y_dist;
  double Cost_local = 0.0;
  for (i = 0; i < num_quadrants; i++)
  {
    for (j = myid; j < chunk; j+=numprocs)
    {
      for (k = j + 1; k < chunk; k++)
      {
        X_dist = (double)Coords[i*chunk + j].X - (double)Coords[i*chunk + k].X;
        Y_dist = (double)Coords[i*chunk + j].Y - (double)Coords[i*chunk + k].Y;
        Cost_local += (double)sqrt(pow((double)X_dist, 2) + pow((double)Y_dist, 2));
        //if (Coords[i*chunk + j].X - Coords[i*chunk + k].X > pow(2,30)) printf("ERROR. %u \n", X_dist);
        
      }      
    }
  }
  MPI_Reduce(&Cost_local, &Cost, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  part_stop = MPI_Wtime();

  // Print results
  /*
  if (myid == 0)
  {
    
    printf("Quadrants computed (Top left, Bottom right):\n");
    
    
    for (i = 0; i < num_quadrants; i++)
    {
      printf("([%d, %d], [%d, %d])\n", Quadrants[i][0], Quadrants[i][1], Quadrants[i][2], Quadrants[i][3]);
    }
    
    printf("RESULTS:\n");
    printf("Total partition cost: %E\n", Cost);
    printf("Total execution time: %f\n", part_stop - quad_start);
    printf("\tQuadrants execution time: %f\n", quad_stop - quad_start);
    printf("\tCost execution time: %f\n", part_stop - part_start);
    
    
  }
  */

  // CHECK NO INTERSECTING QUADRANTS
  /*
  if (myid==0)
  {
    for (i = 0; i < num_quadrants; i++)
    {
      for (j = i+1; j < num_quadrants; j++)
      {
        if ((Quadrants[i][0] >= Quadrants[j][2]) || (Quadrants[j][0] >= Quadrants[i][2])) 
          continue;

        if ((Quadrants[i][1] >= Quadrants[j][3]) || (Quadrants[j][1] >= Quadrants[i][3]))
          continue;

        printf("Intersection! :( %d %d\n", i, j);
        printf("([%d, %d], [%d, %d])\n", Quadrants[i][0], Quadrants[i][1], Quadrants[i][2], Quadrants[i][3]);
        printf("([%d, %d], [%d, %d])\n", Quadrants[j][0], Quadrants[j][1], Quadrants[j][2], Quadrants[j][3]);

      }
    }
  }
  */

}

int main(argc,argv)
  int argc;
 char *argv[];
{
  int num_quadrants;
  int myid, numprocs;
  int  namelen;
  char processor_name[MPI_MAX_PROCESSOR_NAME];
    
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);
  MPI_Get_processor_name(processor_name,&namelen);

  if (argc != 2)
    {
      fprintf (stderr, "Usage: recursive_bisection <#of quadrants>\n");
      MPI_Finalize();
      exit (0);
    }

  fprintf (stderr,"Process %d on %s\n", myid, processor_name);

  num_quadrants = atoi (argv[1]);

  if (myid == 0)
    fprintf (stdout, "Extracting %d quadrants with %d processors \n", num_quadrants, numprocs);

  if (myid == 0)
    {
      int i;

      srand (10000);
      
      for (i = 0; i < NUM_POINTS; i++)
	X_axis[i] = (unsigned int)rand();

      for (i = 0; i < NUM_POINTS; i++)
	Y_axis[i] = (unsigned int)rand();
    }

  MPI_Bcast(&X_axis, NUM_POINTS, MPI_INT, 0, MPI_COMM_WORLD);  
  MPI_Bcast(&Y_axis, NUM_POINTS, MPI_INT, 0, MPI_COMM_WORLD);  

  find_quadrants (num_quadrants);
 
  MPI_Finalize();
  return 0;
}
  

