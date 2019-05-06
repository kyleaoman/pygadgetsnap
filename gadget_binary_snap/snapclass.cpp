#include "snapclass.h"

#include <iostream>

using namespace std;

#ifdef SWIG
Snap S = Snap();
#endif

void Snap::load_snapshot(char* fname, int files)
{
  FILE *fd;
  char buf[200];
  int i, j, k, dummy, ntot_withmasses;
  int t, n, off, pc, pc_new, pc_sph;

#define SKIP fread(&dummy, sizeof(dummy), 1, fd);

  for(i = 0, pc = 1; i < files; i++, pc = pc_new)
    {
      if(files > 1)
	sprintf(buf, "%s.%d", fname, i);
      else
	sprintf(buf, "%s", fname);

      if(!(fd = fopen(buf, "r")))
	{
	  printf("can't open file `%s`\n", buf);
	  exit(0);
	}

      //printf("reading `%s' ...\n", buf);
      fflush(stdout);

      fread(&dummy, sizeof(dummy), 1, fd);
      fread(&header1, sizeof(header1), 1, fd);
      fread(&dummy, sizeof(dummy), 1, fd);

      if(files == 1)
	{
	  for(k = 0, NumPart = 0, ntot_withmasses = 0; k < 6; k++)
	    NumPart += header1.npart[k];
	  Ngas = header1.npart[0];
	}
      else
	{
	  for(k = 0, NumPart = 0, ntot_withmasses = 0; k < 6; k++)
	    NumPart += header1.npartTotal[k];
	  Ngas = header1.npartTotal[0];
	}

      for(k = 0, ntot_withmasses = 0; k < 6; k++)
	{
	  if(header1.mass[k] == 0)
	    ntot_withmasses += header1.npart[k];
	}

      if(i == 0)
	{
	  allocated = allocate_memory();
	}

      SKIP;
      for(k = 0, pc_new = pc; k < 6; k++)
	{
	  for(n = 0; n < header1.npart[k]; n++)
	    {
	      fread(&P[pc_new].Pos[0], sizeof(float), 3, fd);
	      pc_new++;
	    }
	}
      SKIP;

      SKIP;
      for(k = 0, pc_new = pc; k < 6; k++)
	{
	  for(n = 0; n < header1.npart[k]; n++)
	    {
	      fread(&P[pc_new].Vel[0], sizeof(float), 3, fd);
	      pc_new++;
	    }
	}
      SKIP;


      SKIP;
      for(k = 0, pc_new = pc; k < 6; k++)
	{
	  for(n = 0; n < header1.npart[k]; n++)
	    {
	      fread(&Id[pc_new], sizeof(int), 1, fd);
	      pc_new++;
	    }
	}
      SKIP;


      if(ntot_withmasses > 0)
	SKIP;
      for(k = 0, pc_new = pc; k < 6; k++)
	{
	  for(n = 0; n < header1.npart[k]; n++)
	    {
	      P[pc_new].Type = k;

	      if(header1.mass[k] == 0)
		fread(&P[pc_new].Mass, sizeof(float), 1, fd);
	      else
		P[pc_new].Mass = header1.mass[k];
	      pc_new++;
	    }
	}
      if(ntot_withmasses > 0)
	SKIP;


      if(header1.npart[0] > 0)
	{
	  SKIP;
	  for(n = 0, pc_sph = pc; n < header1.npart[0]; n++)
	    {
	      fread(&P[pc_sph].U, sizeof(float), 1, fd);
	      pc_sph++;
	    }
	  SKIP;

	  SKIP;
	  for(n = 0, pc_sph = pc; n < header1.npart[0]; n++)
	    {
	      fread(&P[pc_sph].Rho, sizeof(float), 1, fd);
	      pc_sph++;
	    }
	  SKIP;

	  if(header1.flag_cooling)
	    {
	      SKIP;
	      for(n = 0, pc_sph = pc; n < header1.npart[0]; n++)
		{
		  fread(&P[pc_sph].Ne, sizeof(float), 1, fd);
		  pc_sph++;
		}
	      SKIP;
	    }
	  else
	    for(n = 0, pc_sph = pc; n < header1.npart[0]; n++)
	      {
		P[pc_sph].Ne = 1.0;
		pc_sph++;
	      }
	}
      
#ifdef OUTPUTPOTENTIALS
      SKIP;
      for(k = 0, pc_new = pc; k < 6; k++)
	{
	  for(n = 0; n < header1.npart[k]; n++)
	    {
	      fread(&P[pc_new].Pot, sizeof(float), 1, fd);
	      pc_new++;
	    }
	}
      SKIP;
#endif

#ifdef OUTPUTTIMESTEP
      SKIP;
      for(k = 0, pc_new = pc; k < 6; k++)
	{
	  for(n = 0; n < header1.npart[k]; n++)
	    {
	      fread(&P[pc_new].Tstp, sizeof(float), 1, fd);
	      pc_new++;
	    }
	}
      SKIP;
#endif

      fclose(fd);
    }

  Time = header1.time;
  Redshift = header1.time;

  return;
}

void Snap::unload_snapshot()
{
  if(allocated)
    {
      free(P_free);
      free(Id_free);
    }
  allocated = 0;
  return;
}

int Snap::allocate_memory(void)
{
  //printf("allocating memory...\n");
  
  if(!(P = (particle_data*) malloc(NumPart * sizeof(struct particle_data))))
    {
      fprintf(stderr, "failed to allocate memory.\n");
      exit(0);
    }
  P_free = P;
  P--;				/* start with offset 1 */
  
  
  if(!(Id = (int*) malloc(NumPart * sizeof(int))))
    {
      fprintf(stderr, "failed to allocate memory.\n");
      exit(0);
    }
  Id_free = Id;
  Id--;				/* start with offset 1 */
  
  //printf("allocating memory...done\n");

  return 1;
}

Snap::Snap()
{
  allocated = 0;
  return;
}

Snap::~Snap()
{
  //free allocated memory if allocated
  if(allocated)
    {
      free(P_free);
      free(Id_free);
    }
  return;
}

int Snap::ID(int n)
{
  return Id[n+1]; //0th Id doesn't contain a particle ID, it contains '0', there are nparticle+1 elements in Id[]
}

float Snap::pos(int n, int k)
{
  return P[n+1].Pos[k];
}

float Snap::vel(int n, int k)
{
  return P[n+1].Vel[k];
}

float Snap::mass(int n)
{
  return P[n+1].Mass;
}

int Snap::type(int n)
{
  return P[n+1].Type;
}

float Snap::rho(int n)
{
  return P[n+1].Rho;
}

float Snap::u(int n)
{
  return P[n+1].U;
}

float Snap::temp(int n)
{
  return P[n+1].Temp;
}

float Snap::ne(int n)
{
  return P[n+1].Ne;
}

float Snap::pot(int n)
{
  return P[n+1].Pot;
}

float Snap::tstp(int n)
{
  return P[n+1].Tstp;
}

int Snap::np()
{
  return NumPart;
}

int Snap::ng()
{
  return Ngas;
}

double Snap::t()
{
  return Time;
}

double Snap::z()
{
  return Redshift;
}
