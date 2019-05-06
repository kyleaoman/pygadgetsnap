#ifndef INC_SNAPCLASS
#define INC_SNAPCLASS

//c includes defined in gadget read example
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

struct io_header_1
{
  int npart[6];
  double mass[6];
  double time;
  double redshift;
  int flag_sfr;
  int flag_feedback;
  int npartTotal[6];
  int flag_cooling;
  int num_files;
  double BoxSize;
  double Omega0;
  double OmegaLambda;
  double HubbleParam;
  char fill[256 - 6 * 4 - 6 * 8 - 2 * 8 - 2 * 4 - 6 * 4 - 2 * 4 - 4 * 8];	/* fills to 256 Bytes */
};

struct particle_data
{
  float Pos[3];
  float Vel[3];
  float Mass;
  int Type;

  float Rho, U, Temp, Ne;
  float Pot;
  float Tstp;
};

class Snap
{
 private:
  struct io_header_1 header1;
  int NumPart, Ngas;
  struct particle_data *P, *P_free;
  int *Id, *Id_free;
  double Time, Redshift;
  int allocate_memory(void);
  int allocated;

 public:
  Snap();
  ~Snap();
  void load_snapshot(char*, int);
  void unload_snapshot();

  int ID(int);
  float pos(int, int);
  float vel(int, int);
  float mass(int);
  int type(int);
  float rho(int);
  float u(int);
  float temp(int);
  float ne(int);
  float pot(int);
  float tstp(int);
  int np();
  int ng();
  double t();
  double z();
};

#endif //INC_SNAPCLASS
