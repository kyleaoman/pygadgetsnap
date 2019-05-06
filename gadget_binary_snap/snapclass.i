%module gadget_binary_snap
%{
#include "snapclass.h"
%}

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

%inline %{
  extern Snap S;
%}
