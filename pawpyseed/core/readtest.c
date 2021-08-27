#include "reader.h"
#include "utils.h"
#include <stdio.h>
#include <stdlib.h>

int main() {
  double w[16];
  for (int i = 0; i < 16; i++)
    w[i] = 0.0625;
  read_wavefunctions("charge_0/WAVECAR", w);
  return 0;
}
