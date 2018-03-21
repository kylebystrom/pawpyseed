#ifndef PSEUDOPROJECTOR_H
#define PSEUDOPROJECTOR_H
#include "utils.h"

void vc_pseudoprojection(pswf_t* wf_ref, pswf_t* wf_proj, int BAND_NUM, double* results);

double* pseudoprojection(pswf_t* wf_ref, pswf_t* wf_proj, int BAND_NUM);

#endif