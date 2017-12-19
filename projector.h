#ifndef PROJECTOR_H
#define PROJECTOR_H

void vc_pseudoprojection(pswf_t* wf_ref, pswf_t* wf_proj, int BAND_NUM, double* results);

float* pseudoprojection(pswf_t* wf_ref, pswf_t* wf_proj, int BAND_NUM);

double* read_and_project(int BAND_NUM, double* kpt_weights, char* bulkfile, char* defectfile);

#endif

