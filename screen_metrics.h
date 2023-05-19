#ifndef SCREEN_METRICS_H
#define SCREEN_METRICS_H

void get_cpx(const double cmm, const int cgpx, int* const cpx, size_t* const I, const int verb);

void get_ca_dims(const int ppc, const int gpx, const int gpy, size_t* const nrows, size_t* const ncols, size_t* const nwords, const int verb);

void SetWinNodeco(Display* const dis, Window win);

#endif // SCREEN_METRICS_H
