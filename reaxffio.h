#ifndef _REAXFFIO_H_
#define _REAXFFIO_H_
double write_reaxff_input_files(int obj, int ip);
double write_para_reax();
double write_para_everystep(FILE *fp);
double analysis_for_reaxff_output(char filename[]);
double Initiation4reaxff();
double ReadParas4reaxff();
double ReadMaxdp4reaxff();
double ReadUplimit4reaxff();
double ReadDownlimit4reaxff();
double write_data_md(int obj, int ip);
#endif
