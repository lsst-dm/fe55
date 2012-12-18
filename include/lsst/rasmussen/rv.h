#if !defined(LSST_RASMUSSEN_H)
#define LSST_RASMUSSEN_H
/* rv.h		24 July 89	Roland Vanderspek

   Include file for RV tools software
*/
struct data_str {
    char mode;
    int framenum;
    int chipnum;
    int x,y;
    float data[9];
};

#endif
