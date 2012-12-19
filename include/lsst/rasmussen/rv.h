#if !defined(LSST_RASMUSSEN_RV_H)
#define LSST_RASMUSSEN_RV_H
/* rv.h		24 July 89	Roland Vanderspek

   Include file for RV tools software
*/
struct data_str {
    float data[9];
    int framenum;
    int chipnum;
    int x,y;
    char mode;
};

#endif
