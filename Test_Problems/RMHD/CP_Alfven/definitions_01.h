#define  PHYSICS                 RMHD
#define  DIMENSIONS              2
#define  COMPONENTS              3
#define  GEOMETRY                CARTESIAN
#define  BODY_FORCE              NO
#define  COOLING                 NO
#define  INTERPOLATION           LINEAR
#define  TIME_STEPPING           RK2
#define  DIMENSIONAL_SPLITTING   NO
#define  NTRACER                 0
#define  USER_DEF_PARAMETERS     1

/* -- physics dependent declarations -- */

#define    EOS                     IDEAL
#define    ENTROPY_SWITCH          NO
#define    MHD_FORMULATION         CONSTRAINED_TRANSPORT

/* -- pointers to user-def parameters -- */

#define  SCRH               0

/* -- supplementary constants (user editable) -- */ 

#define  INITIAL_SMOOTHING     NO
#define  WARNING_MESSAGES      YES
#define  PRINT_TO_FILE         YES
#define  INTERNAL_BOUNDARY     NO
#define  SHOCK_FLATTENING      NO
#define  ARTIFICIAL_VISCOSITY  NO
#define  CHAR_LIMITING         NO
#define  LIMITER               MC_LIM
#define  CT_EMF_AVERAGE        ARITHMETIC
#define  CT_EN_CORRECTION      NO
#define  ASSIGN_VECTOR_POTENTIAL  YES
#define  UPDATE_VECTOR_POTENTIAL  NO
