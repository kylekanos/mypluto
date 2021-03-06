#define  PHYSICS                 HD
#define  DIMENSIONS              2
#define  COMPONENTS              2
#define  GEOMETRY                CARTESIAN
#define  BODY_FORCE              NO
#define  COOLING                 NO
#define  INTERPOLATION           LINEAR
#define  TIME_STEPPING           RK2
#define  DIMENSIONAL_SPLITTING   NO
#define  NTRACER                 0
#define  USER_DEF_PARAMETERS     6

/* -- physics dependent declarations -- */

#define    EOS                     IDEAL
#define    ENTROPY_SWITCH          NO
#define    THERMAL_CONDUCTION      SUPER_TIME_STEPPING
#define    VISCOSITY               NO
#define    ROTATING_FRAME          NO

/* -- pointers to user-def parameters -- */

#define  BMAG               0
#define  THETA              1
#define  T_IN               2
#define  T_OUT              3
#define  RHO_IN             4
#define  RHO_OUT            5

/* -- supplementary constants (user editable) -- */ 

#define  INITIAL_SMOOTHING     NO
#define  WARNING_MESSAGES      YES
#define  PRINT_TO_FILE         YES
#define  INTERNAL_BOUNDARY     NO
#define  SHOCK_FLATTENING      NO
#define  ARTIFICIAL_VISCOSITY  NO
#define  CHAR_LIMITING         NO
#define  LIMITER               VANLEER_LIM
#define  STS_nu                0.01
