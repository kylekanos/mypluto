#define  PHYSICS                 HD
#define  DIMENSIONS              1
#define  COMPONENTS              1
#define  GEOMETRY                CARTESIAN
#define  BODY_FORCE              NO
#define  COOLING                 NO
#define  INTERPOLATION           LINEAR
#define  TIME_STEPPING           RK2
#define  DIMENSIONAL_SPLITTING   YES
#define  NTRACER                 0
#define  USER_DEF_PARAMETERS     1

/* -- physics dependent declarations -- */

#define    EOS                     IDEAL
#define    ENTROPY_SWITCH          NO
#define    THERMAL_CONDUCTION      SUPER_TIME_STEPPING
#define    VISCOSITY               NO
#define    ROTATING_FRAME          NO

/* -- pointers to user-def parameters -- */

#define  SCRH               0

/* -- supplementary constants (user editable) -- */ 

#define  INITIAL_SMOOTHING     NO
#define  WARNING_MESSAGES      NO
#define  PRINT_TO_FILE         NO
#define  INTERNAL_BOUNDARY     YES
#define  SHOCK_FLATTENING      NO
#define  ARTIFICIAL_VISCOSITY  NO
#define  CHAR_LIMITING         NO
#define  LIMITER               DEFAULT
#define  STS_nu                0.01
