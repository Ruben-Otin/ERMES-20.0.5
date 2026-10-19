
#if !defined(KRATOS_KRATOS_LINEAR_SOLVERS_H_INCLUDED)
#define  KRATOS_KRATOS_LINEAR_SOLVERS_H_INCLUDED

// External includes 
#include <memory>

// Project includes
#include "../includes/define.h"

// Complex solvers
#include "bicg_complex_solver.h"
#include "qmr_complex_solver.h"

// Complex solvers for cold plasma
#include "qmr_cp_complex_solver.h"
#include "bicg_cp_complex_solver.h"

// Internal direct solver ( sparse LU, all storage modes )
#include "direct_complex_solver.h"

#endif //  KRATOS_KRATOS_LINEAR_SOLVERS_H_INCLUDED defined 


