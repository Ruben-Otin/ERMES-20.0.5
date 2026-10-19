*##############################################################################
*################  Dirichlet elements ( PEC | PMC | TEC )  ####################
*##############################################################################
// Perfect Electric Conductor (PEC) surface elements
*set cond Perfect_Electric_Conductor_(PEC) *elems *canrepeat
*loop elems *onlyInCond
PEC = n([*GlobalNodes(1),*GlobalNodes(2),*GlobalNodes(3)]);
*end elems

// Perfect Magnetic Conductor (PMC) surface elements
*set cond Perfect_Magnetic_Conductor_(PMC) *elems *canrepeat
*loop elems *onlyInCond
PMC = n([*GlobalNodes(1),*GlobalNodes(2),*GlobalNodes(3)]);
*end elems

// Transversal ElectriC field (TEC) surface elements
*set cond Transversal_Electric_Condition_(TEC) *elems *canrepeat
*loop elems *onlyInCond
TEC = n([*GlobalNodes(1),*GlobalNodes(2),*GlobalNodes(3)]);
*end elems