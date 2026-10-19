*##############################################################################
*######################  Solve and print results ##############################
*##############################################################################
*Intformat "%i"
*Realformat "%.18f"
*#######################  Electrostatic solver  ###############################
*if(strcmp(GenData(Problem_mode),"Electrostatic")==0)
// Electrostatic solver
*if(strcmp(GenData(Max_iterations),"custom")==0)
*set var maxit(int)=GenData(custom_Max_it,real)
*else
*set var maxit(int)=GenData(Max_iterations,real)
*endif
*if(strcmp(GenData(Tolerance),"custom")==0)
*set var toler(real)=GenData(custom_Tol,real)
*else
*set var toler(real)=GenData(Tolerance,real)
*endif
*if(strcmp(GenData(Iteration_step),"custom")==0)
*set var stepit(int)=GenData(custom_It_step,real)
*else
*set var stepit(int)=GenData(Iteration_step,real)
*endif
LinearSolver *GenData(Preconditioner) = *GenData(Solver)(*maxit,*stepit,*toler);

// Solve electrostatic problem
ElectromagneticStrategy.Solve(Electromagnetic_Group);

// Print mesh
ProblemType = Print_High_Order_Mesh;

// Calculate gradients
ProblemType = Calculate_Static_E;

// Print results
*if(strcmp(GenData(Electrostaic_potential),"1")==0)
PrintOnNodes(V);
*endif
*if(strcmp(GenData(Electrostatic_E_field),"1")==0) 
ProblemType = Show_Static_E;
*endif
*if(strcmp(GenData(Electrostatic_current_density),"1")==0)
ProblemType = Show_Static_J;
*endif
*if(strcmp(GenData(Electrostatic_Joule_heating),"1")==0)
ProblemType = Show_Static_Joule_Heating;
*endif

*if(strcmp(GenData(Export_electrostatic_currents),"1")==0)
// Export electrostatic currents to file
ProblemType = Export_Static_Currents;

*endif
// Read projection elements
ProblemType = Read_Projector_Elements;

// Print static field projections
ProblemType = Print_Static_Projections;
*else
*##########################  Complex solver  ##################################
// Complex solver
*if(strcmp(GenData(Max_iterations),"custom")==0)
*set var maxit(int)=GenData(custom_Max_it,real)
*else
*set var maxit(int)=GenData(Max_iterations,real)
*endif
*if(strcmp(GenData(Tolerance),"custom")==0)
*set var toler(real)=GenData(custom_Tol,real)
*else
*set var toler(real)=GenData(Tolerance,real)
*endif
*if(strcmp(GenData(Iteration_step),"custom")==0)
*set var stepit(int)=GenData(custom_It_step,real)
*else
*set var stepit(int)=GenData(Iteration_step,real)
*endif
LinearSolver *GenData(Preconditioner) = *GenData(Solver)(*maxit,*stepit,*toler);

// Solve complex problem
ElectromagneticStrategy.Solve(Electromagnetic_Group);

// Calculate fields
ProblemType = Calculate_E;
ProblemType = Calculate_H;
ProblemType = Calculate_B;
ProblemType = Calculate_J;

// Read projection elements
ProblemType = Read_Projector_Elements;

*if(strcmp(GenData(Sweep_frequency_mode),"0")==0)
// Print field projections
ProblemType = Print_Field_Projections;

// Print mesh
ProblemType = Print_High_Order_Mesh;

*if( (strcmp(GenData(E),"1")==0) || (strcmp(GenData(E(t)),"1")==0) )
// Print E field
Print(IMAG_E);
Print(REAL_E);
Calculate(MOD_E);
Print(MOD_E);

*endif
*if( (strcmp(GenData(H),"1")==0) || (strcmp(GenData(H(t)),"1")==0) )
// Print H field
Print(IMAG_H);
Print(REAL_H);
Calculate(MOD_H);
Print(MOD_H);

*endif
*if( (strcmp(GenData(J),"1")==0) || (strcmp(GenData(J(t)),"1")==0) )
// Print J currents
Print(REAL_J);
Print(IMAG_J);
Calculate(MOD_J);
Print(MOD_J);

*endif
*if( (strcmp(GenData(B),"1")==0) || (strcmp(GenData(B(t)),"1")==0) )
// Print B field
Print(IMAG_B);
Print(REAL_B);
Calculate(MOD_B);
Print(MOD_B);

*endif
// Other results
*if( (strcmp(GenData(Complex_frequency_mode),"0")==0) || (strcmp(GenData(Problem_mode),"Cold_plasma")==0) )
*if(strcmp(GenData(Poynting_vector),"1")==0)
ProblemType = Show_Poynting_Vector;
*endif
*if( (strcmp(GenData(Problem),"Full_wave")==0) && (strcmp(GenData(Joule_heating),"1")==0) )
ProblemType = Show_Joule_Heating;
*endif
*if( (strcmp(GenData(Problem),"Full_wave")==0) && (strcmp(GenData(Lorentz_force),"1")==0) )
ProblemType = Show_Lorentz_Force;
*endif
*if( (strcmp(GenData(Problem),"Cold_plasma")==0) && (strcmp(GenData(Electron_density),"1")==0) )
ProblemType = Show_Electron_Density;
*endif
*if( (strcmp(GenData(Problem),"Cold_plasma")==0) && (strcmp(GenData(External_B_field),"1")==0) )
ProblemType = Show_B_External;
*endif
*if( (strcmp(GenData(Problem),"Cold_plasma")==0) && (strcmp(GenData(Permittivity_tensor),"1")==0) )
ProblemType = Show_Permittivity_Tensor;
*endif
*if( (strcmp(GenData(Problem),"Cold_plasma")==0) && (strcmp(GenData(E_parallel),"1")==0) )
ProblemType = Show_E_Parallel;
*endif
*if( (strcmp(GenData(Problem),"Cold_plasma")==0) && (strcmp(GenData(E_perpendicular),"1")==0) )
ProblemType = Show_E_Perpendicular;
*endif
*endif

// Time domain results
*if( (strcmp(GenData(Complex_frequency_mode),"0")==0) || (strcmp(GenData(Problem_mode),"Cold_plasma")==0) )
*if( (strcmp(GenData(E(t)),"1")==0) || (strcmp(GenData(H(t)),"1")==0) || (strcmp(GenData(B(t)),"1")==0) || (strcmp(GenData(J(t)),"1")==0) || ( (strcmp(GenData(Problem),"Cold_plasma")==0) && (strcmp(GenData(E_parallel_(t)),"1")==0) ) || ( (strcmp(GenData(Problem),"Cold_plasma")==0) && (strcmp(GenData(E_perpendicular_(t)),"1")==0) ) )
*set var tfrequency(real) = GenData(Frequency,real)
*set var timeperiod(real) = operation(1/tfrequency)
*set var timeinc(real)    = operation(timeperiod/32)
*set var initime(real)    = timeinc(real)
*set var endtime(real)    = operation(timeperiod)
*format "%18.15f%18.15f%18.15f%18.15f"
for(TIME = *initime; TIME < *endtime + *timeinc; TIME = TIME + *timeinc)
{
    CreateTimeStep();
*if(strcmp(GenData(E(t)),"1")==0)
    Calculate(E);
    Print(E);
*endif
*if(strcmp(GenData(B(t)),"1")==0)
    Calculate(B);
    Print(B);
*endif
*if(strcmp(GenData(H(t)),"1")==0)
    Calculate(H);
    Print(H);
*endif
*if(strcmp(GenData(J(t)),"1")==0)
    Calculate(J);
    Print(J);
*endif
*if( (strcmp(GenData(Problem),"Cold_plasma")==0) && (strcmp(GenData(E_parallel_(t)),"1")==0) )
    ProblemType = Show_Eparallel_t;
*endif
*if( (strcmp(GenData(Problem),"Cold_plasma")==0) && (strcmp(GenData(E_perpendicular_(t)),"1")==0) )
    ProblemType = Show_Eperpendicular_t;
*endif
}
*endif
*endif

*endif
*endif 