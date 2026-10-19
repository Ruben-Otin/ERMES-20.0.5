*##############################################################################
*##########################  Voltages on nodes  ###############################
*##############################################################################
// Voltages applied on nodes
*Realformat "%.18f"
*if(strcmp(GenData(Problem_mode),"Electrostatic")==0)
*Set cond Point_Voltage_Electrostatic *nodes
*Add cond Line_Voltage_Electrostatic *nodes
*Add cond Surface_Voltage_Electrostatic *nodes
*loop nodes *OnlyInCond
No[*NodesNum].V.Fix(*cond(Voltage,real));
*end nodes
*else
*Set cond Point_Voltage_Full_wave *nodes
*Add cond Line_Voltage_Full_wave *nodes
*Add cond Surface_Voltage_Full_wave *nodes
*loop nodes *OnlyInCond
No[*NodesNum].V.FixC([*cond(Module_Voltage,real),*cond(Phase_Voltage,real)]);
*end nodes
*endif