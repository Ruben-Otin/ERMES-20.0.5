*##############################################################################
*####################  List of nodes near a singularity  ######################
*##############################################################################
// List of singular nodes (ungaged layers formulation)
*Set cond Point_Singularity *nodes
*Add cond Line_Singularity *nodes
*Add cond Surface_Singularity *nodes
*loop nodes *OnlyInCond
No[*NodesNum].Sg.Fix(*cond(Ungaged_layers));
*end nodes