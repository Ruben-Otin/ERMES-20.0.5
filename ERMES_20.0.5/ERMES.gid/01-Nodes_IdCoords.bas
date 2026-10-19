*##############################################################################
*#####################  Nodes list ( ID and coordinates)  #####################
*##############################################################################
// List of nodes (Ids and XYZ coordinates)
*Realformat "%.18f"
*if(strcmp(GenData(Length_units),"custom")==0)
*set var pLengthFactor(real) = GenData(custom_L_factor,real)
*elseif(strcmp(GenData(Length_units),"m")==0)
*set var pLengthFactor(real) = 1
*elseif(strcmp(GenData(Length_units),"dm")==0)
*set var pLengthFactor(real) = 0.1
*elseif(strcmp(GenData(Length_units),"cm")==0)
*set var pLengthFactor(real) = 0.01
*elseif(strcmp(GenData(Length_units),"mm")==0)
*set var pLengthFactor(real) = 0.001
*endif
*loop nodes
No[*nodesnum] = p(*operation(pLengthFactor*NodesCoord(1)),*operation(pLengthFactor*NodesCoord(2)),*operation(pLengthFactor*NodesCoord(3)));
*end