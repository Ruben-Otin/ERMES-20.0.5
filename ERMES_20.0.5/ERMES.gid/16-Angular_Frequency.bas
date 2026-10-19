*##############################################################################
*##########################  Angular frequency  ###############################
*##############################################################################
*Realformat "%.18f"
*set var dospi(real) = 6.283185307179586476925286766559
*set var pfrequency(real) = GenData(Frequency,real)
ProblemFrequency = *operation(pfrequency*dospi);