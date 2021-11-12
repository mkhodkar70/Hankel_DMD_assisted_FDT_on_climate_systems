# Hankel_DMD_assisted_FDT_on_climate_systems
The code applies the fluctuation-dissipation theorem to a dry-core general circulation model, whose basis functions for model reduction can be EOFs, DMD modes or Hankel-DMD modes (decided by the user), to predict the responses of the idealized climate system to different thermal and mechanical forcings. 

Combined_FDT_HDMD_Code.m: It's the main code that has to be run to plot the zonal-wind and temperature responses of the system to an external forcing that is set in the forcing segment of the code. The method of mode reduction, the sampling rate and the hyperparameters of the FDT, DMD and EOF (POD) are all determined at the top of the code. 
dmd_sorted: It sorts the DMD moes based on the decay rate of their corresponding eigenvalues. 
my_compress: It converts the multi-block delay-embedded DMD mode to a single-block mode by taking the last block as the basis function (projection mode). It can also be modified (by uncommenting the My way segment) to calculate the basis function as an average of all blocks. The results will change only insubstantially. The number of blocks equals the delay-embedding dimension.
b2r.m: It creates a blue-to-red color map. 
