##
# @file
# This file is part of SeisSol.
#
# @author Alex Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
#
# @section LICENSE
# Copyright (c) 2014, SeisSol Group
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# 
# 1. Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
# 
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
# 
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from this
#    software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#
# @section DESCRIPTION
# Performance model for SeisSol setups.

computeMeshStatistics <- function( i_meshName,
                                   i_pathToCsvFile,
                                   i_performanceParameters,
                                   i_outputFileTable,
                                   i_outputFilePlot ) {
  message( 'reading mesh', i_pathToCsvFile )
  l_meshStatistics <- read.csv( file = i_pathToCsvFile,
                                header = TRUE )

  #
  # compute data to send and receive in MPI-communication
  #
  # mpi-information per time step in GB, regular mpi faces
  l_meshStatistics$mpi_send <- l_meshStatistics$mpi_faces*i_performanceParameters$dofs_size/1024^3

  # mpi-information per time step in GB, dr mpi face come on top
  l_meshStatistics$mpi_send <- l_meshStatistics$mpi_send + (l_meshStatistics$mpi_dr_faces*i_performanceParameters$dofs_size/1024^3)

  # so far everything is the same
  l_meshStatistics$mpi_receive <- l_meshStatistics$mpi_send

  #
  # compute ratios
  #
  l_meshStatistics$ratio_dynamic_rupture_faces_elements <- l_meshStatistics$dynamic_rupture_faces / l_meshStatistics$elements

  # compute timings for
  # * time integraton of all elements
  # * volume integration of all elements
  # * boundary integration of all elements
  # * additional flux computation for dunamic rupture faces
  # * time integration of MPI-elements
  # * communication of MPI-elements
  # * computation of volume integration and non-MPI time integration
  if( !i_performanceParameters$mic ) {
    l_meshStatistics$time_integration                 <- l_meshStatistics$elements*i_performanceParameters$time_integration
    l_meshStatistics$volume_integration               <- l_meshStatistics$elements*i_performanceParameters$volume_integration
  }
  
  l_meshStatistics$boundary_integration             <- l_meshStatistics$elements*i_performanceParameters$boundary_integration  
  l_meshStatistics$time_communication               <- (l_meshStatistics$mpi_send +  l_meshStatistics$mpi_receive) / i_performanceParameters$network_bandwidth  
  l_meshStatistics$dynamic_rupture_fluxes           <- l_meshStatistics$dynamic_rupture_faces*i_performanceParameters$dynamic_rupture  

  if( i_performanceParameters$mic ) {
    l_meshStatistics$time_integration_mpi             <- l_meshStatistics$mpi_elements*i_performanceParameters$time_integration
    l_meshStatistics$pci_communication                <- ( (l_meshStatistics$mpi_faces*i_performanceParameters$dofs_size*2 + l_meshStatistics$mpi_dr_faces*i_performanceParameters$dofs_size) /1024^3) / i_performanceParameters$pci_bandwidth
    l_meshStatistics$volume_time_integration_interior <- l_meshStatistics$elements*i_performanceParameters$volume_integration + l_meshStatistics$interior_elements*i_performanceParameters$time_integration 

    l_meshStatistics$overlapped_comm_time_volume      <- pmax( l_meshStatistics$time_communication + l_meshStatistics$pci_communication, l_meshStatistics$volume_time_integration_interior )
    l_meshStatistics$cpu_dynamic_rupture              <- l_meshStatistics$dynamic_rupture_fluxes - pmax( l_meshStatistics$volume_time_integration_interior -  l_meshStatistics$overlapped_comm_time_volume, 0 )
    l_meshStatistics$overlapped_dyn_rupture_fluxes    <- pmax( l_meshStatistics$cpu_dynamic_rupture + (l_meshStatistics$dynamic_rupture_faces*i_performanceParameters$dofs_size/1024^3) / i_performanceParameters$pci_bandwidth,
                                                               l_meshStatistics$boundary_integration )
  }

#new overlap
#  if( i_performanceParameters$mic ) {
#    l_meshStatistics$time_integration_mpi             <- l_meshStatistics$mpi_elements*i_performanceParameters$time_integration
#    l_meshStatistics$pci_communication                <- ( (l_meshStatistics$mpi_faces*i_performanceParameters$dofs_size*2 + l_meshStatistics$mpi_dr_faces*i_performanceParameters$dofs_size) /1024^3) / i_performanceParameters$pci_bandwidth
#    l_meshStatistics$volume_time_integration_interior <- l_meshStatistics$elements*i_performanceParameters$volume_integration + l_meshStatistics$interior_elements*i_performanceParameters$time_integration + l_meshStatistics$interior_elements*i_performanceParameters$boundary_integration
#    
#    l_meshStatistics$overlapped_comm_time_volume      <- pmax( l_meshStatistics$time_communication + l_meshStatistics$pci_communication, l_meshStatistics$volume_time_integration_interior )
#    l_meshStatistics$cpu_dynamic_rupture              <- l_meshStatistics$dynamic_rupture_fluxes - pmax( l_meshStatistics$volume_time_integration_interior -  l_meshStatistics$overlapped_comm_time_volume, 0 )
#    l_meshStatistics$overlapped_dyn_rupture_fluxes    <- pmax( l_meshStatistics$cpu_dynamic_rupture + (l_meshStatistics$dynamic_rupture_faces*i_performanceParameters$dofs_size/1024^3) / i_performanceParameters$pci_bandwidth,
#                                                               l_meshStatistics$mpi_elements*i_performanceParameters$boundary_integration )
#  }
  
  if( !i_performanceParameters$mic ) {
    l_meshStatistics$time_total                       <- l_meshStatistics$time_integration + l_meshStatistics$volume_integration + l_meshStatistics$boundary_integration + l_meshStatistics$dynamic_rupture_fluxes + l_meshStatistics$time_communication
  }
  else {
    l_meshStatistics$time_total                       <- l_meshStatistics$time_integration_mpi + l_meshStatistics$overlapped_comm_time_volume + l_meshStatistics$overlapped_dyn_rupture_fluxes
  }
  message( 'summary mesh statistics')
  print( summary( l_meshStatistics ) )

  message( 'writing csv' )
  write.csv( x=l_meshStatistics, file=i_outputFileTable)
  
  message( 'plotting information' )
  pdf( file = i_outputFilePlot,
       paper = 'a4r',
       width=20,
       height=10)

  # print title page
  plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
  text(5, 8, i_meshName)
  text(5, 7, Sys.time())

  # overview boxplots
  if( !i_performanceParameters$mic ) {
    par( mfrow=c(1,6) )
    boxplot(l_meshStatistics$time_integration,       main='time integration'     )
    boxplot(l_meshStatistics$volume_integration,     main='volume integration'   )
    boxplot(l_meshStatistics$boundary_integration,   main='boundary integration' )
    boxplot(l_meshStatistics$dynamic_rupture_fluxes, main='dyn. rupt. fluxes'    )
    boxplot(l_meshStatistics$time_communication,     main='communication'        )
    boxplot(l_meshStatistics$time_total,             main='total time'           )
  }
  else {
    par( mfrow=c(1,6) )
    boxplot( l_meshStatistics$time_integration_mpi,             main='time integration of copy layer'                   )
    boxplot( l_meshStatistics$volume_time_integration_interior, main='time & volume int.'                               )
    boxplot( l_meshStatistics$overlapped_comm_time_volume,      main='overlapped comm., int. time and vol. integration' )
    boxplot( l_meshStatistics$boundary_integration,             main='boundary integration'                             )
    boxplot( l_meshStatistics$overlapped_dyn_rupture_fluxes,    main='overlapped dyn. rupt. and bnd. intgration'        )
    boxplot( l_meshStatistics$time_total,                       main='total time'                                       )
  }
  
  # detailed plot of every characteristic value
  #for( l_name in names(l_meshStatistics[-1]) ) {
  #  layout( matrix(c(1,2), 2, 2, byrow = TRUE), 
  #          widths=c(4,1) )
  #  plot( x=l_meshStatistics$partition,
  #        y=l_meshStatistics[,l_name],
  #        xlab="partition",
  #        ylab=l_name )
  #  boxplot(l_meshStatistics[l_name])
  #}
  dev.off()
  
  return( median(l_meshStatistics$time_total) )
}

l_config = list ( mesh_names = list( 'statistics_cube320_400_640_1024',
                                     'statistics_cube400_640_640_2048',
                                     'statistics_cube640_640_800_4096',
                                     'statistics_cube640_900_1280_9216',
                                     'statistics_Landers191M_02_05_02_512',
                                     'statistics_Landers191M_02_05_02_768',
                                     'statistics_Landers191M_02_05_02_1024',
                                     'statistics_Landers191M_02_05_02_1536',
                                     'statistics_Landers191M_02_05_02_3072',
                                     'statistics_Landers191M_02_05_02_2048',
                                     'statistics_Landers191M_02_05_02_4096',
                                     'statistics_Landers191M_02_05_02_6144',
                                     'statistics_Landers191M_02_05_02_9216',
                                     'statistics_Landers191M_02_05_02_12288',
                                     'statistics_Landers191M_02_05_02_24756'),
                  
                  statistics_directory = 'mesh_statistics_gb',
                  
                  performance_parameters = list( supermuc = list( basisFunctions       = 56,
                                                                  quantities           = 9,
                                                                  dofs_size            = 56*9*8,
                                                                  network_bandwidth    = 0.8,                              # TODO: missing reference point
                                                                  mic                  = FALSE,
                                                                  boundary_integration = 389.68 * 1024 / 191098540 / 1000, # reference, landers_1024
                                                                  time_integration     = 169.04 * 1024 / 191098540 / 1000, # reference, landers_1024
                                                                  volume_integration   = 126.88 * 1024 / 191098540 / 1000, # reference, landers_1024
                                                                  #boundary_integration = 43.26 * 9216 / 191098540 / 1000,  # reference, landers_9216
                                                                  #time_integration     = 18.86 * 9216 / 191098540 / 1000,  # reference, landers_9216
                                                                  #volume_integration   = 14.06 * 9216 / 191098540 / 1000,  # reference, landers_9216
                                                                  #time_integration     = 36.33 / 400000 / 100,             # reference, cubes_9216
                                                                  #volume_integration   = 27.18 / 400000 / 100,             # reference, cubes_9216
                                                                  #boundary_integration = 84.86 / 400000 / 100,             # reference, cubes_9216
                                                                  dynamic_rupture      = 109.33 /  22798 / 1000 ),          # TODO: missing reference point
                                                 stampede = list( basisFunctions       = 56,
                                                                  quantities           = 9,
                                                                  dofs_size            = 56*9*8,
                                                                  network_bandwidth    = 0.8,
                                                                  mic                  = TRUE,
                                                                  pci_bandwidth        = 6.3,
                                                                  boundary_integration = 230.68 * 2 / 386518 / 1000 * 56 / 60, # extrapolated from LOH1_2
                                                                  time_integration     = 131.84 * 2 / 386518 / 1000 * 56 / 60, # extrapolated frpm LOH1_2
                                                                  volume_integration   =  88.51 * 2 / 386518 / 1000 * 56 / 60, # extrapolated from LOH1_2
                                                                  dynamic_rupture      = 109.33 /  22798 / 1000 ), # same as SuperMUC, but not relevant here; no impact on TH-2
                                                 tianhe   = list( basisFunctions       = 56,
                                                                  quantities           = 9,
                                                                  dofs_size            = 56*9*8,
                                                                  #network_bandwidth    = 0.367, # fitted via cube-runtimes of ~100.25 seconds
                                                                  #network_bandwidth    = 0.61, # fitted via landers_2048, 3 cards
                                                                  network_bandwidth    = 0.525, #fitted via landers_2048, rank 2474: 4582 mpi-facesm 65.56s comm. time
                                                                  mic                  = TRUE,
                                                                  pci_bandwidth        = 3.2,
                                                                  boundary_integration = 230.68 * 2 / 386518 / 1000,  # extrapolated from LOH1_2
                                                                  time_integration     = 131.84 * 2 / 386518 / 1000,  # extrapolated from LOH1_2
                                                                  volume_integration   =  88.51 * 2 / 386518 / 1000,  # extrapolated from LOH1_2
                                                                  dynamic_rupture      =  26.79     /   2800 / 1000 ) # rank 2931 for 6144 ranks on TH-2
                                                )
                 )

l_expectedSimulationTime = list()

message( 'lets go' )
for( l_machine in list('supermuc', 'stampede', 'tianhe') ) {
  for( l_meshName in l_config$mesh_names ) {
    message( 'analyzing: ', l_meshName )
    l_expectedSimulationTime[[paste(l_meshName,'_',l_machine, sep='')]] =
     computeMeshStatistics( i_meshName              = l_meshName,
                            i_pathToCsvFile         = paste(l_config$statistics_directory,'/',l_meshName,'.csv', sep=''),
                            i_performanceParameters = (l_config$performance_parameters)[[l_machine]],
                            i_outputFileTable       = paste(l_meshName,'_',l_machine,'.csv', sep=''),
                            i_outputFilePlot        = paste(l_meshName,'_',l_machine,'.pdf', sep='') )
  }
}

l_expectedSimulationTime$statistics_Landers191M_02_05_02_768_tianhe = 0.581

for( l_machine in list('supermuc', 'stampede', 'tianhe') ) {
  l_strongScalingNodes = list(768, 1024, 1536, 2048, 3072, 4096, 6144, 9216, 12288, 24756 )
  l_runtime = list()
  for( l_nodes in l_strongScalingNodes ) {
    l_runtime = c( l_runtime, l_expectedSimulationTime[[paste('statistics_Landers191M_02_05_02_', l_nodes, '_', l_machine, sep='')]] * ( l_nodes / l_strongScalingNodes[[1]]) )
  }
  l_runtime <- l_runtime[[1]] / unlist(l_runtime)
  
  print(l_machine)
  print(unlist(l_strongScalingNodes))
  print(l_runtime)
  
  plot( x=l_strongScalingNodes, l_runtime, ylab='parallel efficiency', main=l_machine )
}
