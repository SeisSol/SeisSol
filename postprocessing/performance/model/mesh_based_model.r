##
# Copyright (c) 2015, Intel Corporation
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# 
#     * Redistributions of source code must retain the above copyright notice,
#       this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of Intel Corporation nor the names of its contributors
#       may be used to endorse or promote products derived from this software
#       without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
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
  l_meshStatistics$mpi_send <- (l_meshStatistics$mpi_copy_elements*i_performanceParameters$dofs_size)/1024^3
  # mpi-information per time step in GB, dr mpi face come on top
  #TODO here need a more accurate number, DR is really complex depending on the neighour realitions: 
  #TODO DR-DR coupling, DR-WP coupling, for now let's overestimate (2* for time derivatives instead of time integrated)
  l_meshStatistics$mpi_send <- l_meshStatistics$mpi_send + ((2*l_meshStatistics$mpi_dr_faces*i_performanceParameters$dofs_size)/1024^3)
  # now get time
  l_meshStatistics$mpi_send <- l_meshStatistics$mpi_send/i_performanceParameters$network_bandwidth

  # receiving is asymmetric to sending :-/
  l_meshStatistics$mpi_recv <- (l_meshStatistics$mpi_ghost_elements*i_performanceParameters$dofs_size)/1024^3
  l_meshStatistics$mpi_recv <- l_meshStatistics$mpi_recv + ((2*l_meshStatistics$mpi_dr_faces*i_performanceParameters$dofs_size)/1024^3)
  l_meshStatistics$mpi_recv <- l_meshStatistics$mpi_recv/i_performanceParameters$network_bandwidth

  # compute timings for
  # * local integration copy layer
  # * local integration interior
  # * neigh integration copy layer
  # * neigh integration interior
  # * additional flux computation for dunamic rupture faces
  # TODO: we might want to substract outflow boundary conditions
  # Please note we incease the time by 10% for the copy layers since they are very small 
  # and therefore run at lower performance
  l_meshStatistics$local_integration_copy           <- l_meshStatistics$mpi_copy_elements*i_performanceParameters$local_integration*1.1
  l_meshStatistics$local_integration_interior       <- l_meshStatistics$interior_elements*i_performanceParameters$local_integration
  l_meshStatistics$neigh_integration_copy           <- l_meshStatistics$mpi_copy_elements*i_performanceParameters$neigh_integration*1.1
  l_meshStatistics$neigh_integration_interior       <- l_meshStatistics$interior_elements*i_performanceParameters$neigh_integration
  l_meshStatistics$dynamic_rupture_fluxes           <- l_meshStatistics$dynamic_rupture_faces*i_performanceParameters$dynamic_rupture  
  
  # compute MPI wait time (recv)
  l_meshStatistics$wait_recv                        <- pmax(l_meshStatistics$mpi_recv -
                                                            l_meshStatistics$local_integration_copy -
                                                            l_meshStatistics$local_integration_interior -
                                                            l_meshStatistics$neigh_integration_interior , 0)

  # compute MPI wait time (send)
  # this assumes that waitforCopyLayerSends() is posted before receiveGhostLayer()
  # this is not true, but the error intorudced by this is very minor!
  l_meshStatistics$wait_send                        <- pmax(l_meshStatistics$mpi_send -
                                                            l_meshStatistics$neigh_integration_copy -
                                                            l_meshStatistics$local_integration_interior -
                                                            l_meshStatistics$neigh_integration_interior -
                                                            l_meshStatistics$dynamic_rupture_fluxes -
                                                            l_meshStatistics$wait_recv , 0)

  # let's sum the waits up
  l_meshStatistics$time_communication_wait          <- l_meshStatistics$wait_recv +
                                                       l_meshStatistics$wait_send

  #let's sum everyting up, and we are done
  l_meshStatistics$time_total                       <- l_meshStatistics$local_integration_copy +
                                                       l_meshStatistics$local_integration_interior +
                                                       l_meshStatistics$neigh_integration_copy +
                                                       l_meshStatistics$neigh_integration_interior +
                                                       l_meshStatistics$dynamic_rupture_fluxes +
                                                       l_meshStatistics$wait_recv +
                                                       l_meshStatistics$wait_send
  
  #message( 'summary mesh statistics')
  #print( summary( l_meshStatistics ) )

  message( 'writing csv' )
  write.csv( x=l_meshStatistics, file=i_outputFileTable)
  
  message( 'plotting information' )
  pdf( file = i_outputFilePlot,
       paper = 'a4r',
       width=20,
       height=10)

  # overview boxplots
  par( mfrow=c(1,9) )
  boxplot(l_meshStatistics$local_integration_copy,         main='local copy'      )
  boxplot(l_meshStatistics$local_integration_interior,     main='local interior'  )
  boxplot(l_meshStatistics$neigh_integration_interior,     main='neigh interior'  )
  boxplot(l_meshStatistics$neigh_integration_copy,         main='neigh copy'      )
  boxplot(l_meshStatistics$dynamic_rupture_fluxes,         main='dyn. rupt.'      )
  boxplot(l_meshStatistics$mpi_recv,                       main='recv time'       )
  boxplot(l_meshStatistics$mpi_recv,                       main='send time'       )
  boxplot(l_meshStatistics$time_communication_wait,        main='wait time'       )
  boxplot(l_meshStatistics$time_total,                     main='total time'      )
  
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
  
  return( max(l_meshStatistics$time_total) )
}

printWeakScaling <- function( i_machineList,
                                   i_nodeCounts,
                                   i_meshID ) {
  for( l_machine in i_machineList ) {
    l_runtime = list()
    for( l_nodes in i_nodeCounts ) {
      l_runtime = c( l_runtime, l_expectedSimulationTime[[paste(i_meshID, '_', l_nodes, '_', l_machine, sep='')]] )
    }
    l_runtime <- l_runtime[[1]] / unlist(l_runtime)
  
    print(l_machine)
    print(unlist(i_nodeCounts))
    print(l_runtime)
  
    pdf(paste('scaling_', i_meshID, '_', l_machine,'.pdf', sep=''))
    plot( x=i_nodeCounts, l_runtime, ylab='parallel efficiency (weak-scaling)', main=l_machine )
    dev.off()
  }
}

printStrongScaling <- function( i_machineList,
                                   i_nodeCounts,
                                   i_meshID ) {
  for( l_machine in i_machineList ) {
    l_runtime = list()
    for( l_nodes in i_nodeCounts ) {
      l_runtime = c( l_runtime, l_expectedSimulationTime[[paste(i_meshID, '_', l_nodes, '_', l_machine, sep='')]] * ( l_nodes / i_nodeCounts[[1]]) )
    }
    l_runtime <- l_runtime[[1]] / unlist(l_runtime)
  
    print(l_machine)
    print(unlist(i_nodeCounts))
    print(l_runtime)
  
    pdf(paste('scaling_', i_meshID, '_', l_machine,'.pdf', sep=''))
    plot( x=i_nodeCounts, l_runtime, ylab='parallel efficiency (strong-scaling)', main=l_machine )
    dev.off()
  }
}

l_config = list ( mesh_names = list( #'statistics_cubes_1024',
                                     #'statistics_cubes_2048',
                                     #'statistics_cubes_4096',
                                     #'statistics_cubes_9216',
                                     #'statistics_Landers191M_02_05_02_512',
                                     #'statistics_Landers191M_02_05_02_768',
                                     #'statistics_Landers191M_02_05_02_1024',
                                     #'statistics_Landers191M_02_05_02_1536',
                                     #'statistics_Landers191M_02_05_02_3072',
                                     #'statistics_Landers191M_02_05_02_2048',
                                     #'statistics_Landers191M_02_05_02_4096',
                                     #'statistics_Landers191M_02_05_02_6144',
                                     #'statistics_Landers191M_02_05_02_9216',
                                     #'statistics_Landers191M_02_05_02_12288',
                                     #'statistics_Landers191M_02_05_02_24756',
                                     'statistics_LOH1_small_1',
                                     'statistics_LOH1_small_2',
                                     'statistics_LOH1_small_4',
                                     'statistics_LOH1_small_8',
                                     'statistics_LOH1_small_16',
                                     'statistics_LOH1_small_32',
                                     'statistics_LOH1_small_64',
                                     'statistics_LOH1_small_128',
                                     'statistics_LOH1_small_256',
                                     'statistics_LOH1_small_512',
                                     'statistics_LOH1_small_1024'),
                  
                  statistics_directory = 'mesh_statistics',
                  
                  performance_parameters = list( supermuc = list( basisFunctions       = 56,
                                                                  quantities           = 9,
                                                                  dofs_size            = 56*9*8,
                                                                  network_bandwidth    = 0.4,                   # unidirectionl, TODO: missing reference point
                                                                  neigh_integration = 0.00000131,               # 1250 elements per core, seissol proxy, random
                                                                  local_integration   = 0.00000244,             # 1250 elements per core, seissol proxy, random
                                                                  dynamic_rupture      = 109.33 /  22798 / 1000 )   # TODO: missing reference point
                                                                                                )
                 )

l_expectedSimulationTime = list()

message( 'lets go' )
for( l_machine in list('supermuc') ) {
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

#message( 'generating weak-scaling plot for cubes' )
#printWeakScaling( list('supermuc'), list( 1024, 2048, 4096, 9216 ), 'statistics_cubes' )

#message( 'generating strong-scaling plot for Landers' )
#printStrongScaling( list('supermuc'), list( 768, 1024, 1536, 2048, 3072, 4096, 6144, 9216, 12288, 24756 ), 'statistics_Landers191M_02_05_02' )

message( 'generating strong-scaling plot for LOH1' )
printStrongScaling( list('supermuc'), list( 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024 ), 'statistics_LOH1_small' )


