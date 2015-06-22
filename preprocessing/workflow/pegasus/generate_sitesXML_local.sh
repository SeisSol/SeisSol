#!/bin/bash
##
# @file
# This file is part of SeisSol.
#
# @author Fabio Gratl (f.gratl AT in.tum.de)
# @author Alex Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
#
# @section LICENSE
# Copyright (c) 2013-2014, SeisSol Group
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

if [[ $# == 0 ]]
then 
  echo "using default name 'sites.xml'"
  outputFile="sites.xml"
  echo "using current working directory as path to run folder"
  runFolderDir=$(pwd)/run
  echo "using current working directory as path to output folder"
  outputFolderDir=$(pwd)/output
  echo "using \$HOME/.ssh/id_rsa as path to private ssh key"
  keyPath=$HOME/.ssh/id_rsa
  
elif [[ $# == 4 ]]
then
  runFolderDir=$1
  outputFolderDir=$2
  keyPath=$3
  outputFile=$4
elif [[ $# != 4 ]]
then
  echo "usage: ./generate_sitesXML_local.sh [/absolute/Path/To/Run/Folder] [/absolute/Path/To/Output/Folder] [/absolute/Path/To/Private/Key] [NameForOutputFile]"
  exit
fi

echo "current run directory is ${runFolderDir}"
echo "current output directory is ${outputFolderDir}"
echo "current key path is ${keyPath}"
echo "writing" $outputFile

echo -e "<?xml version=\"1.0\" encoding=\"UTF-8\"?>

<!-- This stuff at the beginning is just boilerplate -->
<sitecatalog xmlns=\"http://pegasus.isi.edu/schema/sitecatalog\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://pegasus.isi.edu/schema/sitecatalog http://pegasus.isi.edu/schema/sc-4.0.xsd\" version=\"4.0\">

    <!-- The local site contains information about the submit host -->
    <!-- The arch and os keywords are used to match binaries in the transformation catalog -->
    <site handle=\"local\" arch=\"x86_64\" os=\"LINUX\">

        <!-- These are the paths on the submit host were Pegasus stores data -->
        <!-- Scratch is where temporary files go -->
		<directory type=\"shared-scratch\" path=\"${runFolderDir}\">
			<file-server operation=\"all\" url=\"file://${runFolderDir}\"/>
        </directory>
        <!-- Storage is where pegasus stores output files -->
		<directory type=\"local-storage\" path=\"${outputFolderDir}\">
			<file-server operation=\"all\" url=\"file://${outputFolderDir}\"/>
        </directory>

        <!-- This profile tells Pegasus where to find the user's private key for SCP transfers -->
        <profile namespace=\"env\" key=\"SSH_PRIVATE_KEY\">${keyPath}</profile>
        <profile namespace=\"env\" key=\"LD_LIBRARY_PATH\" >${LD_LIBRARY_PATH}</profile>
    </site>
</sitecatalog>" > $outputFile
