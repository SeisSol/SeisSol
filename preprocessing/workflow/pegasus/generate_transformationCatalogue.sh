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
  echo "using default name 'tc.dat'"
  outputFile="tc.dat"
  echo "using current working directory as path to input folder"
  inputFolderDir=$(pwd)/
elif [[ $# == 2 ]]
then
  inputFolderDir=$1
  outputFile=$2
elif [[ $# != 2 ]]
then
  echo "usage: ./generate_transformationCatalogue.sh [/absolute/Path/To/Input/Folder/] [NameForOutputFile]"
  exit
fi

echo "current input directory is: ${inputFolderDir}input/"

#trying to find cube3_stat
path_cube3_stat=$(which cube3_stat)
echo "cube3_stat directory is   : $path_cube3_stat"

echo "writing" $outputFile

echo "If you intend to use git-svn in the workflow set the git path in the output file"

echo -e "# This is the transformation catalog. It lists information about each of the
# executables that are used by the workflow.

tr git {
  site local {
    pfn \"<INSERT/GIT/PATH/HERE>\"
    arch \"x86_64\"
    os \"linux\"
    type \"INSTALLED\"
  }
}

tr svn {
  site local {
    pfn \"/lrz/sys/tools/subversion/1.7.9/bin/svn\"
    arch \"x86_64\"
    os \"linux\"
    type \"INSTALLED\"
  }
}

tr buildscript {
  site local {
    pfn \"${inputFolderDir}input/buildscript.sh\"
    arch \"x86_64\"
    os \"linux\"
    type \"INSTALLED\"
  }
}

tr phi_communication {
  site local {
    pfn \"${inputFolderDir}input/phi_communication.sh\"
    arch \"x86_64\"
    os \"linux\"
    type \"INSTALLED\"
  }
}

tr copyMeshToBench {
  site local {
    pfn \"${inputFolderDir}input/copyMeshToBench.sh\"
    arch \"x86_64\"
    os \"linux\"
    type \"INSTALLED\"
  }
}

tr cube3_stat {
  site local {
    pfn \"${path_cube3_stat}\"
    arch \"x86_64\"
    os \"linux\"
    type \"INSTALLED\"
  }
}

tr compare_cube3_stat {
  site local {
    pfn \"${inputFolderDir}input/compare_cube3_stat.sh\"
    arch \"x86_64\"
    os \"linux\"
    type \"INSTALLED\"
  }
}

tr remove_ranks {
  site local {
    pfn \"${inputFolderDir}input/remove_ranks.sh\"
    arch \"x86_64\"
    os \"linux\"
    type \"INSTALLED\"
  }
}

tr gnuplot{
  site local {
    pfn \"/lrz/sys/graphics/gnuplot/4.6.0/bin/gnuplot\"
    arch \"x86_64\"
    os \"linux\"
    type \"INSTALLED\"
  }
}

tr compare_receivers {
  site local {
    pfn \"/tmp/python/Python-2.7.6/build/bin/python ${inputFolderDir}input/compare_receivers.py\"
    arch \"x86_64\"
    os \"linux\"
    type \"INSTALLED\"
  }
}

tr monitorSlurmState {
  site local {
    pfn \"${inputFolderDir}input/monitorSlurmState.sh\"
    arch \"x86_64\"
    os \"linux\"
    type \"INSTALLED\"
  }
}

tr compare_TetraElastic {
  site local {
    pfn \"${inputFolderDir}input/compare_TetraElastic.py\"
    arch \"x86_64\"
    os \"linux\"
    type \"INSTALLED\"
  }
}

tr mv { 
  site local {
    pfn \"/bin/mv\"
    arch \"x86_64\"
    os \"linux\"
    type \"INSTALLED\"
  }
}

tr ls { 
  site local {
    pfn \"/bin/ls\"
    arch \"x86_64\"
    os \"linux\"
    type \"INSTALLED\"
  }
}

tr tar {
  site local {
    pfn \"/bin/tar\"
    arch \"x86_64\"
    os \"linux\"
    type \"INSTALLED\"
  }
}

tr sbatch {
  site local {
    pfn \"/usr/bin/sbatch\"
    arch \"x86_64\"
    os \"linux\"
    type \"INSTALLED\"
  }
}

tr mkdir {
  site local {
    pfn \"/bin/mkdir\"
    arch \"x86_64\"
    os \"linux\"
    type \"INSTALLED\"
  }
}
tr cp {
  site local {
    pfn \"/bin/cp\"
    arch \"x86_64\"
    os \"linux\"
    type \"INSTALLED\"
  }
}

tr rm {
  site local {
    pfn \"/bin/rm\"
    arch \"x86_64\"
    os \"linux\"
    type \"INSTALLED\"
  }
}

tr cpp {
    site local {
      pfn \"/bin/cpp\"
      arch \"x86_64\"
      os \"linux\"
      type \"INSTALLED\"
    }
}

tr generate_parameter_file {
  site local {
    pfn \"${inputFolderDir}input/generate_parameter_file.sh\"
    arch \"x86_64\"
    os \"linux\"
    type \"INSTALLED\"
  }
}
" > $outputFile
