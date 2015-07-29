##
# @file
# This file is part of SeisSol.
#
# @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
#
# @section LICENSE
# Copyright (c) 2015, SeisSol Group
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
#

from PyQt4.QtGui import *
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
import matplotlib.pyplot as plt

import math
import Navigation


class View(QWidget):

  def __init__(self, parent = None):
    super(View, self).__init__(parent)

    self.figure = plt.figure()
    self.canvas = FigureCanvas(self.figure)
    self.canvas.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
    
    layout = QHBoxLayout(self)
    self.navigations = []
    for i in range(0, 2):
      navigation = Navigation.Navigation(self)
      navigation.activeItemChanged.connect(self.plot)
      navigation.setSizePolicy(QSizePolicy.Maximum, QSizePolicy.Minimum)
      layout.addWidget(navigation)
      self.navigations.append(navigation)

    layout.addWidget(self.canvas)

  def plot(self):
    waveforms = []
    numPlots = 0
    names = set()
    for nav in self.navigations:
      for wf in nav.getActiveWaveforms():
        numPlots = max(numPlots, len(wf.names))
        names.update(set(wf.names))
        waveforms.append(wf)
    names = list(names)
    names.sort()

    wf = waveforms[0]

    self.figure.clear()
    numRows = math.ceil(math.sqrt(numPlots));
    numCols = math.ceil(numPlots / numRows)
    subplots = dict()
    for i in range(0, len(names)):
      subplots[ names[i] ] = self.figure.add_subplot(numRows, numCols, i+1)

    for wf in waveforms:
      for name in wf.names:
        p = subplots[name]
        p.plot(wf.time, wf.waveforms[name])
        p.set_xlabel('t')
        p.set_ylabel(name)

    self.figure.tight_layout()
    self.canvas.draw()
