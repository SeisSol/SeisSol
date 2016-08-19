##
# @file
# This file is part of SeisSol.
#
# @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
# @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
#
# @section LICENSE
# Copyright (c) 2015-2016, SeisSol Group
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
try:
	from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
except ImportError:
	from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
import matplotlib.pyplot as plt

import Navigation
import Filters
import Watchdog
import re
import math
import numpy
import os.path
import scipy.fftpack

class View(QWidget):

  def __init__(self, parent = None):
    super(View, self).__init__(parent)
    
    self.__watchdog = Watchdog.Watchdog()
    self.__watchdog.fileChanged.connect(self.refreshAll)

    self.figure = plt.figure()
    self.canvas = FigureCanvas(self.figure)
    self.canvas.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)    
    toolbar = NavigationToolbar(self.canvas, self)
    
    self.navigationLayout = QHBoxLayout()
    
    layout = QHBoxLayout(self)
    self.navigations = []
    self.addNavigation(True)
    
    self.filters = [ Filters.Lowpass(), Filters.Deconvolve(), Filters.Rotate() ]
    filterLayout = QVBoxLayout()    
    for f in self.filters:
      filterLayout.addWidget(f)
      f.filterChanged.connect(self.plot)
    filterLayout.addStretch()
    
    addIcon = QIcon.fromTheme('list-add')
    addNaviButton = QPushButton(addIcon, 'Add navigation', self)
    addNaviButton.clicked.connect(self.addNavigation)
    
    self.maxFreq = QDoubleSpinBox(self)
    self.maxFreq.setValue(10.0)
    self.maxFreq.setVisible(False)
    self.maxFreq.valueChanged.connect(self.plot)
    spectrumIcon = QIcon.fromTheme('network-wireless')
    self.spectrum = QPushButton(spectrumIcon, 'Spectrum', self)
    self.spectrum.setCheckable(True)
    self.spectrum.clicked.connect(self.plot)
    self.spectrum.toggled.connect(self.maxFreq.setVisible)
    self.diff = QPushButton('Diff', self)
    self.diff.setCheckable(True)
    self.diff.clicked.connect(self.plot)
    self.diff.clicked.connect(self.spectrum.setHidden)
    self.spectrum.toggled.connect(self.diff.setHidden)
    
    autoRefresh = QPushButton(QIcon.fromTheme('view-refresh'), 'Auto', self)
    autoRefresh.setCheckable(True)
    autoRefresh.clicked.connect(self.__watchdog.toggle)
    
    saveAll = QPushButton(QIcon.fromTheme('document-save'), '', self)
    saveAll.clicked.connect(self.savePlots)

    toolLayout = QHBoxLayout()
    toolLayout.addWidget(addNaviButton)
    toolLayout.addWidget(self.diff)
    toolLayout.addWidget(self.spectrum)
    toolLayout.addWidget(self.maxFreq)
    toolLayout.addWidget(autoRefresh)
    toolLayout.addWidget(saveAll)
    toolLayout.addWidget(toolbar)
    plotLayout = QVBoxLayout()
    plotLayout.addLayout(toolLayout)
    plotLayout.addWidget(self.canvas)
    layout.addLayout(self.navigationLayout)
    layout.addLayout(plotLayout)
    layout.addLayout(filterLayout)
    
  def addNavigation(self, noclose = False):
    navigation = Navigation.Navigation(noclose)
    navigation.activeItemChanged.connect(self.plot)
    navigation.folderChanged.connect(self.navigationFolderChanged)
    navigation.close.connect(self.closeNavigation)
    navigation.setSizePolicy(QSizePolicy.Maximum, QSizePolicy.Minimum)
    self.navigationLayout.addWidget(navigation)
    self.navigations.append(navigation)

  def navigationFolderChanged(self, oldFolder, newFolder):
    self.__watchdog.removeFolder(oldFolder)
    self.__watchdog.addFolder(newFolder)
    
  def closeNavigation(self, widget):
    self.navigations.remove(widget)
    self.navigationLayout.removeWidget(widget)
    widget.deleteLater()
    self.plot()
    
  def refreshAll(self):
    for navigation in self.navigations:
        navigation.refreshFolder()

  def plot(self):
    wfc = [wf for nav in self.navigations for wf in nav.getActiveWaveforms()]
    for filt in self.filters:
      if filt.isChecked():
        for wf in wfc:
          filt.apply(wf)
          
    if self.diff.isChecked() and len(wfc) > 0:
      wf0 = wfc.pop()
      for nWf, wf in enumerate(wfc):
        wfc[nWf].subtract(wf0)

    names = set([name for wf in wfc for name in wf.waveforms.iterkeys()])
    numPlots = len(names)

    self.figure.clear()
    if numPlots > 0:
      names = list(names)
      names.sort()

      numRows = math.ceil(math.sqrt(numPlots));
      numCols = math.ceil(numPlots / numRows)
      subplots = dict()
      for i in range(len(names)):
        subplots[ names[i] ] = self.figure.add_subplot(numRows, numCols, i+1)

      for nWf, wf in enumerate(wfc):
        for name, waveform in wf.waveforms.iteritems():
          p = subplots[name]
          if self.spectrum.isChecked():
            n = len(waveform)
            dt = wf.time[1]-wf.time[0] # assume equally spaced samples
            f = scipy.fftpack.fftfreq(n, dt)
            W = dt * scipy.fftpack.fft(waveform)
            maxFreqIndices = numpy.argwhere(f > self.maxFreq.value())
            L = maxFreqIndices[0] if len(maxFreqIndices) > 0 else n/2
            p.loglog(f[1:L], numpy.absolute(W[1:L]), label=str(nWf))
            p.set_xlabel('f [Hz]')
          elif self.diff.isChecked():
            p.plot(wf.time, waveform, label='{}-0'.format(nWf+1))
            p.set_xlabel('t (s)')
          else:
            p.plot(wf.time, waveform, label=str(nWf))
            p.set_xlabel('t (s)')
          p.set_ylabel(name)

      self.figure.tight_layout()

    for i in range(len(names)):
      subplots[ names[i] ].legend(prop={'size':8}, frameon=False)
    self.canvas.draw()
  
  def savePlots(self):
    filetypes = self.canvas.get_supported_filetypes_grouped()
    defaultFiletype = self.canvas.get_default_filetype()
    filters = []
    selectedFilter = ''
    for name, extensions in sorted(filetypes.iteritems()):
      filtr = '{0} ({1})'.format(name, ' '.join(['*.{0}'.format(ext) for ext in extensions]))
      if defaultFiletype in extensions:
        selectedFilter = filtr
      filters.append(filtr)
    fileName, filtr = QFileDialog.getSaveFileNameAndFilter(self, 'Choose a save location.', '', ';;'.join(filters), selectedFilter)
    fileName = os.path.splitext(str(fileName))[0]
    extension = re.search(r'\*(\.[a-zA-Z]+)', str(filtr)).group(1)
    
    maxRow = min([nav.numberOfRows() for nav in self.navigations])
    for row in range(maxRow):
      for nav in self.navigations:
        nav.selectWaveformAt(row)
      self.plot()
      self.canvas.print_figure('{0}{1:03}{2}'.format(fileName, row+1, extension))
