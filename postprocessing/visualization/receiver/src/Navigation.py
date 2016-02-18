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

from PyQt4.QtCore import *
from PyQt4.QtGui import *
import os
import copy

import Tecplot
import Waveform

class Navigation(QWidget):
  activeItemChanged = pyqtSignal(name='activeItemChanged')
  close = pyqtSignal(QWidget, name='close')

  def __init__(self, noclose = False, parent = None):
    super(Navigation, self).__init__(parent)
    
    self.currentFolder = ''

    openIcon = QIcon.fromTheme('folder-open')
    openButton = QPushButton(openIcon, '', self)
    openButton.clicked.connect(self.selectFolder)
    refreshIcon = QIcon.fromTheme('view-refresh')
    refreshButton = QPushButton(refreshIcon, '', self)
    refreshButton.clicked.connect(self.refreshFolder)
    if not noclose:
      closeIcon = QIcon.fromTheme('window-close')
      closeButton = QPushButton(closeIcon, '', self)
      closeButton.clicked.connect(self.emitClose)
    
    self.receiverList = QListView(self)
    self.model = QStandardItemModel()
    self.receiverList.setModel(self.model)
    self.receiverList.clicked.connect(self.activeItemChanged)
    
    buttonLayout = QHBoxLayout()
    buttonLayout.addWidget(openButton)
    buttonLayout.addWidget(refreshButton)
    if not noclose:
      buttonLayout.addWidget(closeButton)
    buttonLayout.addStretch()

    layout = QVBoxLayout(self)
    layout.addLayout(buttonLayout)
    layout.addWidget(self.receiverList)

  def selectFolder(self):
    folder = QFileDialog.getExistingDirectory(self, 'Open directory', self.currentFolder, QFileDialog.ShowDirsOnly)
    if not folder.isEmpty():
      self.currentFolder = str(folder)
      self.readFolder(self.currentFolder)
      
  def readFolder(self, folder):
    if len(folder) != 0:
      self.model.clear()
      files = [f for f in os.listdir(folder) if os.path.isfile(os.path.join(folder,f))]
      files.sort()
      for f in files:
        wf = Tecplot.read(os.path.join(folder,f))
        item = QStandardItem(f)
        item.setData(wf)
        self.model.appendRow(item)
    
  def getActiveWaveforms(self):      
    waveforms = []
    for index in self.receiverList.selectedIndexes():
      wf = self.model.itemFromIndex(index).data().toPyObject()
      waveforms.append( copy.deepcopy(wf) )
    return waveforms
    
  def refreshFolder(self):
    self.readFolder(self.currentFolder)
    self.activeItemChanged.emit()
    
  def numberOfRows(self):
    return self.model.rowCount() 
    
  def selectWaveformAt(self, row):
    self.receiverList.selectionModel().select(self.model.index(row, 0), QItemSelectionModel.ClearAndSelect)
  
  def emitClose(self):
    self.close.emit(self)
    


