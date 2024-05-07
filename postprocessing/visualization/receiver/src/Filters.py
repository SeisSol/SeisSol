##
# @file
# This file is part of SeisSol.
#
# @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
#
# @section LICENSE
# Copyright (c) 2016, SeisSol Group
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

from PyQt5.QtCore import *
from PyQt5.QtWidgets import *
import math
import scipy.signal
import numpy
import abc

class Filter(QGroupBox):
  activeItemChanged = pyqtSignal(name='filterChanged')

  def __init__(self, title, parent = None):
    super(Filter, self).__init__(title, parent)
    
    self.setCheckable(True)
    self.setChecked(False)
    
    self.toggled.connect(self.filterChanged)
    
  @abc.abstractmethod
  def apply(self, wf):
    """ Filter waveforms. """
    return

class Lowpass(Filter):
  def __init__(self, parent = None):
    super(Lowpass, self).__init__('Lowpass', parent)

    lowpassOrderLabel = QLabel('Order', self)
    self.lowpassOrder = QSpinBox(self)
    self.lowpassOrder.setValue(4)
    self.lowpassOrder.valueChanged.connect(self.filterChanged)
    cutoffLabel = QLabel('Cutoff (Hz)', self)
    self.cutoff = QDoubleSpinBox(self)
    self.cutoff.setValue(0.5)
    self.cutoff.setSingleStep(0.1)
    self.cutoff.valueChanged.connect(self.filterChanged)

    filterLayout = QFormLayout(self)
    filterLayout.addRow(lowpassOrderLabel, self.lowpassOrder)
    filterLayout.addRow(cutoffLabel, self.cutoff)
    
  def apply(self, wf):
    Fs = 1.0 / (wf.time[1] - wf.time[0])
    cutoff = min(0.5 * Fs, max(1. / Fs, self.cutoff.value()))
    b, a = scipy.signal.butter(self.lowpassOrder.value(), cutoff, 'low', fs=Fs)
    for name in wf.waveforms.keys():
      wf.waveforms[name] = scipy.signal.filtfilt(b, a, wf.waveforms[name])
      
class Deconvolve(Filter):
  def __init__(self, parent = None):
    super(Deconvolve, self).__init__('Deconvolve', parent)
    
    inLabel = QLabel('Input', self)
    inFun = QLabel('t/T^2 * exp(-t/T)', self)
    TLabel = QLabel('T', self)
    self.T = QDoubleSpinBox(self)
    self.T.setValue(0.1)
    self.T.setMaximum(float('inf'))
    self.T.valueChanged.connect(self.filterChanged)

    outLabel = QLabel('Output', self)
    outFun = QLabel('Gauss', self)
    sigmaLabel = QLabel('Sigma', self)
    self.sigma = QDoubleSpinBox(self)
    self.sigma.setValue(0.05)
    self.sigma.setMaximum(float('inf'))
    self.sigma.valueChanged.connect(self.filterChanged)
    
    filterLayout = QFormLayout(self)
    filterLayout.addRow(inLabel, inFun)
    filterLayout.addRow(TLabel, self.T)
    filterLayout.addRow(outLabel, outFun)
    filterLayout.addRow(sigmaLabel, self.sigma)
    
  def deconv(self, waveform, dt):
    y = numpy.zeros(waveform.size)
    T = self.T.value()
    s = self.sigma.value()
    F = lambda t : numpy.exp(-numpy.square(t-4*s) / (2*s**2)) / (math.sqrt(2*math.pi) * s)
    G = lambda t : 1 - 2*T/s**2 * (t-4*s) - T**2/s**2 * (1 - numpy.square(t-4*s) / s**2)
    eps = 1e-16
    th = s * (4 + math.sqrt(-2 * math.log(eps))) # exp(-(t-4*s)^2 / (2*s**2)) is smaller than eps if t >= th

    # convolves waveform with F*T. The integral int f(t) dt is approximated as simple and stupid as sum_i f(t_i) * dt
    for j in range(y.size):
      i = numpy.arange(max(0, int(j-th/dt)), j)
      t = (j - i) * dt
      f = numpy.multiply(F(t), G(t))
      y[j] = numpy.sum(dt * numpy.multiply(waveform[i], f))

    return y
      
  
  def apply(self, wf):
    dt = wf.time[1] - wf.time[0]
    keys = [f'v{i+1}' for i in range(3)]
    for k in keys:
      if k in wf.waveforms:
        wf.waveforms[k] = self.deconv(wf.waveforms[k], dt)

class Rotate(Filter):
  def __init__(self, parent = None):
    super(Rotate, self).__init__('Rotate', parent)
    
    coordsysLabel = QLabel('Coordinates', self)
    self.coordsys = QComboBox(self)
    self.coordsys.addItem('ned')
    self.coordsys.addItem('seu')
    self.coordsys.currentIndexChanged.connect(self.filterChanged)
        
    epicenterXLabel = QLabel('Epicenter X', self)
    self.epicenterX = QDoubleSpinBox(self)
    self.epicenterX.setValue(0.0)
    self.epicenterX.setMaximum(float('inf'))
    self.epicenterX.valueChanged.connect(self.filterChanged)
    epicenterYLabel = QLabel('Epicenter Y', self)
    self.epicenterY = QDoubleSpinBox(self)
    self.epicenterY.setValue(0.0)
    self.epicenterY.setMaximum(float('inf'))
    self.epicenterY.valueChanged.connect(self.filterChanged)
    
    filterLayout = QFormLayout(self)
    filterLayout.addRow(coordsysLabel, self.coordsys)
    filterLayout.addRow(epicenterXLabel, self.epicenterX)
    filterLayout.addRow(epicenterYLabel, self.epicenterY)
    
  def apply(self, wf):
    if 'v1' in wf.waveforms and 'v2' in wf.waveforms and 'v3' in wf.waveforms:
      epicenter = numpy.array([self.epicenterX.value(), self.epicenterY.value(), 0.0])
      radial = wf.coordinates - epicenter
      phi = math.acos(radial[0] / numpy.linalg.norm(radial))
      
      u = wf.waveforms.pop('v1')
      v = wf.waveforms.pop('v2')
      w = wf.waveforms.pop('v3')
      
      if self.coordsys.currentText() == 'seu':
        u = -u
        w = -w
      
      wf.waveforms['radial'] = math.cos(phi) * u + math.sin(phi) * v
      wf.waveforms['transverse'] = -math.sin(phi) * u + math.cos(phi) * v
      wf.waveforms['vertical'] = w
      for new_comp in ['radial', 'transverse', 'vertical']:
        wf.show[new_comp] = True

class Pick(Filter):
  def __init__(self, parent = None):
    super(Pick, self).__init__('Pick components', parent)
    self.cb_widget_list = []
    self.layout = QGridLayout()
    self.setLayout(self.layout)

  def create_checkbox(self, name):
    widget = QCheckBox(name)
    widget.stateChanged.connect(self.filterChanged)
    return widget

  def add_widgts_to_layout(self):
    cols = 3
    nwidget_over_3 = len(self.cb_widget_list) // cols
    rows = nwidget_over_3 + 1
    for k, widget in enumerate(self.cb_widget_list):
      row = k // cols
      col = k - row * cols
      self.layout.addWidget(widget, row, col)

  def create_checkboxes(self, wf):
    assert self.cb_widget_list == []

    for name in wf.waveforms.keys():
      self.cb_widget_list.append(self.create_checkbox(name))

    self.add_widgts_to_layout()

  def update_checkboxes(self, wf):
    existing_widgets = [cb.text() for cb in self.cb_widget_list]
    required_widgets = list(wf.waveforms.keys())

    # remove not needed checkboxes
    # can't remove from list, while iterating over list
    widgets_to_remove = []
    for w in self.cb_widget_list:
      if not w.text() in required_widgets:
        widgets_to_remove.append(w)
    for w in widgets_to_remove:
      w.deleteLater()
      self.cb_widget_list.remove(w)

    # add new checkboxes
    for name in required_widgets:
      if not name in existing_widgets:
        self.cb_widget_list.append(self.create_checkbox(name))

    self.add_widgts_to_layout()

  def apply(self, wf):
    if not self.cb_widget_list:
      self.create_checkboxes(wf)
    else:
      self.update_checkboxes(wf)

    for widget in self.cb_widget_list:
      var_name = widget.text()
      wf.show[var_name] = widget.isChecked()
