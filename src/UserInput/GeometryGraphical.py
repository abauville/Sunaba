import matplotlib.pyplot as plt
import numpy as np
from InputDef import *

## Geometry ##
class Geom_Circle(Geom_Circle):
    def plot(self):
        phi = np.linspace(0,2*pi,30)
        x = self.cx + self.radius*np.cos(phi)
        y = self.cy + self.radius*np.sin(phi)
        plt.plot(x,y)

class Geom_Rect(Geom_Rect):
    def plot(self):
        llx = self.llx
        lly = self.lly
        width = self.width
        height = self.height
        plt.plot([llx,llx+width,llx+width,llx,llx],[lly,lly,lly+height,lly+height,lly])    

class Geom_Line(Geom_Line):
    def plot(self):
        if self.definedFor == "y":
            x = np.array([self.min,self.max])
            y = np.array([self.a*x[0]+self.b, self.a*x[1]+self.b])
        elif self.definedFor == "x":
            y = np.array([self.min,self.max])
            x = np.array([self.a*y[0]+self.b, self.a*y[1]+self.b])

        plt.plot(x,y)

         

class Geom_Sine(Geom_Sine):
    # if definedFor "y":
    # y = base + amplitude*sin(wavelength*x*2*pi+ wavephase)
    def plot(self):
        if self.definedFor == "y":
           x = np.linspace(self.min,self.max,30)
           y = self.base +self.amplitude*np.sin(self.wavelength*x*2*pi + self.wavephase)
        elif self.definedFor == "y":
           y = np.linspace(self.min,self.max,30)
           x = self.base +self.amplitude*np.sin(self.wavelength*y*2*pi + self.wavephase)

        plt.plot(x,y)
        
        
class Geom_Polygon(Geom_Polygon):
    def plot(self):
        plt.plot([x,x[0]],[y,y[0]])
