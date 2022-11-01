# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 12:20:35 2022

@author: Oliver
"""

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

class energyConservation(object):
    def __init__(self,energyOutputRecords):
        self.graphData = energyOutputRecords
        
        

class SolarAnimation(object):

    def __init__(self,motionRecord,patch_List, niter):
        # set data set of planet's position and their graphical representations
        self.dbOfMotion = motionRecord
        self.planets = patch_List

        # set up simulation parameters
        self.niter = niter
    
    def animate(self, i):
        # update the position of the circles based off the previous calculated data set
        for p in range(len(self.planets)):
            currentPlanet = self.dbOfMotion[p]
            xpos = currentPlanet[i][0]
            ypos = currentPlanet[i][1]
            self.planets[p].center = (xpos,ypos)
        return self.planets


    def run(self):
        # create plot elements
        fig = plt.figure()
        ax = plt.axes()
        
        # add circles to axes
        for i in range(0, len(self.planets)):
            ax.add_patch(self.planets[i])

        # set up the axes
        ax.axis('scaled')
        ax.set_xlim(-2.5, 2.5)
        ax.set_ylim(-2.5, 2.5)

        # create the animator
        self.anim = FuncAnimation(fig, self.animate, frames = self.niter, repeat = False, interval = 0.05, blit = True)

        # show the plot
        plt.show()