# -*- coding: utf-8 -*-
"""
Created on Mon Mar 28 13:14:53 2022

@author: Oliver
"""
import Animation as ani
import numpy as np
import matplotlib.pyplot as plt
import itertools

#2000000 3.302

class Simulation(object):
    #Initialising variables
    def __init__(self,t,i):
        self.timestep = t
        self.num_iterations = i
        self.body_list = []
        self.patch_list = []
        #gravitational constant used throughout the program
        self.gravity = 6.673*(10**-11)
        #This class uses a series of arrays to store position information and total energy data as part of an animating over an existing data set
        self.book_of_movement = []
        self.energyGrid_x_Axis = []
        self.energyGrid_y_Axis = []
        #a simplified converted that takes the timestamp a planet needs to orbit the sun then converts said timestamp to earth days
        self.coversionToEarthDaysMultiplier = (365/992.7)
    
    
    #initial Velocity defined by the equation sqrt(G*(mass of sun)/positional vector between the sun and given planet)
    def initialVelocityY(self,m,R_ij):
        return np.sqrt(((self.gravity)*m)/R_ij)
    #reads data from a text file, the process works via a "filtration" system of for loops, first reading all of the raw data and storing into an array, which is then
    #sanitized and organised into respect planets and their information. Information is then assigned data types relating to their function.
    
    #note: all masses are recorded in 10^23 scale
    def read_input_data(self):
        rawInputList = []
        with open("Planets.txt","r") as fi:
            for ln in fi:
                rawInputList.append(ln[3:-1])
        organisedData = []
        temporaryArray = []
        for i in range(len(rawInputList)):
            if (rawInputList[i] == ''):
                organisedData.append(temporaryArray)
                temporaryArray = []
            else:
                temporaryArray.append(rawInputList[i])
        #the sun is filtered first as it is a special case, particularly due to the fact that it's calculation for initial velocity will lead to a division of zero,
        #hence it's initial velocity will be of vector (0 0)
        sun = organisedData[0].copy()
        sunInitialPosition = np.array([float(sun[3]),float(sun[4])])
        sunMass = float(sun[2])
        sunInitialVelocity = np.array([0,0])
        #information regarding calculation is stored in body_list and graphical represent information is stored in patch_list
        self.body_list.append(body(str(sun[0]),str(sun[1]),sunMass,sunInitialPosition,sunInitialVelocity,0,0))
        self.patch_list.append(plt.Circle((float(sun[3]),float(sun[4])),float(sun[5]),color = str(sun[1]), animated = True))  
        organisedData.pop(0)
        #the other planets in the text file are then sanitized and stored into their respective arrays as classes or patches respectively
        for b in range(len(organisedData)):
            currentPlanet = organisedData[b]
            name = str(currentPlanet[0])
            colour = str(currentPlanet[1])
            mass = float(currentPlanet[2])
            initialP = np.array([float(currentPlanet[3]),float(currentPlanet[4])])
            #a previous and current acceleration of 0 is constant for all planets
            previousAcceleration = 0
            currentAcceleration = 0
            vR_ij = initialP - sunInitialPosition
            #x value of the position vector between a given planet and it's sun, used in the calculation for intiial velocity
            xR_ij = vR_ij[0]
            #all planets initial position is on the x-axis of our plot, hence a velocity vector value x of 0 will remain constant
            initialVelocity = np.array([0,(self.initialVelocityY(sunMass,xR_ij))])
            self.body_list.append(body(name,colour,mass,initialP,initialVelocity,currentAcceleration,previousAcceleration))
            #radius is not accurate to real world planet size 
            radius = float(currentPlanet[5])
            self.patch_list.append(plt.Circle((float(currentPlanet[3]),float(currentPlanet[4])),radius,color = str(currentPlanet[1]), animated = True))   
    
    
    #method to output total energy to a text file in the folder, due to it's design before every use please clear content in text file if there is any.
    def OutputtingTotalEnergy(self,t,j):
        f = open("totalEnergyOfSystem.txt", "a")
        output = "time(t): " + str(t) + " | energy (joules): " + str(j) +"\n"
        f.write(output)
        f.close()
    #Potential and kinetic energy are calculated to find total energy of system then appended to a data set for future graphing
    #the for loop ensures all iterations are taken into account with these iterations being scaled by the timestep to graph relevative to the given timestep
    def calc_tot_system_energy(self):
        for i in range(self.num_iterations):
            PE = self.calc_PE()
            KE = self.calc_KE()
            system_tot_energy = self.calc_tot_energy(PE,KE)
            self.energyGrid_y_Axis.append(system_tot_energy)
            xAxis = i * self.timestep
            self.energyGrid_x_Axis.append(xAxis)
            self.OutputtingTotalEnergy(xAxis, system_tot_energy)
            
    #when called will append the satellite to mars to the body and patch lists to be used in launch
    def introduceSatellite(self):
        name = "AFL_Satellite"
        colour = "darkgray"
        #based upon the dry mass of the perservence satellite
        mass = 1.025e-20
        #initial x pos is situated just off the surface of earth
        initialxPos = 1.6
        initialyPos = 0
        initialP = np.array([initialxPos,initialyPos])
        #initial velocity was determined through experimentation
        initialVelocity = np.array([0.0045,0.0076])
        self.body_list.append(body(name,colour,mass,initialP,initialVelocity,0,0))
        self.patch_list.append(plt.Circle((initialxPos,initialyPos),0.01,color = 'darkgray', animated = True))   
        
    #graph of energy findings
    def LineEnergyGraph(self):
        self.calc_tot_system_energy()
        plt.plot(self.energyGrid_x_Axis,self.energyGrid_y_Axis)
        plt.title('Conservation of energy in solar System')
        plt.xlabel('time(t)')
        plt.ylabel('Energy (Joules)')
        plt.show()
    

    
    #run simulation works via a for loop iterating over every planet including the sun and over all iterations to sequencing append position data to the book of motion
    #during this process the read input data method is declared at the top and when experiment 3 is in place, the satellite is introduced.
    #The method ensures new position, new acceleration and velocity are sequentially calculated in that order. 
    #During iterations the method will also check if the planet has made a full orbit via method check_oribitial_period in the body class.
    def run_simulation(self,satelliteLaunch):
        self.read_input_data()
        #a boolean statement as 1 or 0, to determine whether to launch the satellite from earth
        if (satelliteLaunch == 1):
            self.introduceSatellite()
        elif (satelliteLaunch == 0):
            pass
        planets = self.body_list
        planetQuan = len(self.body_list)
        for i in range (planetQuan):
            recOfMotion = []
            currentPlanet = planets[i]
            #produces of copy of other planets to be used in calculating the acceleration of the current planet
            otherPlanets = planets.copy()
            otherPlanets.pop(i)
            initialPosition = getattr(currentPlanet,'position')
            recOfMotion.append(initialPosition)
            #Using the assumption that the previous acceleration of the current acceleration at timestep t will be approximately identical
            currentAcceleration = self.calc_acceleration(planetQuan,currentPlanet,otherPlanets,initialPosition)
            previousAcceleration = currentAcceleration
            #increment counter for total orbits
            totalOrbits = 0
            for x in range (self.num_iterations):
                #Calculations for planet's Velocity and Position
                newPosition = currentPlanet.update_position(self.timestep,currentAcceleration,previousAcceleration)
                recOfMotion.append(newPosition)
                currentPlanetOrbitChecker = currentPlanet.check_orbital_period(recOfMotion[x-1])
                if (currentPlanetOrbitChecker == 1):
                    totalOrbits += currentPlanetOrbitChecker
                else:
                    pass
                newAcceleration = self.calc_acceleration(planetQuan,currentPlanet,otherPlanets,newPosition)
                currentPlanet.update_velocity(newAcceleration,self.timestep)
            print("Body:",getattr(currentPlanet,'name'),"| total solar orbits:", ((1/2) * totalOrbits))
            #positional data of all planets is appended to the array
            self.book_of_movement.append(recOfMotion)

    #Method determines when the satellite was at it's closest to mars returning the timestep it was at this position
    #using a switch system if a new vector magnitude is smaller than the one previous
    def closestToMars(self):
        marsPositions = self.book_of_movement[4]
        satellitePositions = self.book_of_movement[5]
        
        closestMagPosition = (np.linalg.norm(marsPositions[0] - satellitePositions[0]))
        closestTimeStep = 0
        for x in range(self.num_iterations):
            testingMagPosition = (np.linalg.norm(marsPositions[x] - satellitePositions[x]))
            if (testingMagPosition <= closestMagPosition):
                closestMagPosition = testingMagPosition
                closestTimeStep = x
            else:
                pass
        return self.displayTotalEarthDays(self.timestep, closestTimeStep)
    
    #returns the total number of earth days the simulation ran for
    def displayTotalEarthDays(self,t,i):
        return int((t * i) * self.coversionToEarthDaysMultiplier)
    
    #acceleration is calculated via the gravity force law, with parameters of a current planet and all other planets, including the current planet's position. With this information
    #a new acceleration is calculated for the current planet using the masses and positions of all other planets, making a multi body system.
    def calc_acceleration(self,mQuan,currentPlanet,otherplanets,planetPosition):
        massesQuantity = len(otherplanets)
        G = self.gravity
        sigSum = 0
        currentPlanetPosition = planetPosition
        #print(getattr(currentPlanet,'name'))
        for j in range(massesQuantity):
            pairPlanet = otherplanets[j]
            pairPlanetMass = getattr(pairPlanet,'mass') 
            #print(pairPlanetMass)
            pairPlanetPosition = getattr(pairPlanet,'position')                
            diffOfPosition = currentPlanetPosition - pairPlanetPosition
            #Calculation
            magDiffOfPosition = (np.linalg.norm(diffOfPosition))
            sumFrac = (pairPlanetMass/(magDiffOfPosition**3))
            acSum = sumFrac*diffOfPosition
            sigSum += acSum
        return (-G * sigSum)
    
    #cals a class and method from a seperate file to run the animation for organisation reasons
    def animationRunner(self):
        ani.SolarAnimation(self.book_of_movement,self.patch_list, self.num_iterations).run()
    #sum of all kinetic energies of the system, using the planet's mass and current velocity at time t
    def calc_KE(self):
        tKE = 0
        planets = self.body_list
        k = len(planets)
        for i in range(k):            
            Mi = getattr(planets[i],'mass')
            velocityVector = getattr(planets[i],'velocity')
            #Vi calculated by finding the magnitude of the current time velocityVector
            Vi = (np.linalg.norm(velocityVector))
            keSum = (1/2)*Mi*(Vi**2)
            tKE += keSum
        return tKE
    
    
    #used a nested for loop to iterate over all planets of the system to calculate potential energy in pairs between planets
    #itertools.product produces a nested for loop using only a single use, hence the use of it's application here
    def calc_PE(self):
        PE = 0
        k = (len(self.body_list))
        #Mi = getattr(currentPlanet,'mass')
        for i, j in itertools.product(range(k), range (k)):
            #ensures the same planet doesn't find the potential energy of itself and cause an error to the method
            if (i == j):
                pass
            else:   
                iPlanet = self.body_list[i]
                jPlanet = self.body_list[j]
                iMass = getattr(iPlanet,'mass')
                jMass = getattr(jPlanet,'mass')
                iPosition = getattr(iPlanet,'position')
                jPosition = getattr(jPlanet,'position')
                numerator = (self.gravity)*iMass*jMass
                #R_ij is calculated by the difference of the vectors and finding it's magnitude
                denominator = (np.linalg.norm(jPosition - iPosition))
                PE += (numerator/denominator)  
        #the sum is halfed to remove doubling the number of pairs calculated
        return -(1/2)*PE
    #sum of Kinetic and potential energy
    def calc_tot_energy(self,KE,PE):
        return KE + PE   
 
class body():
    def __init__(self,n,c,m,p,v,CA,PA):
        self.name = n
        self.colour = c
        self.mass = m
        self.position = p
        self.velocity = v
        self.current_acceleration = CA
        self.previous_acceleration = PA
    
    #future position is calculated via current position and current velocity, once future position is calculated the current positon will then be updated with this value.
    def update_position(self,Δt,CA,PA):
        Ri = self.position
        Vi = self.velocity
        #equatin to calculate position
        position = Ri + Vi*Δt + (1/6)*((4*CA - PA))*(Δt**2) 
        self.position = position
        return position

    #future velocity takes into account previous,current and future acceleration including the given timestep. Velocity will then be updated with the future velocity
    def update_velocity(self,FA,Δt):
        Vi = self.velocity
        PA = self.previous_acceleration
        CA = self.current_acceleration 
        #two lines of the same equation for debugging purposes
        calPartOne = ((2*FA)+(5*CA)-PA)*(Δt)
        velocity =  Vi + (1/6)*calPartOne
        self.velocity = velocity
        #during the calculation of velocity as it is the final step in the 3 step loop, this step changes current acceleration to previous and future to current acceleration
        self.previous_acceleration = CA
        self.current_acceleration = FA
        
    #check_orbitial_period assesses the signs of the x and y values of posiiton vectors, when the condition is met the future returns a 1 to be incremented as a complete orbit.
    def check_orbital_period(self,previousPosition):
        previousXPos = previousPosition[0]
        previousYPos = previousPosition[1]
        currentPos = self.position
        currentXPos = currentPos[0]
        currentYPos = currentPos[1]
        
        if ((previousXPos > 0) and (previousYPos < 0) and (currentXPos > 0) and (currentYPos > 0)):
            return 1
        else:
            return 0
        
#experimenting class for system
class solarSystem:
    def run(self):
        t = 0.2
        i = (7000)
        planets = Simulation(t,i)
        planets.run_simulation(1)
        print("days to get to mars:",planets.closestToMars())
        print("this simulation ran for approximately", planets.displayTotalEarthDays(t,i),"earth days")
        #animation and line energy graph causes conflict when run together, run either or at a time
        #total energy data will also be written to a file when LineEnergyGraph() is called
        planets.animationRunner()
        #planets.LineEnergyGraph()

def main():
    grid = solarSystem()
    grid.run()
                    
main()
       