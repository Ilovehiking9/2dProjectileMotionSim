import matplotlib.pyplot as plt
import math
import time
import csv
#density of steel
#1.5 diameter
#

def createUnitVector(angleVector):
    #angle vector should be in (magnitude, angle), or polar format
    
    unitVector = []

    x = math.cos(math.radians(angleVector[1]))*angleVector[0]
    y = math.sin(math.radians(angleVector[1]))*angleVector[0]

    unitVector.append(x)
    unitVector.append(y)

    return unitVector


#for 2d world spaces ONLY (if i see any 3d shit im gonna shoot someone)
class object():
    

    def __init__(self, s0 = [0, 0], v0 = [0, 0], a0 = [0, 0], m = 0, dragCoefficient = 0, crossArea = 0):
        
        self.displacement = s0
        self.velocity = v0
        self.acceleration = a0  
        self.mass = m 
        self.dragCoefficient = dragCoefficient
        self.crossArea = crossArea

    def setDisplacement(self, angleVector):
        self.displacement = angleVector

    def setVelocity(self, angleVector):
        self.velocity = angleVector

    def setAcceleration(self, angleVector):
        self.acceleration = angleVector


    def getAirResistance(self, velocityVector):
        airResistance = []
        for dimensionalVelocity in velocityVector:
            directionalAirResistance = -0.5*1.293*self.dragCoefficient*self.crossArea*(dimensionalVelocity**2)
            if dimensionalVelocity < 0:
                directionalAirResistance *= -1
            airResistance.append(directionalAirResistance)

        return airResistance

    def getDisplacement(self, timeInterval, x0, v0, a0):
        
        displacementX = 0.5 * a0[0] * timeInterval ** 2 + v0[0] * timeInterval + x0[0]
        displacementY = 0.5 * a0[1] * timeInterval ** 2 + v0[1] * timeInterval + x0[1]

        return [displacementX, displacementY]

    def getVelocity(self, timeInterval, v0, a0):

        velocityX = a0[0] * timeInterval + v0[0] 
        velocityY = a0[1] * timeInterval + v0[1]

        return [velocityX, velocityY]

    def getAcceleration(self, netForce):
        acceleration = [netForce[0]/self.mass, netForce[1]/self.mass]
        return acceleration

    
    def getForceOfGravity(self):
        return [0, -9.81 * self.mass]

    def runSimulation(self, timeInterval = 0.001, duration = 1):

        simulationData = []

        #initialization
        displacement = self.displacement
        velocity = self.velocity
        acceleration = self.acceleration
        t = 0
        timeStep = 0

        simulationData.append({"time": t, "displacement": displacement, "velocity": velocity, "acceleration": acceleration})
        print(simulationData)

        while t <= duration:
            timeStep += 1
            t += timeInterval

            #calculate forces
            #gravity + air resistance + ##user forces later##
            forces = []

            #forces.append(self.getAirResistance(velocity))


            forces.append(self.getForceOfGravity())
            
            netForce = [0, 0]
            for force in forces:
                netForce[0] += force[0]
                netForce[1] += force[1]

            displacement = self.getDisplacement(timeInterval, simulationData[-1]["displacement"], simulationData[-1]["velocity"], simulationData[-1]["acceleration"])
            velocity = self.getVelocity(timeInterval, simulationData[-1]["velocity"], simulationData[-1]["acceleration"])
            acceleration = self.getAcceleration(netForce)

            simulationData.append({"time": t, "displacement": displacement, "velocity": velocity, "acceleration": acceleration})

            

        time = []
        displacementx = []
        displacementy = []




        with open("data.csv", "w") as data:
            data.writelines(f"Time, X Displacement, Y Displacement, X Velocity, Y Velocity, X Acceleration, Y Acceleration\n")

            for i in simulationData:
                
                data.writelines(f"{i["time"]}, {i["displacement"][0]}, {i["displacement"][1]}, {i["velocity"][0]}, {i["velocity"][1]}, {i["acceleration"][0]}, {i["acceleration"][1]}\n")

                displacementx.append(i["displacement"][0])
                displacementy.append(i["displacement"][1])
        
        plt.plot(displacementx, displacementy)
        plt.xlabel('Horizontal Displacement (m)')
        plt.ylabel("Vertical Displacement (m)")

        plt.show()

            



if __name__ == "__main__":

    steelBall = object(a0 = [0, -0.134397], 
                       v0=createUnitVector([4.878, 45]), 
                       s0 = [0, 0.146 + 0.908], m=13.7/1000, 
                       crossArea=(math.pi*((0.75/100)**2)), 
                       dragCoefficient=0.47)
    

    steelBall.runSimulation()

