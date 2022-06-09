#GRAVITATIONAL CONSTANT
#G = 6.6743 * 10**(-15) #mass in kg, dist in meters, time in seconds
G = 1.0 #???

#LAW OF MOTION
def LOM_Force(m1, m2, r, G=G):
    return G*m1*m2/(r**2)