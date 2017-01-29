# -*- coding: utf-8 -*-
import sys
import argparse
import math
# http://braeunig.us/space/index.htm
# http://braeunig.us/space/comb.htm

# Gas constant (J/(K*kmol))
R = 8314.4598

def exhaustVelocity(chamberP, chamberT, k, molarM, exitP):
	v = math.sqrt((2*k/(k-1))*(R*chamberT/molarM)*(1.0-(exitP/chamberP)**((k-1.0)/k)))
	return v

def throatP(chamberP, k):
	p = chamberP*((1.0 + (k-1)/2.0)**(-k/(k-1.0)))
	return p

def throatT(chamberT, k):
	t = chamberT/(1.0 + (k-1.0)/2.0)
	return t

def MachNum(chamberP, k, exitP):
	Mn = math.sqrt( (2.0/(k-1.0))*((chamberP/exitP)**((k-1.0)/k) - 1.0) )
	return Mn

def exitArea(chamberP, chamberT, k, molarM, exitP, A_t):
	Mn = MachNum(chamberP, k, exitP)
	#print ("MachN: " + str(Mn))
	A_e = (A_t/Mn)*( (1.0 + (k-1.0)/2.0 * Mn*Mn)/((k+1.0)/2.0) )**((k+1.0)/(2.0*(k-1.0)))
	return A_e

def main():

	# Parse commandline
	cmdLine = argparse.ArgumentParser()

	cmdLine.add_argument("massFlow", type=float, help="Thrust chamber mass flow")
	cmdLine.add_argument("chamberPress", type=float, help="Thrust chamber pressure")
	cmdLine.add_argument("temp", type=float, help="Thrust chamber temperature")
	cmdLine.add_argument("k", type=float, help="Spesific heat ratio")
	cmdLine.add_argument("areaRatio", type=float, help="Nozzle area ratio")
	cmdLine.add_argument("molarMass", type=float, help="Exhaust molar mass")

	args = cmdLine.parse_args()

	lowerP = 0.0
	# exit pressure can't be higher than throat pressure
	upperP = throatP(args.chamberPress, args.k)
	testP = (lowerP+upperP)/2.0
	#print ("starting p: " + str(testP))

	# Throat parameters
	P_t = throatP(args.chamberPress, args.k)
	T_t = throatT(args.temp, args.k)
	A_t = args.massFlow/P_t * math.sqrt(R*T_t/(args.molarMass*args.k))

	# Loop to find correct pressure on nozzle-exit
	iters = 0;
	while (1):
		iters = iters + 1
		v_e = exhaustVelocity(args.chamberPress, args.temp, args.k, args.molarMass, testP)
		A_e = exitArea(args.chamberPress, args.temp, args.k, args.molarMass, testP, A_t)

		# Do a binary search to find correct exit-pressure for the given area ratio
		if (A_e > A_t*args.areaRatio):
			lowerP = testP
			testP = (upperP+testP)/2.0
		else:
			upperP = testP
			testP = (lowerP+testP)/2.0

		#print ("New p: " + str(testP) + " area ratio with prev iter: " + str(A_e/A_t))

		# Exit loop when acceptable precision is reached
		if (math.fabs(args.areaRatio-(A_e/A_t)) < 0.01):
			break

	#print ("Iterations: " + str(iters))

	thrust_m = v_e * args.massFlow
	thrust_s = thrust_m + (testP - 101325)*A_e
	thrust_v = thrust_m + testP*A_e

	# 80% bell-nozzle length
	# Length percentage is relative to 15 deg half-angle conical nozzle
	aRatio = A_e/A_t
	n_radius = math.sqrt(A_e/math.pi)
	n_length = n_radius/math.tan(math.radians(15))*(math.sqrt(aRatio)-1)/math.sqrt(aRatio)
	n_length = n_length*0.8


	print ("Exhaust velocity: {0:.1f}".format(v_e) + " m/s")
	print ("Throat: {0:.6f}".format(A_t) + " m² {0:.0f}".format(T_t) + " K")
	print ("Exit: {0:.6f}".format(A_e) + " m² Area ratio: {0:.3f}".format(aRatio))
	print ("Pressure at nozzle exit: {0:.0f}".format(testP) + " Pa")
	print ("\nThrust")
	print (" Sea level: {0:.0f}".format(thrust_s) + " N")
	print (" Vacuum:    {0:.0f}".format(thrust_v) + " N")

	print ("\nBell-nozzle 80%:\n Length:   {0:.3f}".format(n_length) + " m")
	print (" Diameter: {0:.3f}".format(n_radius*2) + " m")


if __name__ == '__main__':
        try:
                main()
        except KeyboardInterrupt:
                sys.exit(1)
