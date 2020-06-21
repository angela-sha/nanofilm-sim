import math
import random
import sys
sys.path = sys.path + ['/home/cmliepold/Dans_GPU/md_engine-master7/build_tmp/python/build/lib.linux-x86_64-2.7']
from Sim import *
import numpy as np

state = State()
state.deviceManager.setDevice(0)
state.bounds = Bounds(state, lo = Vector(0, 0, 0), hi = Vector(110.26, 173.27, 110.26))
state.rCut = 2.5
state.padding = 0.6
state.periodicInterval = 7
state.shoutEvery = 100

state.atomParams.addSpecies(handle='gold', mass=1, atomicNum=1)
nonbond = FixLJCut(state, 'cut')
nonbond.setParameter('eps', 'gold', 'gold', 1)
nonbond.setParameter('sig', 'gold', 'gold', 1)
state.activateFix(nonbond)

xlist = random.sample(range(math.ceil(state.bounds.lo[0]+1), 
	math.floor(state.bounds.hi[0], 0.01)), 10000)
ylist = random.sample(range(math.ceil(state.bounds.lo[1]+1), 
	math.floor(state.bounds.hi[1], 0.01)), 10000)
for i in range(10000):
	state.addAtom('gold', Vector(xlist[i],ylist[i],110.25/2.), 0.0)
# add 10,000 particles - reduced density of ~0.8

InitializeAtoms.initTemp(state, 'all', temp)
fixNVT = FixNoseHoover(state, 'temp', 'all', temp, timeConstant=1)
state.activateFix(fixNVT)

wallDist = 110.26/2.
topWall = FixWallHarmonic(state, handle="topwall", groupHandle='all', origin=Vector(0, 0, state.bounds.hi[2]), forceDir=Vector(0, 0, -1),dist=wallDist, k=100)
bottomWall = FixWallHarmonic(state, handle="bottomwall", groupHandle='all',origin=Vector(0, 0,state.bounds.lo[2]),forceDir=Vector(0, 0, 1), dist=wallDist, k=1000)
state.activateFix(topWall)
state.activateFix(bottomWall)
deform = FixDeform(state, handle='def', groupHandle='all', deformRate=1, multiplier=Vector(-0.166, 0, 0), applyEvery=1)

pressure = FixPressureBerendsen(state, 'constP', .2, 10, 1)
integVerlet = IntegratorVerlet(state)

writeconfig = WriteConfig(state, fn="Sim1", writeEvery=100, format='xyz', handle='writer')
state.activateWriteConfig(writeconfig)

writeconfig.writeEvery = 100
integVerlet.run(4000)

state.activateFix(deform)
writeconfig.writeEvery = 100
integVerlet.run(53125)

