#t_gcc_procyte: 713.48 s

import time
sim = theSimulator.createStepper('SpatiocyteStepper', 'SS')
sim.VoxelRadius = 2.5e-9 
sim.SearchVacant = 0

theSimulator.rootSystem.StepperID = 'SS'
theSimulator.createEntity('Variable', 'Variable:/:GEOMETRY').Value = 0
theSimulator.createEntity('Variable', 'Variable:/:LENGTHX').Value = 1.15e-6
theSimulator.createEntity('Variable', 'Variable:/:LENGTHY').Value = 1.15e-6
theSimulator.createEntity('Variable', 'Variable:/:LENGTHZ').Value = 1.15e-6
theSimulator.createEntity('Variable', 'Variable:/:VACANT')
theSimulator.createEntity('Variable', 'Variable:/:XYPLANE').Value = 0
theSimulator.createEntity('Variable', 'Variable:/:XZPLANE').Value = 0
theSimulator.createEntity('Variable', 'Variable:/:YZPLANE').Value = 0

theSimulator.createEntity('System', 'System:/:Surface').StepperID = 'SS'
theSimulator.createEntity('Variable', 'Variable:/Surface:DIMENSION').Value = 2
theSimulator.createEntity('Variable', 'Variable:/Surface:VACANT')

theSimulator.createEntity('Variable', 'Variable:/:A').Value = 160000

#logger = theSimulator.createEntity('VisualizationLogProcess', 'Process:/:logger')
#logger.VariableReferenceList = [['_', 'Variable:/Surface:VACANT'], ['_', 'Variable:/:A']]

populator = theSimulator.createEntity('MoleculePopulateProcess', 'Process:/:pop')
populator.VariableReferenceList = [['_', 'Variable:/:A']]

diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffuseA')
diffuser.VariableReferenceList = [['_', 'Variable:/:A']]
diffuser.D = 1e-12

run(0.00005)
start = time.time()
run(0.1)
end = time.time()
duration = end-start
print duration
