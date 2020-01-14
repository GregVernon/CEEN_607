import os
import sys
import numpy

def main():
	## Operate on the currently displayed object
	currentViewport = session.viewports[session.currentViewportName]
	activeObject = currentViewport.displayedObject
	
	if activeObject.__class__.__name__ == "Part":
		exportPartMesh(activeObject)
	elif activeObject.__class__.__name__ == "Assembly":
		exportAssemblyMesh(activeObject)
	

def exportPartMesh(activeObject):
	## PART MESH
	NODES = activeObject.nodes
	ELEM = activeObject.elements
	
	## Get number of nodes in each element
	elemNodes = numpy.zeros(len(ELEM),dtype='int')
	for e in range(0,len(ELEM)):
		elemNodes[e] = len(ELEM[e].connectivity)
	
	xyz = numpy.zeros(shape=(len(NODES),3),dtype='double')
	eCONN = numpy.zeros(shape=(len(ELEM),max(elemNodes)),dtype='int')
	
	for n in range(0,len(NODES)):
		nodeLabel = NODES[n].label
		xyz[nodeLabel-1,:] = NODES[n].coordinates
	
	for e in range(0,len(ELEM)):
		eCONN[e,0:elemNodes[e]] = numpy.array(ELEM[e].connectivity) + 1
	
	fName = activeObject.name + ".nodes"
	numpy.savetxt(fName, xyz,delimiter=",")
	
	fName = activeObject.name + ".econn"
	numpy.savetxt(fName, eCONN, fmt="%i", delimiter=",")

def exportAssemblyMesh(activeObject):
	INSTANCE = activeObject.instances
	for i in range(0,len(INSTANCE)):
		exportPartMesh(INSTANCE[i])

if __name__ == "__main__":
	main()