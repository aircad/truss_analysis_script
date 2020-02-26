from sapy import displmethod
from sapy import element
from sapy import gmsh
from sapy import structure
from sapy import plotter
import math

"""
-----------
DESCRIPTION
-----------
"""

"""
Will calculate forces, required link cross sections and lengths, and total truss mass and efficiency

Positions in csv are in cm
mass returned in kg
forces returned in N
dimensions returned in m units

"""

#link material attributes
HOOKES = 3390000000  #(Pa)
ULT_STREN = 24860000  #(Pa)
SAFE_FAC = 1
DENSITY = 0.141 #(g/cm^3)
GEO_FILE_NAME = "test"
MAX_HEIGHT = 9 # (cm)
MAX_LENGTH = 39.5
COMPRESSION_COMPENSATION = 4
ITERATION_SIZE = 1 # (cm)
NUM_ITERATION = 1

"""
---------
FUNCTIONS
---------
"""

# Creates .geo file for truss calc
def filePrep(nodeData = [], linkData = []):
    #converting into .geo file
    test_file = open(GEO_FILE_NAME + ".geo","w+")

    #points:
    index = 1
    for point in nodeData:    
        test_file.write("Point(" + str(index) + ") = {" + point[0] + "," + point[1] + "," + point[2] + "};\n")
        index += 1
        
    index = 1
    for line in linkData:
        test_file.write("Line(" + str(index) + ") = {" + str(int(line[0]) + 1) + "," + str(int(line[1]) + 1) + "};\n")
        index += 1
   
    test_file.close()

# Processes Truss and outputs forces
def processTruss(endData = [], forceData = {}, linkData = []):
    #processing the truss
    file = GEO_FILE_NAME

    #setting supports
    bound = {}
    pinnedSupport = False
    for support in endData:
        if not pinnedSupport:
            bound[int(support)] = [1,1] # 1 = not free to move, 0 = free to move along axis
            pinnedSupport = True
        else:
            bound[int(support)] = [0,1]

    #bound = {0:[1,1], 1:[0,1]}
            
    #applying load
    load = {}
    forceDictKeys = forceData.keys()

    for key in forceDictKeys:
        load[int(key)] = [int(forceData[key][0]), int(forceData[key][1])] # applies load on node n-1 from .geo file aka node n from .csv
    
    #load = {2: [1,0]}

    #print(bound)
    #print(load)
                     
    ele = element.Data()

    for i in range(len(linkData)):
        ele.E[i] = 1.
        ele.A[i] = 1.
        ele.TYPE[i] = 'Truss'

    mesh = gmsh.Parse(file)
    model = structure.Builder(mesh, ele, bound) 

    U, Q = displmethod.solver(mesh,model,ele,load)
 
    #print("Forces in Links (N) " ,Q) # units = newtons
    return Q

#distance finding function
def findDist(p,q):
    return math.sqrt(sum((px-qx)**2 for px,qx in zip(p,q)))

def linkEvaluation(linkForces = [], nodeData = [], linkData = [], forceData = []):
    
    linkArea = []
    for link in linkForces:
        if link > 0:
            linkArea.append(link*SAFE_FAC*10000/ULT_STREN)
        else:
            linkArea.append(-link*SAFE_FAC*COMPRESSION_COMPENSATION*10000/ULT_STREN)
    
    #linkArea = abs(linkForces)*SAFE_FAC * 10000/ULT_STREN

    #print("Link Areas (cm^2) ", linkArea) # units = cm^2

    linkVolume = []
    linkLength = []

    index = 0
    for link in linkData:
        node1 = nodeData[int(link[0])]
        node2 = nodeData[int(link[1])]
        #length in cm
        
        #length = math.dist([float(node1[0]),float(node1[1]),float(node1[2])], [float(node2[0]),float(node2[1]),float(node2[2])])
        length = findDist([float(node1[0]),float(node1[1]),float(node1[2])], [float(node2[0]),float(node2[1]),float(node2[2])])
        
        #print("Node1:", node1, " Node2:", node2, " length:", length)        
        linkLength.append(round(length,4))
        linkVolume.append(round(linkArea[index]*length,4))

    #print("Link Lengths (cm) ", linkLength) # units = cm
    #print("Link Volumes (cm^3) ", linkVolume) # units = cm^3

    # finding mass of links and truss
    
    linkMass = []
    for link in linkVolume:
        linkMass.append(link*DENSITY) # units = g

    trussMass = 0 # units = g

    for m in linkMass:
        trussMass += m

    # finding load
    temp = forceData.items()
    appLoad = 0
    for i in temp:
        appLoad += abs(int(i[1][1])) # only interested in Y forces

    strengthWeightRatio = round(appLoad*1000/9.81 / trussMass,4) # (g) = 1N*1000/(9.81m/s^2)

    #print("Mass (g) ", trussMass)
    #print("Strength to Mass", strengthWeightRatio)

    return linkArea, linkLength, trussMass, strengthWeightRatio

# outputs a list of possible node configurations
def NodeWiggler(nodes = [], ends = [], forceNode = 0):
    nodeID = []
    index = 0
    # does not allow end wiggling
    for node in nodes:
        if str(index) not in ends:
            nodeID.append(index)
        index+=1
    nodesArray = [[]]
    #accomodating for the possibilities
    for ID in nodeID:
        newNodesArray = []
        # checking for forceNode
        if ID != forceNode:
            for data in nodesArray:
                for x in range((2*NUM_ITERATION+1)**2):
                    newNodesArray.append(data.copy())
            nodesArray = newNodesArray
            #print(nodesArray)
            xIter = -NUM_ITERATION
            yIter = -NUM_ITERATION
            
            for nodeSeries in nodesArray:
                nodeSeries.append((str(float(nodes[ID][0])+(xIter*ITERATION_SIZE)), str(float(nodes[ID][1])+(yIter*ITERATION_SIZE)), nodes[ID][2]))
                xIter += 1
                if xIter == NUM_ITERATION+1:
                    xIter = -NUM_ITERATION
                    yIter += 1
                if yIter == NUM_ITERATION+1:
                    yIter = -NUM_ITERATION
        else:
            for data in nodesArray:
                for x in range(2*NUM_ITERATION+1):
                    newNodesArray.append(data.copy())        
            nodesArray = newNodesArray
            #print(nodesArray)
            yIter = -NUM_ITERATION
            
            for nodeSeries in nodesArray:
                nodeSeries.append((nodes[ID][0], str(float(nodes[ID][1])+(yIter*ITERATION_SIZE)), nodes[ID][2]))
                yIter += 1
                if yIter == NUM_ITERATION+1:
                    yIter = -NUM_ITERATION
                    
    for each in nodesArray:
        assert(nodesArray.count(each) == 1), "oof not unique generated nodes"

    removalList = []
    for nodeSeries in nodesArray:
        # Assumes ends are in ascending order
        for end in ends:
            nodeSeries.insert(int(end),nodes[int(end)])
        
        #removing possibilities that violate constraints
        xMax = -10000
        xMin = 10000
        yMax = -10000
        yMin = 10000
        for node in nodeSeries:
            if float(node[0]) >= xMax:
                xMax = float(node[0])
            if float(node[0]) <= xMin:
                xMin = float(node[0])
            if float(node[1]) >= yMax:
                yMax = float(node[1])
            if float(node[1]) <= yMin:
                yMin = float(node[1])
        if xMax - xMin > MAX_LENGTH or yMax - yMin > MAX_HEIGHT:
            removalList.append(nodeSeries)
    #removing unwanted:
    for each in removalList:
        nodesArray.remove(each)
    #print(nodesArray)
    #for nodeSeries in nodesArray:
        
    return nodesArray
    
"""
----
Main
----
"""

#reading from csv
try: node_Data = open("Nodes.csv", "r")
except: print("File Not Found")

#splitting along lines
data = node_Data.read().strip("ï»¿").splitlines()

#prep for parsing
forceDataStarted = False
linkDataStarted = False
endsDataStarted = False

nodeData = []
linkData = []
forceData = {}
endData = []

for d in data:
    #adding nodes
    #print(d)
    if not forceDataStarted and not linkDataStarted and not endsDataStarted:
        if d.find("Links") == -1:
            temp = d.split(",")
            if temp[0] != "":
                #putting nodes in proper format for processing
                nodeData.append((temp[1],temp[2],temp[3]))
        else:
            linkDataStarted = True
    #adding links
    elif not forceDataStarted and not endsDataStarted:
        if d.find("Forces") == -1:
            temp = d.split(",")
            if temp[0] != "":
                #putting links in proper format for processing
                linkData.append((temp[0],temp[1]))
        else:
            forceDataStarted = True
    #adding forces
    elif not endsDataStarted:
        if d.find("Ends") == -1:
            temp = d.split(",")
            if temp[0] != "":
                #putting forces in proper format for processing
                forceData[temp[0]] = (temp[1],temp[2],temp[3])
        else:
            endsDataStarted = True
    #adding end connections
    else:
        temp = d.split(",")
        endData.append(temp[0])

node_Data.close()

#print(nodeData)
#print(linkData)
#print(forceData)
#print(endData)

#getting force node id
temp = forceData.keys()
forceNode = 0
for first in temp:
    forceNode = int(first)
nodesArray = NodeWiggler(nodeData, endData, forceNode)

print("""
        ---------
        ORIGINAL:
        ---------
        """)
#creating truss processing file
filePrep(nodeData, linkData)

#Finding forces in links
linkForces = processTruss(endData, forceData, linkData)

#optimizing truss cross sections and calculating efficiency
linkArea, linkLen, trussMass, str_to_mass = linkEvaluation(linkForces, nodeData, linkData, forceData)

#outputting results
print("Forces in Links (N) " ,linkForces) # units = newtons
print("Link Areas (cm^2) ", linkArea) # units = cm^2
print("Link Lengths (cm) ", linkLen) # units = cm
print("Mass (g) ", trussMass) # units = g
print("Strength to Mass", str_to_mass) # g/g

#checking iterations:
bestSeries = []
bestForces = []
bestArea = []
bestLen = []
bestMass = 1000000
bestRatio = 0
for nodeSeries in nodesArray:
    #creating truss processing file
    filePrep(nodeSeries, linkData)
    #Finding forces in links
    linkForces = processTruss(endData, forceData, linkData)
    #optimizing truss cross sections and calculating efficiency
    linkArea, linkLen, trussMass, str_to_mass = linkEvaluation(linkForces, nodeSeries, linkData, forceData)
    if trussMass < bestMass:
        bestSeries = nodeSeries
        bestForces = linkForces
        bestArea = linkArea
        bestLen = linkLen
        bestMass = trussMass
        bestRatio = str_to_mass

print("""
        -----
        BEST:
        -----
        """)
index = 0
for node in bestSeries:
    print("node " + str(index) + ": (" + node[0] + "," + node[1] + "," + node[2] + ")") #printing out nodes of best result
    index += 1
    
print("Forces in Links (N) " ,bestForces) # units = newtons
print("Link Areas (cm^2) ", bestArea) # units = cm^2
print("Link Lengths (cm) ", bestLen) # units = cm
print("Mass (g) ", bestMass) # units = g
print("Strength to Mass", bestRatio) # g/g
