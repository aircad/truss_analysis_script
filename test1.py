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
MAX_HEIGHT = 9.5 # (cm)
MAX_LENGTH = 40

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
 
    print("Forces in Links (N) " ,Q) # units = newtons
    return Q

def linkEvaluation(linkForces = [], nodeData = [], linkData = [], forceData = []):
    linkArea = abs(linkForces)*SAFE_FAC * 10000/ULT_STREN

    print("Link Areas (cm^2) ", linkArea) # units = cm^2

    linkVolume = []
    linkLength = []

    index = 0
    for link in linkData:
        node1 = nodeData[int(link[0])]
        node2 = nodeData[int(link[1])]
        #length in cm
        length = math.dist([float(node1[0]),float(node1[1]),float(node1[2])], [float(node2[0]),float(node2[1]),float(node2[2])])
        #print("Node1:", node1, " Node2:", node2, " length:", length)        
        linkLength.append(round(length,4))
        linkVolume.append(round(linkArea[index]*length,4))

    print("Link Lengths (cm) ", linkLength) # units = cm
    print("Link Volumes (cm^3) ", linkVolume) # units = cm^3

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

    print("Mass (g) ", trussMass)
    print("Strength to Mass", strengthWeightRatio)

    return linkArea, linkLength, trussMass, strengthWeightRatio


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
"""
print(nodeData)
print(linkData)
print(forceData)
print(endData)
"""

#creating truss processing file
filePrep(nodeData, linkData)

#Finding forces in links
linkForces = processTruss(endData, forceData, linkData)

#optimizing truss cross sections and calculating efficiency
linkArea, linkLen, trussMass, str_to_mass = linkEvaluation(linkForces, nodeData, linkData, forceData)

