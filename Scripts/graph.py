import numpy as np
from tqdm import tqdm
import re
def levenshtein(seq1, seq2):
    size_x = len(seq1) + 1
    size_y = len(seq2) + 1
    matrix = np.zeros ((size_x, size_y))
    for x in xrange(size_x):
        matrix [x, 0] = x
    for y in xrange(size_y):
        matrix [0, y] = y

    for x in xrange(1, size_x):
        for y in xrange(1, size_y):
            if seq1[x-1] == seq2[y-1]:
                matrix [x,y] = min(
                    matrix[x-1, y] + 1,
                    matrix[x-1, y-1],
                    matrix[x, y-1] + 1
                )
            else:
                matrix [x,y] = min(
                    matrix[x-1,y] + 1,
                    matrix[x-1,y-1] + 1,
                    matrix[x,y-1] + 1
                )
    return (matrix[size_x - 1, size_y - 1])

patterns = []

alphabests= [0,1,2]

def genSeq(prefix, k):
	if (k==0):
		patterns.append(prefix)
		return

	for i in range(0,len(alphabests)):
		newPrefix = prefix + [alphabests[i]]
		genSeq(newPrefix,k-1)
	return




# LEnodes = genfromtxt('LowH_Hyp.csv', delimiter=',')
#         LEnodes = sigs[1:,:-1]



class Node:
  def __init__(self, pval=np.nan, ACE=np.nan, n= np.nan, isLowEnt = False):
    self.pval = pval
    self.ACE = ACE
    self.n=n    
    self.isLowEnt = isLowEnt  

def genrateAllNodes(n=4):
    genSeq([],n)
    dic={}
    for idx,pattern in enumerate(patterns):
        dic[re.sub("[^0-9]", "", str(pattern))]= Node()
    return dic

def readExecutedLowEntNodes(nodes):
    with open("result.txt") as f:
        for line in f:
            parts = line.split(",")
            pattern = re.sub("[^0-9]", "", parts[0].split(":")[1])
            pval =  parts[1].split("=")[1].strip()
            ACE= parts[2].split("=")[1].strip()
            n= parts[3].split("=")[1].strip()
            nodes [ pattern] = Node(pval, ACE, n, True)
    return nodes

def createGraph(nodes):
    graph={}
    for i, pattern_i in nodes.iteritems(): 
        adjList = list()
        for j, pattern_j in nodes.iteritems():
            if (levenshtein(i, j)==1):
                adjList.append(j)
        graph[i] = adjList
    return graph

def writeGraphtoFile(graph, nodes):
    with open("graph.gv", "w") as f:
        f.write("strict graph G {\n")
        for pattern, node in graph.iteritems():
            if not (nodes[pattern].isLowEnt):
                continue
            if  nodes[pattern].isLowEnt:   
                f.write('\t{} [color="red",label="{}, pval:{}"];\n'.format(pattern, pattern, nodes[pattern].pval))
            else:
                f.write('\t{} [label="{}, pval:{}"];\n'.format(pattern, pattern, nodes[pattern].pval))
        for pattern, adjList in graph.iteritems(): 
            for adjPattern in adjList:
                # print pattern, adjPattern
                # print 
                if not (nodes[pattern].isLowEnt) or not  (nodes[adjPattern].isLowEnt) :
                    continue
                f.write("\t{} -- {};\n".format(pattern, adjPattern))
        f.write("}")


def main():
    nodes = genrateAllNodes(3)
    # A  = np.array(patterns)
    # np.savetxt("LowH_Hyp.csv", A.astype(int), delimiter=",",fmt='%s')
    nodes = readExecutedLowEntNodes(nodes)
    graph = createGraph(nodes)
    writeGraphtoFile(graph, nodes)



main()
# with open("graph.gv", "w") as f:
#     f.write("graph G {\n")
#     for idx, pattern in dic.iteritems():
#         f.write('\t{} [label="{}"];\n'.format(idx, pattern))
#     for i in tqdm(range(0,len(nodes))): 
#         for j in range(i,len(nodes)):
#             str1 = re.sub("[^0-9]", "", str(dic[i]))
#             str2 = re.sub("[^0-9]", "", str(dic[j]))
#             if (levenshtein(str1,str2)==1):
#                 f.write("\t{} -- {};\n".format(i,j))
#     f.write("}")




