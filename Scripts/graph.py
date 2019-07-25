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

nodes = []

alphabests= [0,1,2]

def genSeq(prefix, k):
	if (k==0):
		nodes.append(prefix)
		return

	for i in range(0,len(alphabests)):
		newPrefix = prefix + [alphabests[i]]
		genSeq(newPrefix,k-1)
	return


# genSeq([],2)
# print nodes
# print len(nodes)


# print levenshtein("1222", "1221")



LEnodes = genfromtxt('LowH_Hyp.csv', delimiter=',')
        LEnodes = sigs[1:,:-1]

genSeq([],4)
dic={}
for idx,node in enumerate(nodes):
    dic[idx]=node

with open("graph.gv", "w") as f:
    f.write("graph G {\n")
    for idx, pattern in dic.iteritems():
        f.write('\t{} [label="{}"];\n'.format(idx, pattern))
    for i in tqdm(range(0,len(nodes))): 
        for j in range(i,len(nodes)):
            str1 = re.sub("[^0-9]", "", str(dic[i]))
            str2 = re.sub("[^0-9]", "", str(dic[j]))
            if (levenshtein(str1,str2)==1):
                f.write("\t{} -- {};\n".format(i,j))
    f.write("}")

