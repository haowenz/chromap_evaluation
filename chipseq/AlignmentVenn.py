#!/usr/bin/env python3

import sys

if (len(sys.argv) == 1):
	print("usage: a.py bwa.sam bowtie.sam chromap.paf > intersection_summary.out")
	sys.exit(1)

def IsMatch(a, b):
	if (a[1] == -1 or b[1] == -1):
		return False
	if (a[0] != b[0]):
		return False
	if (a[1] + 10 < b[1] or a[1] - 10 > b[1]):
		return False
	return True

#fp = open(sys.argv[1])
#chrId = {}
#for line in fp:
#	if (line[1:3] != "SQ"):
#		continue
#	chrName = line.split()[1].split(":")[1]
#	chrId[chrName] = len(chrId)
#fp.close()

result = [0, 0, 0, 0, 0, 0, 0, 0] # bit encoding. 0-bwa, 1-bowtie, 2-chromap
fps = []
for tool in [1, 2, 3]:
	fps.append(open(sys.argv[tool]))
alignments = {}
for tool in [0, 1]:
	for line in fps[tool]:
		if (line[0] == "@"):
			continue

		cols = line.rstrip().split()
		flag = int(cols[1])
		if (flag & 0x900 != 0):
			continue
		if (flag & 0x4 != 0):
			continue
		readid = cols[0]
		if (flag & 0x40 != 0):
			readid += "/1"
		else:
			readid += "/2"
		chrid = cols[2]
		pos = int(cols[3]) - 1
		if (readid not in alignments):
			alignments[readid] = [["-1", -1], ["-1", -1], ["-1", -1]]
		alignments[readid][tool] = [chrid, pos]
	print("#Done with aligner " + str(tool))
hasSuffix = True
lineCnt = 0
for line in fps[2]:
	cols = line.rstrip().split()
	readid = cols[0]
	chrid = cols[5]
	pos = int(cols[7])
	if (hasSuffix == False):
		if (lineCnt%2 == 0):
			readid += "/2"
		else:
			readid += "/1"
	if (readid not in alignments):
		alignments[readid] = [["-1", -1], ["-1", -1], ["-1", -1]]
	alignments[readid][2] = [chrid, pos]
	lineCnt += 1
print("#Done with aligners")

for readid in alignments:
	a = alignments[readid]
	hasAlign = []
	for i in range(3):
		if (a[i][1] != -1):
			hasAlign.append(i)

	flag = 0
	if (len(hasAlign) == 0):
		flag = 0 # ignore the case that none of the aligner mapped the read. Shouldn't happen though
	elif (len(hasAlign) == 1):
		flag = (1<<hasAlign[0])
		result[flag] += 1
	elif (len(hasAlign) == 2):
		if (IsMatch(a[hasAlign[0]], a[hasAlign[1]])):
			flag = (1<<hasAlign[0]) | (1<<hasAlign[1])
			result[flag] += 1
		else:
			result[1<<hasAlign[0]] += 1
			result[1<<hasAlign[1]] += 1
	else:
		for i, j in [[0, 1], [1, 2], [0, 2]]:
			if (IsMatch(a[i], a[j])):
				flag |= (1<<hasAlign[i]) | (1<<hasAlign[j])
		if (flag == 0):
			result[1] += 1
			result[2] += 1
			result[4] += 1
		elif (flag == 7):
			result[flag] += 1
		else: # two of the aligners agree with each other
			result[flag] += 1
			result[7^flag] += 1
	#print(readid, flag)	

print("\t".join([str(x) for x in result]))

for fp in fps:
	fp.close()
