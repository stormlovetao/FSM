
fp = open("/n/data2/bwh/nrlgy/scherzer/twang/FSM/RolladVidal")

id_dic = {}
ID = 1
for line in fp:
	linesplit = line.strip().split("\t")
	protein1 = linesplit[0]
	protein2 = linesplit[1]
	protein1_split = protein1.split("|")
	protein2_split = protein2.split("|")
	switch = 0
	for item in protein1_split:
		if item not in id_dic:
		 	id_dic[item] = ID
		 	switch = 1
	if switch ==1:
		ID += 1
		switch = 0
	for item in protein2_split:
		if item not in id_dic:
		 	id_dic[item] = ID
		 	switch = 1
	if switch ==1:
		ID += 1
		switch = 0

fp.close()
fp = open("/n/data2/bwh/nrlgy/scherzer/twang/FSM/RolladVidal")
outfile = open("/n/data2/bwh/nrlgy/scherzer/twang/FSM/RolladVidal.reformat.txt", 'w')
outfile.write(str(ID-1)+"\n") 
for line in fp:
	linesplit = line.strip().split("\t")
	protein1 = linesplit[0]
	protein2 = linesplit[1]
	protein1_split = protein1.split("|")
	protein2_split = protein2.split("|")
	outline = str(id_dic[protein1_split[0]]) + " " + str(id_dic[protein2_split[0]]) + "\n"
	outfile.write(outline)

outfile.close()


