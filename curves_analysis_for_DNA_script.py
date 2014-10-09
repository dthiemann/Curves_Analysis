#!/usr/bin/python
import os, math, sys, random, subprocess
import commands
import time

#Author: Dylan Thiemann
#Email: dylanthiemann@gmail.com

#Get command line arguments
if len(sys.argv) == 3:
	analysisStartTime = sys.argv[1]
	analysisEndTime = sys.argv[2]
elif len(sys.argv) == 2:
	analysisStartTime = sys.argv[1]
	analysisEndTime = None
elif len(sys.argv) == 1:
	analysisStartTime = 5350
	analysisEndTime = sys.maxint
else:
	print("You have too many arguements, program terminated")
	quit() 

#Pre-defined functions for later use
def meanstdv(x):
    from math import sqrt
    n, mean, std = len(x), 0, 0
    for a in x:
        mean = mean + a
    mean = mean / float(n)
    for a in x:
        std = std + (a - mean)**2
    std = sqrt(std / float(n-1))
    return std

analyze_BP = False
analyze_specifics = raw_input("Would you like to analyze specific nucleotide pairs? (y/n): ")
if analyze_specifics == "y":
	analyze_BP = True

startTime = []
lis_file = open("total_DNA_properties_for_every_snapshot.lis","r")
sequence = []
sequence2 = []

BP_Axis = []
Intra_BP_parameters = []
Inter_BP = []

BP_Axis_std = []
Intra_BP_std = []
Inter_BP_std = []

Backbone_params_strand1 = []
Backbone_params_strand2 = []
numFiles = 0

print("Obtaining data from total_DNA_properties_for_every_snapshot.lis")
line = lis_file.readline()
line_list = line.split()
line_list.append(False)
while (True):
	if line_list[0] == "TITLE":
		if (int(float(line_list[5])) == int(analysisStartTime)):
			break;
		else:
			line = lis_file.readline()
			line_list = line.split()
			line_list.append(False)
	else:
                line = lis_file.readline()
                line_list = line.split()
                line_list.append(False)


exitLoop = False
endReached = False
while (line and not exitLoop):						#Runs untill we reach the end of the file or we hit the end marker
	line_list = line.split()
	line_list.append(False)
	if ((line_list[0] == "TITLE") and line_list[5] == (str(int(analysisEndTime))+".00000")):
		endReached = True
	if (len(startTime) < 2) and (line_list[0] == "TITLE"):
		startTime.append(line_list[5])
	#Handles the BP-Axis data
	if line_list[0] == "(A)":
		numFiles += 1
		line = lis_file.readline() 		#Skips the blank line
		line = lis_file.readline()		 #First line of data
		if not BP_Axis: 			#If the BP-Axis data structure is currently empty
			tempList = line.split()
			tempList.append(False)
			counter0 = 0
			while tempList[0]:	#Continue only untill we hit a blank line
				BP_Axis.append([])
				BP_Axis_std.append([])
				for element in tempList:
					if element == "---":
						BP_Axis_std[counter0].append([0.0])
						BP_Axis[counter0].append(0.0)
					elif element == False:
						continue
					else:
						try:
							BP_Axis_std[counter0].append([float(element)])
							BP_Axis[counter0].append(float(element))
							
						except ValueError:
							BP_Axis_std[counter0].append([element])
							BP_Axis[counter0].append(element)
				tempList = lis_file.readline().split()
				tempList.append(False)
				counter0 += 1
		else: 					#indices of interest = 4, 5 ,6 ,7 ,8
			counter1 = 0
			tempList = line.split()
			tempList.append(False)
			if tempList[4] == "NaN":
				break
			while counter1 < len(BP_Axis): #tempList[0]
				BP_Axis[counter1][4] += float(tempList[4])
				BP_Axis_std[counter1][4].append(float(tempList[4]))
				BP_Axis[counter1][5] += float(tempList[5])
				BP_Axis_std[counter1][5].append(float(tempList[5]))
				BP_Axis[counter1][6] += float(tempList[6])
				BP_Axis_std[counter1][6].append(float(tempList[6]))
				BP_Axis[counter1][7] += float(tempList[7])
				BP_Axis_std[counter1][7].append(float(tempList[7]))
				if tempList[8] == "---":
					BP_Axis[counter1][8] += 0.0
					BP_Axis_std[counter1][8].append(0.0)
				else: 
					BP_Axis[counter1][8] += float(tempList[8])
					BP_Axis_std[counter1][8].append(tempList[8])
				tempList = lis_file.readline().split()
				tempList.append(False)
				counter1 = counter1 + 1

	#Handles the Intra-BP parameters
	if line_list[0] == "(B)":
		line = lis_file.readline() #Read the blank lines inbetween the title and the actual data
		line = lis_file.readline()
		line = lis_file.readline()
		line = lis_file.readline()
		if not Intra_BP_parameters:                         #If the BP-Axis data structure is currently empty
                        tempList = line.split()
                        tempList.append(False)
                        counter0 = 0
                        while tempList[0]:      #Continue only untill we hit a blank line
                                Intra_BP_parameters.append([])
				Intra_BP_std.append([])
                                for element in tempList:
                                        if element == False:
                                                continue
                                        else:
                                                try:
                                                        Intra_BP_parameters[counter0].append(float(element))
							Intra_BP_std[counter0].append([float(element)])
                                                except ValueError:
                                                        Intra_BP_parameters[counter0].append(element)
							Intra_BP_std[counter0].append([element])
                                tempList = lis_file.readline().split()
                                tempList.append(False)
                                counter0 += 1
                else:                                   #indices of interest = 4, 5 ,6 ,7 ,8
                        counter1 = 0
                        tempList = line.split()
                        tempList.append(False)
                        if tempList[4] == "NaN":
                                break
                        while counter1 < len(Intra_BP_parameters): #tempList[0]
                                Intra_BP_parameters[counter1][4] += float(tempList[4])
                                Intra_BP_parameters[counter1][5] += float(tempList[5])
                                Intra_BP_parameters[counter1][6] += float(tempList[6])
                                Intra_BP_parameters[counter1][7] += float(tempList[7])
                                Intra_BP_parameters[counter1][8] += float(tempList[8])
				Intra_BP_parameters[counter1][9] += float(tempList[9])
				Intra_BP_std[counter1][4].append(tempList[4])
				Intra_BP_std[counter1][5].append(tempList[5])
				Intra_BP_std[counter1][6].append(tempList[6])
				Intra_BP_std[counter1][7].append(tempList[7])
				Intra_BP_std[counter1][8].append(tempList[8])
				Intra_BP_std[counter1][9].append(tempList[9])
                                tempList = lis_file.readline().split()
                                tempList.append(False)
                                counter1 = counter1 + 1

	#Handles Inter-BP data
	if line_list[0] == "(C)":
		line = lis_file.readline()
		line = lis_file.readline()
		if not Inter_BP:
			tempList = line.split()
			tempList.append(False) # Fail safe incase the list is empty
			counter0 = 0
			while tempList[0]:
				Inter_BP.append([])
				Inter_BP_std.append([])
                                for element in tempList:
                                        if element == False:
                                                continue
                                        else:
                                                try:
                                                        Inter_BP[counter0].append(float(element))
							Inter_BP_std[counter0].append([float(element)])
                                                except ValueError:
                                                        Inter_BP[counter0].append(element)
							Inter_BP_std[counter0].append([element])
                                tempList = lis_file.readline().split()
                                tempList.append(False)
                                counter0 += 1
		else:                                   #indices of interest = 4, 5 ,6 ,7 ,8
                        counter1 = 0
                        tempList = line.split()
                        tempList.append(False)
                        if tempList[4] == "NaN":
                                break
                        while counter1 < len(Inter_BP): #tempList[0]
                                Inter_BP[counter1][4] += float(tempList[4])
                                Inter_BP[counter1][5] += float(tempList[5])
                                Inter_BP[counter1][6] += float(tempList[6])
                                Inter_BP[counter1][7] += float(tempList[7])
                                Inter_BP[counter1][8] += float(tempList[8])
                                Inter_BP[counter1][9] += float(tempList[9])
				Inter_BP[counter1][10] += float(tempList[10])
				Inter_BP[counter1][11] += float(tempList[11])
				Inter_BP_std[counter1][4].append(float(tempList[4]))
				Inter_BP_std[counter1][5].append(float(tempList[5]))
				Inter_BP_std[counter1][6].append(float(tempList[6]))
				Inter_BP_std[counter1][7].append(float(tempList[7]))
				Inter_BP_std[counter1][8].append(float(tempList[8]))
				Inter_BP_std[counter1][9].append(float(tempList[9]))
				Inter_BP_std[counter1][10].append(float(tempList[10]))
				Inter_BP_std[counter1][11].append(float(tempList[11]))
                                tempList = lis_file.readline().split()
                                tempList.append(False)
                                counter1 = counter1 + 1
	
	#Handles the Backbone Parameters
	if line_list[0] == "(D)":
		if endReached == True:
			exitLoop = True
		line = lis_file.readline()
		line = lis_file.readline() 		#Strand 1 line
		line = lis_file.readline()
		line = lis_file.readline()
		
		if not Backbone_params_strand1:
			tempList = line.split()
			tempList.append(False)
			counter0 = 0
			while tempList[0]:
				Backbone_params_strand1.append([])
				for element in tempList:
					if element != False:
						Backbone_params_strand1[counter0].append(element)
				tempList = lis_file.readline().split()
				tempList.append(False)
				counter0 += 1
			line = lis_file.readline()
			line = lis_file.readline()
			line = lis_file.readline()
			tempList2 = line.split() #Start of strand 2 analysis
			tempList2.append(False)
			counter1 = 0
			while tempList2[0]:
				Backbone_params_strand2.append([])
				for element in tempList2:
					if element != False:
						Backbone_params_strand2[counter1].append(element)	#Adding the values as strings, not ints (IMPORTANT)!
				tempList2 = lis_file.readline().split()
				tempList2.append(False)
				counter1 += 1
		else:
			tempList = line.split()
			tempList.append(False)
			counter0 = 0
			while tempList[0]:
				index = 3
				while index < (len(Backbone_params_strand1[counter0])):
					if tempList[index] == "----":
						newString = " " + tempList[index]
                                                Backbone_params_strand1[counter0][index] += newString
						index += 1
					else:
						newString = " " + tempList[index]
						Backbone_params_strand1[counter0][index] += newString
						index += 1
				tempList = lis_file.readline().split()
				tempList.append(False)
				counter0 += 1
			line = lis_file.readline()
			line = lis_file.readline()
			line = lis_file.readline()
			
			#Continue strand2 analysis
			tempList2 = line.split()
			tempList2.append(False)
			counter1 = 0
			while tempList2[0]:
				index2 = 3
				while index2 < (len(Backbone_params_strand2[counter1])):
					if tempList2[index2] == "----":
						newString = " " + tempList2[index2]
						Backbone_params_strand2[counter1][index2] += newString
						index2 += 1
					else:
						newString = " " + tempList2[index2]
						Backbone_params_strand2[counter1][index2] += newString
						index2 += 1
				tempList2 = lis_file.readline().split()
				tempList2.append(False)
				counter1 += 1
	line = lis_file.readline()
		
numFiles = numFiles - 1 		#There is one blank file at the end which we don't need to account for
lis_file.close()

######################################################
	#Calculating Data for Output File #1
######################################################
BP_Axis_std_final = []
Intra_BP_std_final = []
Inter_BP_std_final = []

counter_BP = 0
for element in BP_Axis_std:
	BP_Axis_std_final.append([])
	BP_Axis_std_final[counter_BP].append(element[0][0])
	BP_Axis_std_final[counter_BP].append(element[1][0])
	BP_Axis_std_final[counter_BP].append(element[2][0])
	BP_Axis_std_final[counter_BP].append(element[3][0])
	index = 4
	while index < len(element):
		try:
			std = meanstdv(map(float,element[index]))
			BP_Axis_std_final[counter_BP].append(std)
			index += 1
		except TypeError:
			BP_Axis_std_final[counter_BP].append(0.0)
			index += 1
	counter_BP += 1

counter_Intra_BP = 0
for element in Intra_BP_std:
        Intra_BP_std_final.append([])
        Intra_BP_std_final[counter_Intra_BP].append(element[0][0])
        Intra_BP_std_final[counter_Intra_BP].append(element[1][0])
        Intra_BP_std_final[counter_Intra_BP].append(element[2][0])
        Intra_BP_std_final[counter_Intra_BP].append(element[3][0])
        index = 4
        while index < len(element):
                try:
                        std = meanstdv(map(float,element[index]))
                        Intra_BP_std_final[counter_Intra_BP].append(std)
                        index += 1
                except TypeError:
                        Intra_BP_std_final[counter_Intra_BP].append(0.0)
                        index += 1
        counter_Intra_BP += 1

counter_Inter_BP = 0
for element in Inter_BP_std:
        Inter_BP_std_final.append([])
        Inter_BP_std_final[counter_Inter_BP].append(element[0][0])
        Inter_BP_std_final[counter_Inter_BP].append(element[1][0])
        Inter_BP_std_final[counter_Inter_BP].append(element[2][0])
        Inter_BP_std_final[counter_Inter_BP].append(element[3][0])
        index = 4
        while index < len(element):
                try:
                        std = meanstdv(element[index])
                        Inter_BP_std_final[counter_Inter_BP].append(std)
                        index += 1
                except TypeError:
                        Inter_BP_std_final[counter_Inter_BP].append(0.0)
                        index += 1
        counter_Inter_BP += 1

for element in BP_Axis:
	sequence.append(element[1])
	index = 4
	while index < len(element):
		element[index] /= float(numFiles)
		index = index + 1
for element in Intra_BP_parameters:
	index = 4
	while index < len(element):
		element[index] /= numFiles
		index = index + 1
for element in Inter_BP:
	index = 4
	while index < len(element):
		element[index] /= numFiles
		index = index + 1

#Gathering backbone parameter data for each strand
#Pucker located at index 12 
print("Calculating pucker frequencies...")
puckerStrand1 = []
puckerStrand2 = []
puckCounter1 = 0
for nucleotide in Backbone_params_strand1:
	puckerStrand1.append([])
	puckerStrand1[puckCounter1].append(nucleotide[1])
	puckerStrand1[puckCounter1].append([])
	line = nucleotide[12].split()
	uniqueElems = set(line)
	for puckers in uniqueElems:
		numPuckers = line.count(puckers)
		newString = puckers + " " + str(numPuckers)
		puckerStrand1[puckCounter1][1].append(newString)
	puckCounter1 += 1

puckCounter2 = 0
for nucleotide in Backbone_params_strand2:
        puckerStrand2.append([])
        sequence2.append(nucleotide[1]) 		#find the sequence of the complimentary strand
        puckerStrand2[puckCounter2].append(nucleotide[1])
	puckerStrand2[puckCounter2].append([])
        line = nucleotide[12].split()
        uniqueElems = set(line)
        for puckers in uniqueElems:
                numPuckers = line.count(puckers)
                newString = puckers + " " + str(numPuckers)
                puckerStrand2[puckCounter2][1].append(newString)
        puckCounter2 += 1
print("Done with the pucker data!")

################################################################################
		#Writing the data to the output file#
################################################################################

print("Writing output file!")
sequenceName = ""
outputFileName = ""
for nuc in sequence:
	sequenceName += nuc
if analysisStartTime and (analysisEndTime and analysisEndTime != sys.maxint):
	outputFile = open("curves+_average_DNA_analysis_for_" + sequenceName + "_from_" + str(analysisStartTime) + "_to_" + str(analysisEndTime) + ".txt", "w")
	outputFileName = "curves+_average_DNA_analysis_for_" + sequenceName + "_from_" + str(analysisStartTime) + "_to_" + str(analysisEndTime) + ".txt"
	print("Output file is called: curves+_average_DNA_analysis_for_" + sequenceName + "_from_" + str(analysisStartTime) + "_to_" + str(analysisEndTime) + ".txt")
elif analysisStartTime and analysisEndTime == sys.maxint:
	outputFile = open("curves+_average_DNA_analysis_for_" + sequenceName + "_from_" + str(analysisStartTime) + "_to_end" + ".txt", "w")
	outputFileName = "curves+_average_DNA_analysis_for_" + sequenceName + "_from_" + str(analysisStartTime) + "_to_end" + ".txt"
	print("Output file is called: curves+_average_DNA_analysis_for_" + sequenceName + "_from_" + str(analysisStartTime) + "_to_end" + ".txt")
elif not analysisStartTime and analysisEndTime:
	outputFile = open("curves+_average_DNA_analysis_for_" + sequenceName + "_from_5350_to_" + str(analysisEndTime) + ".txt", "w")
	outputFileName = "curves+_average_DNA_analysis_for_" + sequenceName + "_from_5350_to_" + str(analysisEndTime) + ".txt"
	print("Output file is called: curves+_average_DNA_analysis_for_" + sequenceName + "_from_5350_to_" + str(analysisEndTime) + ".txt")
else:
	outputFile = open("curves+_average_DNA_analysis_for_" + sequenceName + "_from_5350_to_end.txt", "w")
	outputFileName = "curves+_average_DNA_analysis_for_" + sequenceName + "_from_5350_to_end.txt"
	print("Output file is called: curves+_average_DNA_analysis_for_" + sequenceName + "_from_5350_to_end.txt")

currentDate = time.strftime("%m/%d/%Y")
outputFile.write("\t****Curves Average Analysis****\t\n")
outputFile.write("\t\t" + currentDate + "\t\t\n")
outputFile.write("\n")
outputFile.write("\n")

outputFile.write("The primary sequence is: \t\t" + sequenceName + "\n")
outputFile.write("The number of files being analyzed: \t" + str(numFiles) + "\n")
outputFile.write("\n")

#Write for BP-Axis data
outputFile.write("(A) BP-Axis\n")
outputFile.write("\t\t\t\t" + "Xdisp\t" + "Ydisp\t" + "Incline\t" + "  Tip\t" + "Ax-bend\n")
#outputFile.write("\n")
for BP_data in BP_Axis:
	for data in BP_data:
		if (type(data) == int) or (type(data) == float):
			data = round(data,3)
			if str(data)[0] == "-":
				outputFile.write(str(data) + "\t")
			else:
				outputFile.write(" " + str(data) + "\t") 
		else:
			outputFile.write(str(data) + "\t")
	outputFile.write("\n")
outputFile.write("\n")

outputFile.write("BP_Axis Standard Deviations\n")
for std in BP_Axis_std_final:
	for data in std:
		if (type(data) == int) or (type(data) == float):
			data = round(data,3)
			outputFile.write(" " + str(data) + "\t")
		else:
			outputFile.write(" " + str(data) + "\t")
	outputFile.write("\n")
outputFile.write("\n")

#Write for Intra-BP Parameters
outputFile.write("(B) Intra-BP\n")
outputFile.write("\t\t\t\t" + "Shear\t" + "Stretch\t" + "Stagger\t" + "Buckle\t" + "Propel\t" + "Opening\n")
#outputFile.write("\n")
for BP_data in Intra_BP_parameters:
	for data in BP_data:
		if (type(data) == int) or (type(data) == float):
			data = round(data,3)
                        if str(data)[0] == "-":
                                outputFile.write(str(data) + "\t")
                        else:
                                outputFile.write(" " + str(data) + "\t")
		else:
			outputFile.write(data + "\t")
	outputFile.write("\n")
outputFile.write("\n")

outputFile.write("Intra_BP Standard Deviations\n")
for std in Intra_BP_std_final:
        for data in std:
                if (type(data) == int) or (type(data) == float):
                        data = round(data,3)
                        outputFile.write(" " + str(data) + "\t")
                else:
                        outputFile.write(" " + str(data) + "\t")
        outputFile.write("\n")
outputFile.write("\n")

#Write for Inter-BP data
outputFile.write("(C) Inter-BP\n")
outputFile.write("\t\t\t\t" + "Shift\t" + "Slide\t" + "Rise\t" + "Tilt\t" + "Roll\t" + "Twist\t" + "H-Ris\t" + "H-Twi\n")
#outputFile.write("\n")
for BP_data in Inter_BP:
	for data in BP_data:
		if (type(data) == int) or (type(data) == float):
			data = round(data,3)
                        if str(data)[0] == "-":
                                outputFile.write(str(data) + "\t")
                        else:
                                outputFile.write(" " + str(data) + "\t")
		else:
			outputFile.write(data + "\t")
	outputFile.write("\n")
outputFile.write("\n")

outputFile.write("Inter_BP Standard Deviations\n")
for std in Inter_BP_std_final:
        for data in std:
                if (type(data) == int) or (type(data) == float):
                        data = round(data,3)
                        outputFile.write(" " + str(data) + "\t")
                else:
                        outputFile.write(" " + str(data) + "\t")
        outputFile.write("\n")
outputFile.write("\n")

puckerColumns = ["C1'en","C1'ex","C2'en","C2'ex","C3'en","C3'ex","C4'en","C4'ex","O1'en","O1'ex"]

#Sort the puckered list into the order determined by the above puckerColumns list
puckerStrand1[0][1].sort() 
outputFile.write("(D) Backbone Parameters (Pucker Only)\n")
outputFile.write("Strand1\n")
counterStrand1 = 0
for puckerS1 in puckerStrand1:
	outputFile.write(str(counterStrand1 + 1) + ")\t" + puckerS1[0] + "\t")
	puckerS1[1].sort()
	columnCounter = 0
        i = 0
        while i <  len(puckerS1[1]):
                if puckerS1[1][i][:5] != puckerColumns[columnCounter]:
                        outputFile.write(puckerColumns[columnCounter] + " 0       ")
                        columnCounter += 1
                elif len(puckerS1[1][i]) == 7:
                        outputFile.write(puckerS1[1][i] + "    " + "   ")
                        columnCounter += 1
                        i += 1
                elif len(puckerS1[1][i]) == 8:
                        outputFile.write(puckerS1[1][i] + "   " + "   ")
                        columnCounter += 1
                        i += 1
                elif len(puckerS1[1][i]) == 9:
                        outputFile.write(puckerS1[1][i] + "  " + "   ")
                        columnCounter += 1
                        i += 1
                elif len(puckerS1[1][i]) == 10:
                        outputFile.write(puckerS1[1][i] + " " + "   ")
                        columnCounter += 1
                        i += 1
                else:
                        outputFile.write(puckerS1[1][i] + "   ")
                        columnCounter += 1
                        i += 1
		if (columnCounter == 9):
			found = False
			for conform in puckerS1[1]:
				if conform[:5] == "O1'ex":
					found = True
			if found == False:
				outputFile.write(puckerColumns[columnCounter] + " 0       ")
	counterStrand1 += 1
	outputFile.write("\n")

outputFile.write("\n")
outputFile.write("Strand2\n")
counterStrand2 = 0
for puckerS2 in puckerStrand2:
        outputFile.write(str(counterStrand2 + 1) + ")\t" + puckerS2[0] + "\t")
        puckerS2[1].sort()
	columnCounter = 0
	i = 0
        while i <  len(puckerS2[1]):
		if puckerS2[1][i][:5] != puckerColumns[columnCounter]:
			outputFile.write(puckerColumns[columnCounter] + " 0       ")
			columnCounter += 1
		elif len(puckerS2[1][i]) == 7:
                        outputFile.write(puckerS2[1][i] + "    " + "   ")
			columnCounter += 1
			i += 1
                elif len(puckerS2[1][i]) == 8:
                        outputFile.write(puckerS2[1][i] + "   " + "   ")
			columnCounter += 1
                        i += 1
                elif len(puckerS2[1][i]) == 9:
                        outputFile.write(puckerS2[1][i] + "  " + "   ")
			columnCounter += 1
                        i += 1
                elif len(puckerS2[1][i]) == 10:
                        outputFile.write(puckerS2[1][i] + " " + "   ")
			columnCounter += 1
                        i += 1
                else:
                        outputFile.write(puckerS2[1][i] + "   ")
			columnCounter += 1
                        i += 1
		if (columnCounter == 9):
                        found = False
                        for conform in puckerS2[1]:
                                if conform[:5] == "O1'ex":
                                        found = True
                        if found == False:
                                outputFile.write(puckerColumns[columnCounter] + " 0       ")

        counterStrand2 += 1
	outputFile.write("\n")

outputFile.close()
print("Done with first file")

##############################################################
	#Editting another file for Backbone Parameters
##############################################################
newFileName = ""
if analysisStartTime and (analysisEndTime and analysisEndTime != sys.maxint):
        newFile = open("backbone_parameters_DNA_for_" + sequenceName + "_from_" + str(analysisStartTime) + "_to_" + str(analysisEndTime) + ".txt", "w")
        newFileName = "backbone_parameters_DNA_for_" + sequenceName + "_from_" + str(analysisStartTime) + "_to_" + str(analysisEndTime) + ".txt"
	print("Second output file will be called: backbone_parameters_DNA_for_" + sequenceName + "_from_" + str(analysisStartTime) + "_to_" + str(analysisEndTime) + ".txt")
elif analysisStartTime and analysisEndTime == sys.maxint:
        newFile = open("backbone_parameters_DNA_for_" + sequenceName + "_from_" + str(analysisStartTime) + "_to_end" + ".txt", "w")
	newFileName = "backbone_parameters_DNA_for_" + sequenceName + "_from_" + str(analysisStartTime) + "_to_end" + ".txt"
        print("Second output file will be called: backbone_parameters_DNA_for_" + sequenceName + "_from_" + str(analysisStartTime) + "_to_end" + ".txt")
elif not analysisStartTime and analysisEndTime:
        newFile = open("backbone_parameters_DNA_for_" + sequenceName + "_from_5350_to_" + str(analysisEndTime) + ".txt", "w")
        newFileName = "backbone_parameters_DNA_for_" + sequenceName + "_from_5350_to_" + str(analysisEndTime) + ".txt"
	print("Second output file will be called: backbone_parameters_DNA_for_" + sequenceName + "_from_5350_to_" + str(analysisEndTime) + ".txt")
else:
        newFile = open("backbone_parameters_DNA_for_" + sequenceName + "_from_5350_to_end.txt", "w")
	newFileName = "backbone_parameters_DNA_for_" + sequenceName + "_from_5350_to_end.txt"
        print("Second output file will be called: backbone_parameters_DNA_for_" + sequenceName + "_from_5350_to_end.txt")

print("Working on second output file")

newFile.write("\t\t\tThis file contains the remainder of the Backbone Parameter data\n\n")
newFile.write("Strand 1\n\n")
index1 = 0
timeDif = int(float(startTime[1])) - int(float(startTime[0]))
while index1 < len(Backbone_params_strand1):
	nucleotideData = Backbone_params_strand1[index1]
	newFile.write("Nucleotide:\t" + nucleotideData[1] + "\n")
	newFile.write("Time    "+ "Alpha     " + "Beta      " + "Gamma     " + "Delta     " + "Epsilon   " + "Zeta      " + "Chi       " + " Phase    " + "Amplitude ")
	newFile.write("\n")

	#Gathering Data	
	Alpha = nucleotideData[3].split()
	Beta = nucleotideData[4].split()
	Gamma = nucleotideData[5].split()
	Delta = nucleotideData[6].split()
	Epsil = nucleotideData[7].split()
	Zeta = nucleotideData[8].split()
	Chi = nucleotideData[9].split()	
	Phase = nucleotideData[10].split()
	Amp = nucleotideData[11].split()	

	time = int(float(startTime[0])) 			#Start time of the simulation
	subIndex = 0
	while subIndex < numFiles:
		#Adding the time stamp :-)
		if len(str(time)) == 4:
			newFile.write(str(time) + "    ")
		elif len(str(time)) == 5:
			newFile.write(str(time) + "   ")
		elif len(str(time)) == 6:
			newFile.write(str(time) + "  ")
		elif len(str(time)) == 7:
			newFile.write(str(time) + " ")
		elif len(str(time)) == 8:
			newFile.write(str(time))
		
		#Adding components
		if len(Alpha[subIndex]) == 3:
			newFile.write(Alpha[subIndex] + "       ")
		elif len(Alpha[subIndex]) == 4:
			newFile.write(Alpha[subIndex] + "      ")
		elif len(Alpha[subIndex]) == 5:
			newFile.write(Alpha[subIndex] + "     ")
		elif len(Alpha[subIndex]) == 6:
			newFile.write(Alpha[subIndex] + "    ")

                if len(Beta[subIndex]) == 3:
                        newFile.write(Beta[subIndex] + "       ")
                elif len(Beta[subIndex]) == 4:
                        newFile.write(Beta[subIndex] + "      ")
                elif len(Beta[subIndex]) == 5:
                        newFile.write(Beta[subIndex] + "     ")
                elif len(Beta[subIndex]) == 6:
                        newFile.write(Beta[subIndex] + "    ")

                if len(Gamma[subIndex]) == 3:
                        newFile.write(Gamma[subIndex] + "       ")
                elif len(Gamma[subIndex]) == 4:
                        newFile.write(Gamma[subIndex] + "      ")
                elif len(Gamma[subIndex]) == 5:
                        newFile.write(Gamma[subIndex] + "     ")
                elif len(Gamma[subIndex]) == 6:
                        newFile.write(Gamma[subIndex] + "    ")

                if len(Delta[subIndex]) == 3:
                        newFile.write(Delta[subIndex] + "       ")
                elif len(Delta[subIndex]) == 4:
                        newFile.write(Delta[subIndex] + "      ")
                elif len(Delta[subIndex]) == 5:
                        newFile.write(Delta[subIndex] + "     ")
                elif len(Delta[subIndex]) == 6:
                        newFile.write(Delta[subIndex] + "    ")

                if len(Epsil[subIndex]) == 3:
                        newFile.write(Epsil[subIndex] + "       ")
                elif len(Epsil[subIndex]) == 4:
                        newFile.write(Epsil[subIndex] + "      ")
                elif len(Epsil[subIndex]) == 5:
                        newFile.write(Epsil[subIndex] + "     ")
                elif len(Epsil[subIndex]) == 6:
                        newFile.write(Epsil[subIndex] + "    ")
            
	   	if len(Zeta[subIndex]) == 3:
                        newFile.write(Zeta[subIndex] + "       ")
                elif len(Zeta[subIndex]) == 4:
                        newFile.write(Zeta[subIndex] + "      ")
                elif len(Zeta[subIndex]) == 5:
                        newFile.write(Zeta[subIndex] + "     ")
                elif len(Zeta[subIndex]) == 6:
                        newFile.write(Zeta[subIndex] + "    ")

		if len(Chi[subIndex]) == 3:
                        newFile.write(Chi[subIndex] + "       ")
                elif len(Chi[subIndex]) == 4:
                        newFile.write(Chi[subIndex] + "      ")
                elif len(Chi[subIndex]) == 5:
                        newFile.write(Chi[subIndex] + "     ")
                elif len(Chi[subIndex]) == 6:
                        newFile.write(Chi[subIndex] + "    ")

		if len(Phase[subIndex]) == 3:
                        newFile.write(Phase[subIndex] + "       ")
                elif len(Phase[subIndex]) == 4:
                        newFile.write(Phase[subIndex] + "      ")
                elif len(Phase[subIndex]) == 5:
                        newFile.write(Phase[subIndex] + "     ")
                elif len(Phase[subIndex]) == 6:
                        newFile.write(Phase[subIndex] + "    ")

                if len(Amp[subIndex]) == 3:
                        newFile.write(Amp[subIndex] + "       " + "\n")
                elif len(Amp[subIndex]) == 4:
                        newFile.write(Amp[subIndex] + "      " + "\n")
                elif len(Amp[subIndex]) == 5:
                        newFile.write(Amp[subIndex] + "     " + "\n")
                elif len(Amp[subIndex]) == 6:
                        newFile.write(Amp[subIndex] + "    " + "\n")
		
		subIndex += 1
		time += timeDif
	index1 += 1
	newFile.write("\n\n")

newFile.write("Strand 2\n\n")
index1 = 0
timeDif = int(float(startTime[1])) - int(float(startTime[0]))
while index1 < len(Backbone_params_strand2):
        nucleotideData = Backbone_params_strand2[index1]
        newFile.write("Nucleotide:\t" + nucleotideData[1] + "\n")
        newFile.write("Time    "+ "Alpha     " + "Beta      " + "Gamma     " + "Delta     " + "Epsil     " + "Zeta      " + "Chi       " + " Phase    " + "Amplitude ")
	newFile.write("\n")

        #Gathering Data
        Alpha = nucleotideData[3].split()
        Beta = nucleotideData[4].split()
        Gamma = nucleotideData[5].split()
        Delta = nucleotideData[6].split()
        Epsil = nucleotideData[7].split()
        Zeta = nucleotideData[8].split()
        Chi = nucleotideData[9].split()
	Phase = nucleotideData[10].split()
        Amp = nucleotideData[11].split()

        time = int(float(startTime[0]))                     #Start time of the simulationi
        subIndex = 0
        while subIndex < numFiles:
                #Adding the time stamp :-)
                if len(str(time)) == 4:
                        newFile.write(str(time) + "    ")
                elif len(str(time)) == 5:
                        newFile.write(str(time) + "   ")
                elif len(str(time)) == 6:
                        newFile.write(str(time) + "  ")
                elif len(str(time)) == 7:
                        newFile.write(str(time) + " ")
                elif len(str(time)) == 8:
                        newFile.write(str(time))

                #Adding components
                if len(Alpha[subIndex]) == 3:
                        newFile.write(Alpha[subIndex] + "       ")
                elif len(Alpha[subIndex]) == 4:
                        newFile.write(Alpha[subIndex] + "      ")
                elif len(Alpha[subIndex]) == 5:
                        newFile.write(Alpha[subIndex] + "     ")
                elif len(Alpha[subIndex]) == 6:
                        newFile.write(Alpha[subIndex] + "    ")

                if len(Beta[subIndex]) == 3:
                        newFile.write(Beta[subIndex] + "       ")
                elif len(Beta[subIndex]) == 4:
                        newFile.write(Beta[subIndex] + "      ")
                elif len(Beta[subIndex]) == 5:
                        newFile.write(Beta[subIndex] + "     ")
                elif len(Beta[subIndex]) == 6:
                        newFile.write(Beta[subIndex] + "    ")

                if len(Gamma[subIndex]) == 3:
                        newFile.write(Gamma[subIndex] + "       ")
                elif len(Gamma[subIndex]) == 4:
                        newFile.write(Gamma[subIndex] + "      ")
                elif len(Gamma[subIndex]) == 5:
                        newFile.write(Gamma[subIndex] + "     ")
                elif len(Gamma[subIndex]) == 6:
                        newFile.write(Gamma[subIndex] + "    ")
		
		if len(Delta[subIndex]) == 3:
                        newFile.write(Delta[subIndex] + "       ")
                elif len(Delta[subIndex]) == 4:
                        newFile.write(Delta[subIndex] + "      ")
                elif len(Delta[subIndex]) == 5:
                        newFile.write(Delta[subIndex] + "     ")
                elif len(Delta[subIndex]) == 6:
                        newFile.write(Delta[subIndex] + "    ")

                if len(Epsil[subIndex]) == 3:
                        newFile.write(Epsil[subIndex] + "       ")
                elif len(Epsil[subIndex]) == 4:
                        newFile.write(Epsil[subIndex] + "      ")
                elif len(Epsil[subIndex]) == 5:
                        newFile.write(Epsil[subIndex] + "     ")
                elif len(Epsil[subIndex]) == 6:
                        newFile.write(Epsil[subIndex] + "    ")

                if len(Zeta[subIndex]) == 3:
                        newFile.write(Zeta[subIndex] + "       ")
                elif len(Zeta[subIndex]) == 4:
                        newFile.write(Zeta[subIndex] + "      ")
                elif len(Zeta[subIndex]) == 5:
                        newFile.write(Zeta[subIndex] + "     ")
                elif len(Zeta[subIndex]) == 6:
                        newFile.write(Zeta[subIndex] + "    ")

		if len(Zeta[subIndex]) == 3:
                        newFile.write(Chi[subIndex] + "       ")
                elif len(Zeta[subIndex]) == 4:
                        newFile.write(Chi[subIndex] + "      ")
                elif len(Zeta[subIndex]) == 5:
                        newFile.write(Chi[subIndex] + "     ")
                elif len(Zeta[subIndex]) == 6:
                        newFile.write(Chi[subIndex] + "    ")

                if len(Zeta[subIndex]) == 3:
                        newFile.write(Phase[subIndex] + "       ")
                elif len(Zeta[subIndex]) == 4:
                        newFile.write(Phase[subIndex] + "      ")
                elif len(Zeta[subIndex]) == 5:
                        newFile.write(Phase[subIndex] + "     ")
                elif len(Zeta[subIndex]) == 6:
                        newFile.write(Phase[subIndex] + "    ")

                if len(Chi[subIndex]) == 3:
                        newFile.write(Amp[subIndex] + "       " + "\n")
                elif len(Chi[subIndex]) == 4:
                        newFile.write(Amp[subIndex] + "      " + "\n")
                elif len(Chi[subIndex]) == 5:
                        newFile.write(Amp[subIndex] + "     " + "\n")
                elif len(Chi[subIndex]) == 6:
                        newFile.write(Amp[subIndex] + "    " + "\n")

                subIndex += 1
                time += timeDif
        index1 += 1
        newFile.write("\n")
newFile.close()

###############################################################################
	#New file for frequency histogram of helical parameters
###############################################################################
helic_freqFileName = ""
if analysisStartTime and (analysisEndTime and analysisEndTime != sys.maxint):
        helic_freqFile = open("helical_param_freq_hist_DNA_for_" + sequenceName + "_from_" + str(analysisStartTime) + "_to_" + str(analysisEndTime) + ".txt", "w")
        helic_freqFileName = "helical_param_freq_hist_DNA_for_" + sequenceName + "_from_" + str(analysisStartTime) + "_to_" + str(analysisEndTime) + ".txt"
	print("Final output file will be called: helical_param_freq_hist_DNA_for_" + sequenceName + "_from_" + str(analysisStartTime) + "_to_" + str(analysisEndTime) + ".txt")
elif analysisStartTime and analysisEndTime == sys.maxint:
        helic_freqFile = open("helical_param_freq_hist_DNA_for_" + sequenceName + "_from_" + str(analysisStartTime) + "_to_end" + ".txt", "w")
	helic_freqFileName = "helical_param_freq_hist_DNA_for_" + sequenceName + "_from_" + str(analysisStartTime) + "_to_end" + ".txt"
        print("Final output file will be called: helical_param_freq_hist_DNA_for_" + sequenceName + "_from_" + str(analysisStartTime) + "_to_end" + ".txt")
elif not analysisStartTime and analysisEndTime:
        helic_freqFile = open("helical_param_freq_hist_DNA_for_" + sequenceName + "_from_5350_to_" + str(analysisEndTime) + ".txt", "w")
	helic_freqFileName = "helical_param_freq_hist_DNA_for_" + sequenceName + "_from_5350_to_" + str(analysisEndTime) + ".txt"
        print("Final output file will be called: helical_param_freq_hist_DNA_for_" + sequenceName + "_from_5350_to_" + str(analysisEndTime) + ".txt")
else:
        helic_freqFile = open("helical_param_freq_hist_DNA_for_" + sequenceName + "_from_5350_to_end.txt", "w")
	helic_freqFile = "helical_param_freq_hist_DNA_for_" + sequenceName + "_from_5350_to_end.txt"
        print("Final output file will be called: helical_param_freq_hist_DNA_for_" + sequenceName + "_from_0_to_end.txt")

#This portion of the file uses to compute frequencies: 
	#BP_Axis_std
	#Intra_BP_std
	#Inter_BP_std
incline_bin = []
tip_bin = []
buckle_bin = []
propel_bin = []
opening_bin = []
tilt_bin = []
roll_bin = []
twist_bin = []
for numBase in range(len(sequence)):
	#Each base pair gets its own set of bins in the overal bin list, indiciated by numBase
	incline_bin.append([])
	tip_bin.append([])
	buckle_bin.append([])
	propel_bin.append([])
	opening_bin.append([])
	tilt_bin.append([])
	roll_bin.append([])
	twist_bin.append([])

	for i in range(0,73):			#Create the bins for the helical rotation parameters
		incline_bin[numBase].append(0)
		tip_bin[numBase].append(0)
		buckle_bin[numBase].append(0)
		propel_bin[numBase].append(0)
		opening_bin[numBase].append(0)
		tilt_bin[numBase].append(0)
		roll_bin[numBase].append(0)
		twist_bin[numBase].append(0)
	indexBP = 6
	while indexBP < len(BP_Axis_std[numBase]):
		for num in BP_Axis_std[numBase][indexBP]:
			try:
				temp = float(num)
				temp += 180
				tempIndex = int(temp/5) #5 for the size of each bin
				if indexBP == 6:
					incline_bin[numBase][tempIndex] += 1
				if indexBP == 7:
					tip_bin[numBase][tempIndex] += 1
			except ValueError:
				break
		indexBP += 1

	indexIntra = 7
	while indexIntra < len(Intra_BP_std[numBase]):
		for num in Intra_BP_std[numBase][indexIntra]:
			try:
				temp = float(num)
				temp += 180
				tempIndex = int(temp/5)
				if indexIntra == 7:
					buckle_bin[numBase][tempIndex] += 1
				elif indexIntra == 8:
					propel_bin[numBase][tempIndex] += 1
				elif indexIntra == 9:
					opening_bin[numBase][tempIndex] += 1
			except ValueError:
				break
		indexIntra += 1

	indexInter = 7
	while (numBase != (len(sequence)-1)) and (indexInter < len(Inter_BP_std[numBase])):
		for num in Inter_BP_std[numBase][indexInter]:
			try:
				temp = float(num)
				temp += 180
				tempIndex = int(temp/5)
				if indexInter == 7:
					tilt_bin[numBase][tempIndex] += 1
				elif indexInter == 8:
					roll_bin[numBase][tempIndex] += 1
				elif indexInter == 9:
					twist_bin[numBase][tempIndex] += 1
			except ValueError:
				break
		indexInter += 1

#Write to the new file
helic_freqFile.write("\t\tThis file contains the frequencies of the helical parameters\n\n")
helic_freqFile.write("Number of files analyzed:\t" + str(numFiles)+"\n")
helic_freqFile.write("The sequence analyzed:\t")
for base in sequence:
	helic_freqFile.write(base)
helic_freqFile.write("\n")
#print incline_bin
for i in range(len(sequence)):
	helic_freqFile.write("Nucleotide:\t" + sequence[i]+"\n")
	helic_freqFile.write("Angle   " + "Incline   "+ "Tip     " + "Buckle  " + "Propel    " + "Opening   " + "Tilt      " + "Roll      " + "Twist     ")
	helic_freqFile.write("\n")
	bin_index = 0
	while bin_index < 72:
	        if len(str(bin_index*5 - 180)) == 6:
                        helic_freqFile.write(str(bin_index*5 - 180) + "  ")
                elif len(str(bin_index*5 - 180)) == 5:
                        helic_freqFile.write(str(bin_index*5 - 180) + "   ")
                elif len(str(bin_index*5 - 180)) == 4:
                        helic_freqFile.write(str(bin_index*5 - 180) + "    ")
                elif len(str(bin_index*5 - 180)) == 3:
                        helic_freqFile.write(str(bin_index*5 - 180) + "     ")
                elif len(str(bin_index*5 - 180)) == 2:
                        helic_freqFile.write(str(bin_index*5 - 180) + "      ")
                elif len(str(bin_index*5 - 180)) == 1:
                        helic_freqFile.write(str(bin_index*5 - 180) + "       ")

		len_incline = 10 - len(str(incline_bin[i][bin_index]))
		len_tip = 10 - len(str(tip_bin[i][bin_index]))
		len_buckle = 10 - len(str(buckle_bin[i][bin_index]))
		len_propel = 10 - len(str(propel_bin[i][bin_index]))
		len_opening = 10 - len(str(opening_bin[i][bin_index]))
		len_tilt = 10 - len(str(tilt_bin[i][bin_index]))
		len_roll = 10 - len(str(roll_bin[i][bin_index]))
		len_twist = 10 - len(str(twist_bin[i][bin_index]))

		space_incline = ""
		for space in range(len_incline):
			space_incline += " "
		helic_freqFile.write(str(incline_bin[i][bin_index]) + space_incline)

                space_tip = ""
                for space in range(len_tip):
                        space_tip += " "
                helic_freqFile.write(str(tip_bin[i][bin_index]) + space_tip)

                space_buckle = ""
                for space in range(len_buckle):
                        space_buckle += " "
                helic_freqFile.write(str(buckle_bin[i][bin_index]) + space_buckle)

                space_propel = ""
                for space in range(len_propel):
                        space_propel += " "
                helic_freqFile.write(str(propel_bin[i][bin_index]) + space_propel)

                space_opening = ""
                for space in range(len_opening):
                        space_opening += " "
                helic_freqFile.write(str(opening_bin[i][bin_index]) + space_opening)

                space_tilt = ""
                for space in range(len_tilt):
                        space_tilt += " "
                helic_freqFile.write(str(tilt_bin[i][bin_index]) + space_tilt)

                space_roll = ""
                for space in range(len_roll):
                        space_roll += " "
                helic_freqFile.write(str(roll_bin[i][bin_index]) + space_roll)

                space_twist = ""
                for space in range(len_twist):
                        space_twist += " "
                helic_freqFile.write(str(twist_bin[i][bin_index]) + space_twist + "\n")
		bin_index += 1
	helic_freqFile.write("\n")
helic_freqFile.close()	

###############################################################################
	#New file for frequency histograms for backbone angles
###############################################################################
freqFileName = ""
if analysisStartTime and (analysisEndTime and analysisEndTime != sys.maxint):
        freqFile = open("backbone_freq_hist_DNA_for_" + sequenceName + "_from_" + str(analysisStartTime) + "_to_" + str(analysisEndTime) + ".txt", "w")
        freqFileName = "backbone_freq_hist_DNA_for_" + sequenceName + "_from_" + str(analysisStartTime) + "_to_" + str(analysisEndTime) + ".txt"
	print("Final output file will be called: backbone_freq_hist_DNA_for_" + sequenceName + "_from_" + str(analysisStartTime) + "_to_" + str(analysisEndTime) + ".txt")
elif analysisStartTime and analysisEndTime == sys.maxint:
        freqFile = open("backbone_freq_hist_DNA_for_" + sequenceName + "_from_" + str(analysisStartTime) + "_to_end" + ".txt", "w")
	freqFileName = "backbone_freq_hist_DNA_for_" + sequenceName + "_from_" + str(analysisStartTime) + "_to_end" + ".txt"
        print("Final output file will be called: backbone_freq_hist_DNA_for_" + sequenceName + "_from_" + str(analysisStartTime) + "_to_end" + ".txt")
elif not analysisStartTime and analysisEndTime:
        freqFile = open("backbone_freq_hist_DNA_for_" + sequenceName + "_from_5350_to_" + str(analysisEndTime) + ".txt", "w")
	freqFileName = "backbone_freq_hist_DNA_for_" + sequenceName + "_from_5350_to_" + str(analysisEndTime) + ".txt"
        print("Final output file will be called: backbone_freq_hist_DNA_for_" + sequenceName + "_from_5350_to_" + str(analysisEndTime) + ".txt")
else:
        freqFile = open("backbone_freq_hist_DNA_for_" + sequenceName + "_from_5350_to_end.txt", "w")
	freqFileName = "backbone_freq_hist_DNA_for_" + sequenceName + "_from_5350_to_end.txt"
        print("Final output file will be called: backbone_freq_hist_DNA_for_" + sequenceName + "_from_0_to_end.txt")

freqFile.write("\t\t\tThis file contains the the frequencies of the backbone parameters\n\n")
freqFile.write("Strand 1\n\n")
freqFile.write("Number of files analyzed:\t " + str(numFiles) + "\n")
index1 = 0
print("Calculating Frequencies")
while index1 < len(Backbone_params_strand1):
        nucleotideData = Backbone_params_strand1[index1]
        freqFile.write("Nucleotide:\t" + nucleotideData[1] + "\n")
        freqFile.write("Angle   "+ "Alpha     " + "Beta      " + "Gamma     " + "Delta     " + "Epsilon   " + "Zeta      " + "Chi       " + " Phase    ")
        freqFile.write("\n")
	
	alpha_bins = []
	beta_bins = []
	gamma_bins = []
	delta_bins = []
	epsil_bins = []
	zeta_bins = []
	chi_bins = []
	phase_bins = []
	#make all the bins (5 angles per bin --> 72 bins total)
	for bin in range(0,73):
		alpha_bins.append(0)
		beta_bins.append(0)
		gamma_bins.append(0)
		delta_bins.append(0)
		epsil_bins.append(0)
		zeta_bins.append(0)
		chi_bins.append(0)
		phase_bins.append(0)

        Alpha = nucleotideData[3].split()
        Beta = nucleotideData[4].split()
        Gamma = nucleotideData[5].split()
        Delta = nucleotideData[6].split()
        Epsil = nucleotideData[7].split()
        Zeta = nucleotideData[8].split()
        Chi = nucleotideData[9].split()
        Phase = nucleotideData[10].split()
	
	for element in Alpha:
		try:
			temp = float(element)
			temp += 180
 			tempIndex = int(temp/5) #5 for the size of each bin
			alpha_bins[tempIndex] += 1 
		except ValueError:
			break

        for element in Beta:
                try:
                        temp = float(element)
                        temp += 180
                        tempIndex = int(temp/5) #5 for the size of each bin
                        beta_bins[tempIndex] += 1
                except ValueError:
                        break

        for element in Gamma:
                try:
                        temp = float(element)
                        temp += 180
                        tempIndex = int(temp/5) #5 for the size of each bin
			#print(temp, tempIndex)
                        gamma_bins[tempIndex] += 1
                except ValueError:
                        break

        for element in Delta:
                try:
                        temp = float(element)
                        temp += 180
                        tempIndex = int(temp/5) #5 for the size of each bin
                        delta_bins[tempIndex] += 1
                except ValueError:
                        break

        for element in Epsil:
                try:
                        temp = float(element)
                        temp += 180
                        tempIndex = int(temp/5) #5 for the size of each bin
                        epsil_bins[tempIndex] += 1
                except ValueError:
                        break

        for element in Zeta:
                try:
                        temp = float(element)
                        temp += 180
                        tempIndex = int(temp/5) #5 for the size of each bin
                        zeta_bins[tempIndex] += 1
                except ValueError:
                        break

        for element in Chi:
                try:
                        temp = float(element)
                        temp += 180
                        tempIndex =int(temp/5) #5 for the size of each bin
                        chi_bins[tempIndex] += 1
                except ValueError:
                        break

        for element in Phase:
                try:
                        temp = float(element)
                        temp += 180
                        tempIndex = int(temp/5) #5 for the size of each bin
                        phase_bins[tempIndex] += 1
                except ValueError:
                        break

	bin_index = 0
	while bin_index < 72:
		if len(str(bin_index*5 - 180)) == 6:
			freqFile.write(str(bin_index*5 - 180) + "  ")
                elif len(str(bin_index*5 - 180)) == 5:
                        freqFile.write(str(bin_index*5 - 180) + "   ")
                elif len(str(bin_index*5 - 180)) == 4:
                        freqFile.write(str(bin_index*5 - 180) + "    ")
                elif len(str(bin_index*5 - 180)) == 3:
                        freqFile.write(str(bin_index*5 - 180) + "     ")
                elif len(str(bin_index*5 - 180)) == 2:
                        freqFile.write(str(bin_index*5 - 180) + "      ")
                elif len(str(bin_index*5 - 180)) == 1:
                        freqFile.write(str(bin_index*5 - 180) + "       ")


		len_alpha = 10 - len(str(alpha_bins[bin_index]))
                len_beta = 10 - len(str(beta_bins[bin_index]))
                len_gamma = 10 - len(str(gamma_bins[bin_index]))
                len_delta = 10 - len(str(delta_bins[bin_index]))
                len_epsil = 10 - len(str(epsil_bins[bin_index]))
                len_zeta = 10 - len(str(zeta_bins[bin_index]))
                len_chi = 10 - len(str(chi_bins[bin_index]))
                len_phase = 10 - len(str(phase_bins[bin_index]))

		space_alpha = ""
		for space in range(len_alpha):
			space_alpha += " "
		freqFile.write(str(alpha_bins[bin_index]) + space_alpha)		

                space_beta = ""
                for space in range(len_beta):
                        space_beta += " "
		freqFile.write(str(beta_bins[bin_index]) + space_beta)

                space_gamma = ""
                for space in range(len_gamma):
                        space_gamma += " "
		freqFile.write(str(gamma_bins[bin_index]) + space_gamma)

                space_delta = ""
                for space in range(len_delta):
                        space_delta += " "
		freqFile.write(str(delta_bins[bin_index]) + space_delta)

                space_epsil = ""
                for space in range(len_epsil):
                        space_epsil += " "
		freqFile.write(str(epsil_bins[bin_index]) + space_epsil)

                space_zeta = ""
                for space in range(len_zeta):
                        space_zeta += " "
		freqFile.write(str(zeta_bins[bin_index]) + space_zeta)

                space_chi = ""
                for space in range(len_chi):
                        space_chi += " "
		freqFile.write(str(chi_bins[bin_index]) + space_chi)
                
		space_phase = ""
                for space in range(len_phase):
                        space_phase += " "
		freqFile.write(str(phase_bins[bin_index]) + space_phase + "\n")
		bin_index += 1
	index1 += 1 
	freqFile.write("\n")


freqFile.write("Strand 2\n\n")
index2 = 0
while index2 < len(Backbone_params_strand2):
        nucleotideData = Backbone_params_strand2[index2]
        freqFile.write("Nucleotide:\t" + nucleotideData[1] + "\n")
        freqFile.write("Angle   "+ "Alpha     " + "Beta      " + "Gamma     " + "Delta     " + "Epsilon   " + "Zeta      " + "Chi       " + " Phase    ")
        freqFile.write("\n")
	
	alpha_bins = []
	beta_bins = []
	gamma_bins = []
	delta_bins = []
	epsil_bins = []
	zeta_bins = []
	chi_bins = []
	phase_bins = []
	#make all the bins (5 angles per bin --> 72 bins total)
	for bin in range(0,73):
		alpha_bins.append(0)
		beta_bins.append(0)
		gamma_bins.append(0)
		delta_bins.append(0)
		epsil_bins.append(0)
		zeta_bins.append(0)
		chi_bins.append(0)
		phase_bins.append(0)

        Alpha = nucleotideData[3].split()
        Beta = nucleotideData[4].split()
        Gamma = nucleotideData[5].split()
        Delta = nucleotideData[6].split()
        Epsil = nucleotideData[7].split()
        Zeta = nucleotideData[8].split()
        Chi = nucleotideData[9].split()
        Phase = nucleotideData[10].split()
	
	for element in Alpha:
		try:
			temp = float(element)
			temp += 180
 			tempIndex = int(temp/5) #5 for the size of each bin
			alpha_bins[tempIndex] += 1 
		except ValueError:
			break

        for element in Beta:
                try:
                        temp = float(element)
                        temp += 180
                        tempIndex = int(temp/5) #5 for the size of each bin
                        beta_bins[tempIndex] += 1
                except ValueError:
                        break

        for element in Gamma:
                try:
                        temp = float(element)
                        temp += 180
                        tempIndex = int(temp/5) #5 for the size of each bin
			#print(temp, tempIndex)
                        gamma_bins[tempIndex] += 1
                except ValueError:
                        break

        for element in Delta:
                try:
                        temp = float(element)
                        temp += 180
                        tempIndex = int(temp/5) #5 for the size of each bin
                        delta_bins[tempIndex] += 1
                except ValueError:
                        break

        for element in Epsil:
                try:
                        temp = float(element)
                        temp += 180
                        tempIndex = int(temp/5) #5 for the size of each bin
                        epsil_bins[tempIndex] += 1
                except ValueError:
                        break

        for element in Zeta:
                try:
                        temp = float(element)
                        temp += 180
                        tempIndex = int(temp/5) #5 for the size of each bin
                        zeta_bins[tempIndex] += 1
                except ValueError:
                        break

        for element in Chi:
                try:
                        temp = float(element)
                        temp += 180
                        tempIndex =int(temp/5) #5 for the size of each bin
                        chi_bins[tempIndex] += 1
                except ValueError:
                        break

        for element in Phase:
                try:
                        temp = float(element)
                        temp += 180
                        tempIndex = int(temp/5) #5 for the size of each bin
                        phase_bins[tempIndex] += 1
                except ValueError:
                        break

	bin_index = 0
	while bin_index < 72:
		if len(str(bin_index*5 - 180)) == 6:
			freqFile.write(str(bin_index*5 - 180) + "  ")
                elif len(str(bin_index*5 - 180)) == 5:
                        freqFile.write(str(bin_index*5 - 180) + "   ")
                elif len(str(bin_index*5 - 180)) == 4:
                        freqFile.write(str(bin_index*5 - 180) + "    ")
                elif len(str(bin_index*5 - 180)) == 3:
                        freqFile.write(str(bin_index*5 - 180) + "     ")
                elif len(str(bin_index*5 - 180)) == 2:
                        freqFile.write(str(bin_index*5 - 180) + "      ")
                elif len(str(bin_index*5 - 180)) == 1:
                        freqFile.write(str(bin_index*5 - 180) + "       ")


		len_alpha = 10 - len(str(alpha_bins[bin_index]))
                len_beta = 10 - len(str(beta_bins[bin_index]))
                len_gamma = 10 - len(str(gamma_bins[bin_index]))
                len_delta = 10 - len(str(delta_bins[bin_index]))
                len_epsil = 10 - len(str(epsil_bins[bin_index]))
                len_zeta = 10 - len(str(zeta_bins[bin_index]))
                len_chi = 10 - len(str(chi_bins[bin_index]))
                len_phase = 10 - len(str(phase_bins[bin_index]))

		space_alpha = ""
		for space in range(len_alpha):
			space_alpha += " "
		freqFile.write(str(alpha_bins[bin_index]) + space_alpha)		

                space_beta = ""
                for space in range(len_beta):
                        space_beta += " "
		freqFile.write(str(beta_bins[bin_index]) + space_beta)

                space_gamma = ""
                for space in range(len_gamma):
                        space_gamma += " "
		freqFile.write(str(gamma_bins[bin_index]) + space_gamma)

                space_delta = ""
                for space in range(len_delta):
                        space_delta += " "
		freqFile.write(str(delta_bins[bin_index]) + space_delta)

                space_epsil = ""
                for space in range(len_epsil):
                        space_epsil += " "
		freqFile.write(str(epsil_bins[bin_index]) + space_epsil)

                space_zeta = ""
                for space in range(len_zeta):
                        space_zeta += " "
		freqFile.write(str(zeta_bins[bin_index]) + space_zeta)

                space_chi = ""
                for space in range(len_chi):
                        space_chi += " "
                freqFile.write(str(chi_bins[bin_index]) + space_chi)

                space_phase = ""
                for space in range(len_phase):
                        space_phase += " "
                freqFile.write(str(phase_bins[bin_index]) + space_phase + "\n")
                bin_index += 1
	index2 += 1 
	freqFile.write("\n")


freqFile.close()


########################################
#Handles the evaulation of individual nucleotides given by the users
########################################
print("\n\n")

if analyze_BP == True:
	print("Creating output file for individual nucleotides")
	nucleotides_to_analyze = []
	x = int(raw_input("how many bases would you like to analyze? "))
	baseNumber = 0
	while baseNumber < x:
		base = raw_input("Enter base # (from pdb) of the nucleotide you would like to analyze: ")
		nucleotides_to_analyze.append(int(base))
		baseNumber += 1
	nucleotides_to_analyze.sort() 			#Sort so we can easily traverse the DS
	#Open up each complete output file and copy/paste relevant information to a new output file
	#Files needed:
	#outputFileName --> curve+ helical parameters file
	#newFileName --> raw data for freqFileName
	#helic_freqFileName --> rotational data in hist form of outputFileName
	#freqFileName --> data from newFileName displayed in hist form
	
	#Start with outPutFileName
	baseIndex = 0
	paramHelicalFile = open(outputFileName,"r")
	newFileNameForIndNucs = "individual_curves_data_for_"
	baseIndex1 = 0
	indSequence = ""
	while baseIndex1 < len(sequence):
		if baseIndex1 + 1 in nucleotides_to_analyze:
			newFileNameForIndNucs += sequence[baseIndex1]
			indSequence += sequence[baseIndex1]
		else:
			newFileNameForIndNucs += sequence[baseIndex1].lower()
			indSequence += sequence[baseIndex1].lower()
		baseIndex1 += 1
	newFileNameForIndNucs += "_all_in_one.txt"
	
	indFile = open(newFileNameForIndNucs,"w")
	indFile.write("\t\tIndividual Data for " + indSequence + "\n")

	#Starting analysis on helical parameters -- means and standard deviations
	line = paramHelicalFile.readline()
	while line:
		line_list = line.split()
		line_list.append(False)
		if line_list[0] == "(A)":
			#print "Made it to AAAAA"
			indFile.write(line)
			line = paramHelicalFile.readline()
			indFile.write(line)
			i = 0
			while i < len(nucleotides_to_analyze):
				element = str(nucleotides_to_analyze[i]) + ")"
				line = paramHelicalFile.readline()
				line_list = line.split()
				line_list.append(False)
				while line_list[0] != element:
					#REID FIXED THIS SHIT
					line = paramHelicalFile.readline()
					line_list = line.split()
					line_list.append(False)
				indFile.write(line)
				i += 1
			line = paramHelicalFile.readline()
			line_list = line.split()
			line_list.append(False)
			indFile.write("\n")
			while line_list[0] != "BP_Axis":
				line = paramHelicalFile.readline()
				line_list = line.split()
				line_list.append(False)
			indFile.write(line.lstrip())
			i = 0
			while i < len(nucleotides_to_analyze):
				element = str(nucleotides_to_analyze[i]) + ")"
                                line = paramHelicalFile.readline()
                                line_list = line.split()
				line_list.append(False)
                                while line_list[0] != element:
                                        line = paramHelicalFile.readline()
					line_list = line.split()
					line_list.append(False)
                                indFile.write(line.lstrip())
                                i += 1
			indFile.write("\n")
			#print "made it out of A"
			line = paramHelicalFile.readline()
		elif line_list[0] == "(B)":
			#print "B made it"
                        indFile.write(line.lstrip())
                        line = paramHelicalFile.readline()
                        indFile.write(line)
                        i = 0
                        while i < len(nucleotides_to_analyze):
                                element = str(nucleotides_to_analyze[i]) + ")"
                                line = paramHelicalFile.readline()
                                line_list = line.split()
				line_list.append(False)
                                while line_list[0] != element:
                                        line = paramHelicalFile.readline()
					line_list = line.split()
					line_list.append(False)
                                indFile.write(line.lstrip())
                                i += 1
                        line = paramHelicalFile.readline()
                        line_list = line.split()
			line_list.append(False)
			indFile.write("\n")
                        while line_list[0] != "Intra_BP":
				line = paramHelicalFile.readline()
                                line_list = line.split()
				line_list.append(False)
                        indFile.write(line)
                        i = 0
                        while i < len(nucleotides_to_analyze):
                                element = str(nucleotides_to_analyze[i]) + ")"
                                line = paramHelicalFile.readline()
                                line_list = line.split()
				line_list.append(False)
                                while line_list[0] != element:
					line = paramHelicalFile.readline()
                                        line_list = line.split()
					line_list.append(False)
                                indFile.write(line.lstrip())
                                i += 1
			indFile.write("\n")

		elif line_list[0] == "(C)":
			#print "Made it to C"
			newNucleotides = []
			for element in nucleotides_to_analyze:
				if element == 1 and element != len(sequence):
					newNucleotides.append(1)
					newNucleotides.append(2)
				if element == len(sequence) and element != 1:
					if element - 1 not in newNucleotides:
						newNucleotides.append(element - 1)
					newNucleotides.append(element)
				else:
					if element - 1 not in newNucleotides:
						newNucleotides.append(element -1)
					if element not in newNucleotides:
						newNucleotides.append(element)
					if element + 1 not in newNucleotides:
						newNucleotides.append(element+1)
                        indFile.write(line.lstrip())
                        line = paramHelicalFile.readline()
                        indFile.write(line)
                        i = 0
                        while i < len(newNucleotides):
                                element = str(newNucleotides[i]) + ")"
                                line = paramHelicalFile.readline()
                                line_list = line.split()
				line_list.append(False)
                                while line_list[0] != element:
					line = paramHelicalFile.readline()
                                        line_list = line.split()
					line_list.append(False)
                                indFile.write(line.lstrip())
                                i += 1
                        line = paramHelicalFile.readline()
                        line_list = line.split()
			line_list.append(False)
                        indFile.write("\n")
                        while line_list[0] != "Inter_BP":
				line = paramHelicalFile.readline()
                                line_list = line.split()
				line_list.append(False)
                        indFile.write(line)
                        i = 0
                        while i < len(newNucleotides):
                                element = str(newNucleotides[i]) + ")"
                                line = paramHelicalFile.readline()
                                line_list = line.split()
				line_list.append(False)
                                while line_list[0] != element:
					line = paramHelicalFile.readline()
                                        line_list = line.split()
					line_list.append(False)
                                indFile.write(line.lstrip())
                                i += 1
                        indFile.write("\n")

		elif line_list[0] == "(D)":
			#print "Made it"
			indFile.write(line.lstrip())
			#For strand 1
			line = paramHelicalFile.readline()
			indFile.write(line)
			i = 0
			while i < len(nucleotides_to_analyze):
                                element = str(nucleotides_to_analyze[i]) + ")"
                                line = paramHelicalFile.readline()
                                line_list = line.split()
				line_list.append(False)
                                while line_list[0] != element:
					line = paramHelicalFile.readline()
                                        line_list = line.split()
					line_list.append(False)
                                indFile.write(line.lstrip())
                                i += 1
			indFile.write("\n")
			#For strand 2
			while line_list[0] != "Strand2":
				line = paramHelicalFile.readline()
				line_list = line.split()
				line_list.append(False)
			indFile.write(line.lstrip())
			i = 0
                        while i < len(nucleotides_to_analyze):
                                element = str(nucleotides_to_analyze[i]) + ")"
                                line = paramHelicalFile.readline()
                                line_list = line.split()
				line_list.append(False)
                                while line_list[0] != element:
					line = paramHelicalFile.readline()
                                        line_list = line.split()
					line_list.append(False)
                                indFile.write(line.lstrip())
                                i += 1
			indFile.write("\n")
			indFile.write("\n")
			paramHelicalFile.close()
			break
		else:
			line = paramHelicalFile.readline()
			
			
	#Working on Helical Frequency histogram file
	helical_param = open(helic_freqFileName,"r")
	indFile.write("\t\tFrequency histogram of helical parameters for " + indSequence + "\n")
	line = helical_param.readline()
	line_list = line.split()
	line_list.append(False)
	basesCovered = 1
	i = 0
	while i < len(nucleotides_to_analyze):
		if line_list[0] == "Nucleotide:" and basesCovered == nucleotides_to_analyze[i]:
			indFile.write("Base number: " + str(nucleotides_to_analyze[i]) + "\n")
			indFile.write(line)
			line = helical_param.readline()
			line_list = line.split()
			line_list.append(False)
			while line_list[0] != "Nucleotide:":
				indFile.write(line)
				line = helical_param.readline()
				line_list = line.split()
				line_list.append(False)
			basesCovered += 1
			i += 1
			#print i, nucleotides_to_analyze[i], basesCovered
			indFile.write("\n")
		elif line_list[0] == "Nucleotide:" and basesCovered != nucleotides_to_analyze[i]:
			basesCovered += 1
			line = helical_param.readline()
			line_list = line.split()
			line_list.append(False)
			while line_list[0] != "Nucleotide:":
				line = helical_param.readline()
				line_list = line.split()
				line_list.append(False)
		else:
			line = helical_param.readline()
			line_list = line.split()
			line_list.append(False) 
		#print nucleotides_to_analyze, basesCovered, line_list[0], basesCovered == nucleotides_to_analyze[i]
		#print i			

print("Finally all done!!!")


print("                       , ; ,   .-'''''-.   , ; ,")
print("                       \\|/  .'         '.  \|//")
print("                        \-;-/   ()   ()   \-;-/ ")
print("                        // ;               ; \\ ")
print("                       //__; :.         .; ;__\\")
print("                      `-----\'.'-.....-'.'/-----'")
print("                             '.'.-.-,_.'.'")
print("                               '(  (..-'")
print("                                 '-'")












