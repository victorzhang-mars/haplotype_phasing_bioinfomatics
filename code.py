import time
import numpy as np

def unmask_with_2(line):
    list_line = list(line)

    for index in range(0, len(list_line)):
        if list_line[index] == '*':
            list_line[index] = '2'

    line = ''.join(list_line)

    return line

def load_genome(file_path):
    fileInput = open(file_path, 'r')
    firstLine = fileInput.readline()

    firstLine = unmask_with_2(firstLine)
    firstLine = map(int, firstLine.strip('\n').split(' '))

    genotypes = [ [] for individual in range(0, len(firstLine)) ]
    for individual in range(0, len(firstLine)):
        genotypes[individual].append(firstLine[individual])
    for line in fileInput:
        line = unmask_with_2(line)
        line = map(int, line.strip('\n').split(' '))
        for individual in range(0, len(line)):
            genotypes[individual].append(line[individual])

    fileInput.close()
    return genotypes


def phaseGenotypeRecurse(genotype, i):
    if i == 0:
        if genotype[i] == 2:
            return [[1]]
        elif genotype[i] == 0:
            return [[0]]
        else:
            return [[1], [0]]
    elif i < 0:
        raise ValueError("Reached a negative index value in phaseGenotypeRecurse()")

    precedingPhases = phaseGenotypeRecurse(genotype, i-1)


    currentSNP = genotype[i]
    phases = []
    if (currentSNP == 2):
        for haplotype in range(0, len(precedingPhases)):
            newPhase = precedingPhases[haplotype]
            newPhase.append(1)
            phases.append(newPhase)

    elif (currentSNP == 0):
        for haplotype in range(0, len(precedingPhases)):
            newPhase = precedingPhases[haplotype]
            newPhase.append(0)
            phases.append(newPhase)

    elif (currentSNP == 1):
        for haplotype in range(0, len(precedingPhases)):
            phase1 = list(precedingPhases[haplotype])
            phase1.append(1)
            phase2 = precedingPhases[haplotype]
            phase2.append(0)
            phases.append(phase1)
            phases.append(phase2)
    else:
        raise ValueError("Genotype at index " + str(i)
                         + " is not a valid symbol from the set {0, 1, 2}! "
                           "Check that the provided genotype is correct. "
                           "Error in phaseGenotypeRecurse().")
    return phases


def phaseGenotype(genotype):
    lastIndex = len(genotype) - 1
    phases = phaseGenotypeRecurse(genotype, lastIndex)
    return phases


def haplotypeComplement(haplotype, genotype):

    if not isinstance(haplotype, list):
        raise ValueError("Haplotype provided to haplotypeComplement() is not a list!")

    if len(haplotype) != len(genotype):
        raise ValueError("The lengths of Haplotype and Genotype arguments "
                         "provided to haplotypeComplement() are not equal!")

    complement = []
    for i in range(0, len(haplotype)):
        genotypeVal = genotype[i]
        haplotypeVal = haplotype[i]
        if ((genotypeVal == 2 and haplotypeVal != 1) or (genotypeVal == 0 and haplotypeVal != 0)):
            if not isinstance(haplotype[i], int):
                raise ValueError("Haplotype at SNP index " + str(i)
                                 + " is not of type int! "
                                   "Check that the provided haplotype is a single haplotype "
                                   "and not a list of all possible genotype phases. Error in haplotypeComplement().")
            else:
                raise ValueError("Haplotype is not compatible with genotype at index " + str(i)
                                 + " in function haplotypeComplement()!")
        else:
            complement.append(genotypeVal - haplotypeVal)
    return complement


def calcHaplotypeProbability(expectation, numGenotypes):
    return expectation / (2 * numGenotypes)



def calcHaplotypePairingExpectation(haplotypePairProbability, totalProbability):
    return haplotypePairProbability / totalProbability



def calcHaplotypePairingProbability(haplotype_1, haplotype_2, haplotypeProbabilityDict):
    haplotype_1_str = ''.join(map(str, haplotype_1))
    haplotype_2_str = ''.join(map(str, haplotype_2))
    haplo_1_prob = haplotypeProbabilityDict[haplotype_1_str]
    haplo_2_prob = haplotypeProbabilityDict[haplotype_2_str]
    return haplo_1_prob * haplo_2_prob




def expectation_maximization(genotypes, numberOfIterations):

    if isinstance(genotypes, list):
        if len(genotypes) == 0:
            return []
        if isinstance(genotypes[0], list):
            if not isinstance(genotypes[0][0], int):
                raise ValueError("Genotype list provided to expectation_maximization() is not of the correct type! "
                                 "Type should be list of genotypes, which is a list of lists (of SNP=int).")
        else:
            raise ValueError("Genotype list provided to expectation_maximization() is not of the correct type! "
                             "Type should be a list of many genotypes, which is a list of lists (of SNP=int).")
    else:
        raise ValueError("Genotype list provided to expectation_maximization() is not of the correct type! "
                         "Type should be list of genotypes, which is a list of lists (of SNP=int).")


    genotypePhases = []
    genotypePhases_complement = []
    for geno in genotypes:
        possibleHaplotypes = []
        possibleHaplotypes_complement = []
        discoveredHaploDict = dict()
        phases = phaseGenotype(geno)
        for phase in phases:
            complement = haplotypeComplement(phase, geno)
            phase_str = ''.join(map(str, phase))
            complement_str = ''.join(map(str, complement))
            if phase_str not in discoveredHaploDict:
                discoveredHaploDict[phase_str] = 1
                discoveredHaploDict[complement_str] = 1
                possibleHaplotypes.append(phase)
                possibleHaplotypes_complement.append(complement)
        genotypePhases.append(possibleHaplotypes)
        genotypePhases_complement.append(possibleHaplotypes_complement)


    haplotypeProbDict = dict()
    numHaplotypes = 0
    for geno in range(0, len(genotypePhases)):
        for haplo in range(0, len(genotypePhases[geno])):
            haplotype = genotypePhases[geno][haplo]
            haplotype_str = ''.join(map(str, haplotype))
            if haplotype_str not in haplotypeProbDict:
                numHaplotypes += 1
                haplotypeProbDict[haplotype_str] = 0.0
            complement = genotypePhases_complement[geno][haplo]
            complement_str = ''.join(map(str, complement))
            if complement_str not in haplotypeProbDict:
                numHaplotypes += 1
                haplotypeProbDict[complement_str] = 0.0


    start = time.time()
    initial_haplo_prob = (1.0 / numHaplotypes)
    for haplo in haplotypeProbDict.iterkeys():
        haplotypeProbDict[haplo] = initial_haplo_prob


    haplotypePairProbs = []
    for geno in genotypePhases:
        numPhases = len(geno)
        pairProbabilities = [ (1.0 / numPhases) for t in range(0, numPhases) ]
        haplotypePairProbs.append(pairProbabilities)
    start = time.time()
    for iteration in range(1, numberOfIterations+1):
        for haplotypeOfInterest in haplotypeProbDict.keys():
            haplo_list = list(map(int, haplotypeOfInterest))
            sumOfPairProbs = 0.0
            for geno in range(0, len(genotypePhases)):
                for haplo in range(0, len(genotypePhases[geno])):
                    haplotype = genotypePhases[geno][haplo]
                    complement = genotypePhases_complement[geno][haplo]
                    if (haplotype == haplo_list) or (complement == haplo_list):
                        sumOfPairProbs += haplotypePairProbs[geno][haplo]
            haplotypeProbDict[haplotypeOfInterest] = calcHaplotypeProbability(sumOfPairProbs, len(genotypes))


        for geno in range(0, len(genotypePhases)):
            pairingProbs = [ 0.0 for t in range(0, len(genotypePhases[geno]))]
            totalProb = 0.0
            for haplo in range(0, len(genotypePhases[geno])):
                haplotype = genotypePhases[geno][haplo]
                genotype = genotypes[geno]
                complement = haplotypeComplement(haplotype, genotype)
                pairingProbs[haplo] = calcHaplotypePairingProbability(haplotype, complement, haplotypeProbDict)
                totalProb += pairingProbs[haplo]

            for haplo in range(0, len(genotypePhases[geno])):
                haplotypePairProbs[geno][haplo] = calcHaplotypePairingExpectation(pairingProbs[haplo], totalProb)


    maxProbabilityPhasesList = []
    for genotype in range(0, len(genotypePhases)):
        maxProbabilityPhase = []
        maxProbability = 0.0
        for phase in range(0, len(genotypePhases[genotype])):
            phase_list = genotypePhases[genotype][phase]
            probability = haplotypePairProbs[genotype][phase]
            if probability > maxProbability:
                maxProbability = probability
                maxProbabilityPhase = phase_list
        maxProbabilityPhaseComplement = haplotypeComplement(maxProbabilityPhase, genotypes[genotype])
        maxProbabilityPhasesList.append(maxProbabilityPhase)
        maxProbabilityPhasesList.append(maxProbabilityPhaseComplement)

    np_result = np.array(maxProbabilityPhasesList).transpose()

    end = time.time()
    print("EM took", end - start, "seconds.")

    return np_result





def em_windowed(genotypes, windowSize, numberOfIterations, fileOutputName):

    start = time.time()

    startPos = 0
    finishedWithAllWindows = False

    chunkAnswers = []
    while not finishedWithAllWindows:

        endOfWindow = min(startPos + windowSize, len(genotypes[0]))
        print(endOfWindow)

        if endOfWindow == len(genotypes[0]):
            finishedWithAllWindows = True

        genotypesToPhase = [ geno[startPos:endOfWindow] for geno in genotypes ]

        answer = expectation_maximization(genotypesToPhase, numberOfIterations)

        chunkAnswers.append(answer)

        startPos = endOfWindow

    end = time.time()
    print("This algorithm took", end - start, "seconds to run.")


    fileOutput = open(fileOutputName, 'w')
    for chunk in chunkAnswers:
        for haplo in chunk:
            stringHaplo = ' '.join(map(str, haplo))
            fileOutput.write(stringHaplo + '\n')
    fileOutput.flush()
    fileOutput.close()




def testDatasetFromClass():
    genotypesToPhase = [[1, 1, 1, 1, 2], [1, 0, 2, 2, 1], [2, 2, 0, 2, 1]]
    numIterations = 16
    result = expectation_maximization(genotypesToPhase, numIterations)
    print(result)




if __name__ == "__main__":

    #testDatasetFromClass()

    genotypes = load_genome("./test_data_masked.txt")
    em_windowed(genotypes, 8, 16, "test_data_sol.txt")

    #genotypes = load_genome("./example_data_1_masked.txt")
    #em_windowed(genotypes, 8, 16, "example_data_1_predicted_sol.txt")
