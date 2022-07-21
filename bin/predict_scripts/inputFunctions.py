import numpy as np

from predict_scripts.components2tensor import components2tensor
#from components2tensor import components2tensor

def getTensorFromEuler (user_input):
    inputs = user_input.split(",")
    inputs = list(map(float, inputs))
    [Aax, Arh, a, b, g] = inputs
    susceptibility_tensor = components2tensor(Aax, Arh, a, b, g)
    return susceptibility_tensor

def getResidueNumbers (inputType):
    # This function handles manual input of residue numbers
    # There are three input types
    # inputType = 1 will accept an input like "1:4,7:9" for "1,2,3,4,7,8,9" (ranges)
    # inputType = 2 will accept two inputs - one like "1:9" and one like "5,6" for "1,2,3,4,7,8,9" (range and exclusion list)
    # inputType = 3 will accept a list of inputs like "1,2,3,4,7,8,9" (list)
    # The function will process the input and return a list of the input values like "[1, 2, 3, 4, 7, 8, 9]"
    # Note that spacing doesn't matter for the above comma separated lists, so "1, 2, 3" is the same as "1,2,3"

    output = []

    if (inputType == 1):
        userInput = input().split(",")
        for inputRange in userInput:
            inputRangeList = inputRange.split(":")
            beginning = int(inputRangeList[0])
            end = int(inputRangeList[1]) + 1 # since Python ranges are exclusive
            output.append(range(beginning, end))

    elif (inputType == 2):
        rangeInput = input().split(":")
        exclusionList = map(int, input().split(","))
        beginning = int(rangeInput[0])
        end = int(rangeInput[1]) + 1 # since Python ranges are exclusive
        output.extend(num for num in range(beginning, end) if num not in exclusionList)
    
    elif (inputType == 3):
        userInput = input().split(",")
        output = map(int, userInput)

    else:
        raise Exception("Invalid input type selected.")

    return (output)


def getSusceptibilityTensor (inputType = 2):
    # This function handles manual input of the susceptibility tensor
    # There are two input types
    # inputType = 1 will accept five numerical inputs
    # inputType = 2 will accept the axial and rhombic components and three euler angles

    if (inputType == 1):
        susceptibility_tensor = np.array([])
        for _ in range (5):
            susceptibility_tensor = np.append(susceptibility_tensor, float(input()))

    if (inputType == 2):
        inputs = []
        for _ in range (5):
            inputs.append(input())
        inputs = list(map(float, inputs))
        [Aax, Arh, a, b, g] = inputs
        susceptibility_tensor = components2tensor(Aax, Arh, a, b, g)
    
    return susceptibility_tensor
