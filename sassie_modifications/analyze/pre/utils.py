import numpy as np
import itertools

def hard_replace (arr, a, b): # replaces all occurences of {a} with {b} in an array {arr} of any dimension
    flat = arr.flatten()
    new_arr = np.where(flat == a, b, flat).reshape(arr.shape)
    return new_arr

def digits2nums (digitlist): # turns a list of digits into an array of numbers (e.g. ['0' '5' '2' '.' '9' '0' '-' '9' '0' '.' '0' '1'] => [52.9, -90.01])
    digitlist = list(digitlist)
    numlist = [list(y) for x, y in itertools.groupby(digitlist, lambda z: z == '0') if not x]
    print (numlist)
    numlist = list(map(lambda x: float("".join(x)), numlist))
    return np.asarray(numlist)


def findstr (str1, str2):
    if not isinstance(str1, str): str1 = "".join(str1.tolist())
    if not isinstance(str2, str): str2 = "".join(str2.tolist())
    if (len(str2) > len(str1)):
        return findstr(str2, str1) # ensures that str1 is the bigger element
    indices = []
    for i in range(len(str1) - len(str2) + 1):
        strcheck = "".join(str1[i: i + len(str2)])
        if strcheck == str2:
            indices.append(i)

    return np.array(indices)


def readasci(filename):
    with open(filename) as my_file:
        lines = my_file.readlines()

    length = max(list(map(lambda x: len(list(x)), lines)))
    lines = list(map(lambda x: x.ljust(length), lines))
    lines = list(map(lambda x: list(x), lines))
    text = np.asarray(lines, dtype = "str")
    print (text.shape)
    return text

def str2num (array, integer = False): # python equivalent of the str2num function in Matlab
    if integer == True:
        return (list(map(lambda x: int(eval("".join(x).strip())), array)))
    else:
        return (list(map(lambda x: float(eval("".join(x).strip())), array)))
