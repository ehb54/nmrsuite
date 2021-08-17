import numpy as np

def findstr (str1, str2):
    if not isinstance(str1, str): str1 = "".join(str1.tolist())
    if not isinstance(str2, str): str2 = "".join(str2.tolist())
    if (len(str2) > len(str1)):
        return findstr(str2, str1) # Ensures that str1 is the bigger element
    indices = []
    for i in range(len(str1) - len(str2) + 1):
        strcheck = "".join(str1[i: i + len(str2)])
        if strcheck == str2:
            indices.append(i)

    return np.array(indices)

def hard_replace (arr, a, b): # Replaces all occurences of {a} with {b} in an array {arr} of any dimension
    flat = arr.flatten()
    new_arr = np.where(flat == a, b, flat).reshape(arr.shape)
    return new_arr

def str2num (array, integer = False): # Python equivalent of the str2num function in Matlab
    if integer == True:
        return (list(map(lambda x: int(eval("".join(x).strip())), array)))
    else:
        return (list(map(lambda x: float(eval("".join(x).strip())), array)))