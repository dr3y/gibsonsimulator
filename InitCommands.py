from copy import deepcopy
guttDNA = "C:\\Users\\Andrey\\SkyDrive\\Guttlab\\DNA\\"
guttPython = "C:\\Users\\Andrey\\SkyDrive\\Guttlab\\Python\\"
guttJar = "C:\\Users\\Andrey\\SkyDrive\\Guttlab\\jarfiles"

def dc(inpt):
    """wrapper for deep copy"""
    return deepcopy(inpt)
def rc(seq):
    """reverse complement of the sequence!!"""
    revdict = {"A":"T",
               "T":"A",
               "G":"C",
               "C":"G",
               "a":"t",
               "t":"a",
               "g":"c",
               "c":"g"}
    return ''.join([revdict[a] for a in seq][::-1])
