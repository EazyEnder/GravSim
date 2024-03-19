import numpy as np

def deserializeArray(string):
    """Transform a string to an array
        >Args: string: The string like x00,x01;x10;x11;...#\n
        >Return an array of tuple [[x00,x01],[x10,x11],...]"""
    array = string.split("#")[0]
    result = []
    for c in array.split(";"):
        result.append([float(c.split(",")[0]),float(c.split(",")[1])])
    return np.array(result)

def serializePositions(positions):
     """Serialize a list of positions to string
        >Args: positions: the list [[x00,x01],[x10,x11],...]
        >Return a string like : x00,x01;x10,x11;...#"""
     string = ""
     for i in range(len(positions)):
          pos = positions[i]
          if(i != 0):
            string += ";"
          string += str(round(pos[0],6))
          string += ","
          string += str(round(pos[1],6))
     return string + "#"