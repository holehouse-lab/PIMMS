## Function to create an n-dimensional array
 
# Main function to call
# typeOfitem 
#     Should be a class which can generate objects
#     e.g. float, int, complex, or any other type, such as 
#     myCoolClass
#
# dimensions
#     value for a 1D array, or a list or tuple defining the
#     dimensions, for higher order arrays. e.g. a 3D array 
#     might be [2,3,4]
#
def nDarray(typeOfitem, dimensions):
  depth = 0
  if type(dimensions) == int:
            dimensions = [dimensions]
  return(_recursiveAllocator(typeOfitem,dimensions,depth))

 
# recursive internal function
def _recursiveAllocator(basetype, dimensionList, depth):
  
  # Base case
  if depth == len(dimensionList)-1:
    currentDimension = dimensionList[depth]
    array = []
    for i in range(0,currentDimension):
      array.append(basetype())
    return array
 
  # Recursive case 
  else:
    array=[]
    currentDimension = dimensionList[depth]
 
    # for each element in each dimension recursivly
    # call the function
    for i in range(0,currentDimension):
      array.append(recursiveAllocator(basetype,dimensionList, depth+1))
    return array


