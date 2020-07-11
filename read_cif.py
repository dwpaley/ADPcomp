from math import cos, radians
def read_cif(cifName, gij, xyzDict, UijDict, sUijDict):
    '''takes a cifName string and the given objects (gij is a Numpy array, all
    others are dictionaries). Constructs the gij tensor and sets the values
    in the array. For all the others, puts the data in the dictionaries with
    atom names as keys and lists [x,y,z] or [U11, U22, U33, U23, U13, U12] as
    values.
    '''
    with open(cifName, 'r') as cifFile:
        cifLines = cifFile.readlines()





    


    def _parse_entry(line, start, end):
        '''Takes a single line from a cif file as a string. For the entries 
        [start:end] (regular slice notation), checks for an esd and then writes
        the value and the esd (if found) to a pair of lists. When all the 
        requested entries in the line have been read, the values and esd's lists
        are returned as a 2-tuple of lists.
        '''

        nums, sNums = [], []
        for entry in line.split()[start:end]:
            if '(' in entry:
                numStr = entry[:entry.find('(')]
                sNumStr = entry[ (1+entry.find('(')) : (entry.find(')')) ]
                signif = -(len(numStr) - numStr.find('.') - 1) if '.' in numStr\
                            else 0 #this is the exponent of the smallest sig fig
                nums.append(float(numStr))
                sNums.append(float(sNumStr) * 10**signif)
            else: 
                nums.append(float(entry))
                sNums.append(0.0)
        return (nums, sNums)

    #Find and read the cell parameters; construct the metric tensor and write
    #it to gij.
    gij.dtype = float
    line, cellPars = '', []
    cifIter = iter(cifLines)
    while not line.startswith('_cell_length_a'): line = next(cifIter)
    for i in range(6):
        cellPars.append(_parse_entry(line, 1, 2)[0][0])
        line = next(cifIter)
    a, b, c = cellPars[0:3]
    cal, cbe, cga = [cos(radians(angle)) for angle in cellPars[3:6]]
    gij[0] = [a**2, a*b*cga, a*c*cbe]
    gij[1] = [b*a*cga, b**2, b*c*cal]
    gij[2] = [c*a*cbe, c*b*cal, c**2]

    #Find the atoms sites and add each atom x,y,z to the xyzDict as a list.
    line = ''
    cifIter = iter(cifLines)
    while not line.strip().startswith('_atom_site_label'):
        line = next(cifIter)
    while line.strip().startswith('_'):
        line = next(cifIter)
    while line.split():
        xyzDict[line.lower().split()[0]] = _parse_entry(line, 2, 5)[0]
        line=next(cifIter)

    #find the atom aniso_U entries and add U11...U12 and sU11...sU12
    #to the UijDict and sUijDict as lists.
    line = ''
    cifIter = iter(cifLines)
    while not line.strip().startswith('_atom_site_aniso_label'):
        line = next(cifIter)
    while line.strip().startswith('_'):
        line = next(cifIter)
    while line.split():
        values = _parse_entry(line, 1, 7)
        UijDict[line.lower().split()[0]] = values[0]
        sUijDict[line.lower().split()[0]] = values[1]
        line = next(cifIter)

