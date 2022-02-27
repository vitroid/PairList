def load(file):
    """
    read a gro file
    """
    line = file.readline()
    natom = int(file.readline())
    atoms = dict()
    for i in range(natom):
        line = file.readline()
        atomname = line[10:15].replace(' ', '')
        if atomname not in atoms:
            atoms[atomname] = []
        x, y, z = [float(x) for x in line[20:].split()]
        atoms[atomname].append((x, y, z))
    cell = [float(x) for x in file.readline().split()]
    return atoms, cell
