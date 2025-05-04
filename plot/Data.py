from numpy import loadtxt,fromfile,float64,int32

def read_parameters(prmts):
    """------------------------------------------------------------------------
    Read system's parameters from input file.
    ---------------------------------------------------------------------------
    Prmts:
    0       # Nb: Number of bodies in the galaxy
    1       # M : Mass of Galaxy
    2       # r : Radius of Galaxy
    3       # dt: Time step
    4       # steps: Evolution steps
    5       # jump: Data storage interval
    ------------------------------------------------------------------------"""
    data = loadtxt("../input", max_rows=6)
    return data[prmts]

def read_full_data(file_path, N):
    with open(file_path, "rb") as f:
        Pos = fromfile(f, dtype=float64, count=3 * N).reshape(N, 3)
        Vec = fromfile(f, dtype=float64, count=3 * N).reshape(N, 3)
        Mass = fromfile(f, dtype=float64, count=N)       
        Ids = fromfile(f, dtype=int32, count=N)       
    return Pos, Vec, Mass, Ids
