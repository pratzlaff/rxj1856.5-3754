def qdp_files(year):
    trans = { 'neg' : 'm', 'pos' : 'p' }
    return { o : './fits/qe_N0013/rxj1856_bbody_{}_{}_orders_5.qdp'.format(year, trans[o]) for o in ('neg', 'pos') }

def qdp_ratios(qdpfile):
    wav, orders = read_qdp(qdpfile)
    ratio = orders[0] / orders.sum(axis=0)
    return wav[::-1], ratio[::-1]

def read_qdp(qdpfile):
    data = np.genfromtxt(qdpfile, skip_header=3, loose=True, unpack=True)
    n = int((data.shape[1]-1)/2)
    return data[0,:n], data[5:,:n]

def get_ratios(year):
    files = qdp_files(year)
    ratios = { }
    for order in files:
        wav, ratio = qdp_ratios(files[order])
        ratios[order] = { 'wav': wav, 'ratio' : ratio }
    return ratios
