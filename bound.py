from math import prod
from mpmath import *
mp.dps=15

def prwbc(ks, B):
    aA_s = [[],[(1,1)]]
    bBm_s = [(0,1)]
    for k in ks:
        bBm_s.append((1,k))
    bB_s = [bBm_s,[]]
    #print("prwbc", B)

    #return foxh(aA_s, bB_s, B * B, zeroprec=1000, maxterms=100000, maxprec=5000)
    return foxh(aA_s, bB_s, B * B, maxterms=1000000, maxprec=10000)

def test_prwbc_1():
    # Degenerate to Rayleigh distribution
    ks = [1]
    B = 10

    lhs = prwbc(ks, B)
    rhs = exp(-B * B)

    print(lhs, rhs)

    assert almosteq(lhs, rhs)

def findB(ks, eps, Bhint):
    # find a B such that prwbc(ks, B) < eps

    Bmin = Bhint
    Bmax = Bhint * 2
    while prwbc(ks, Bmax) > eps:
        Bmax *= 2
    while Bmax - Bmin > 0.1 and (Bmax - Bmin) / Bmax > 0.01:
        B = (Bmin + Bmax) / 2
        if prwbc(ks, B) < eps:
            Bmax = B
        else:
            Bmin = B

    return Bmax

def totalIndep(n, eps, Bhint):
    # for n terms
    # s^{n/2} * \prod^n e_i
    assert n % 2 == 0
    n2 = n // 2
    ks = [1] * n

    # find a B such that prwbc(ks, B) < eps
    B = findB(ks, eps, Bhint)
    return B, prwbc(ks, B)

def indep(n, eps, Bhint):
    # for 2n terms
    # s^n * \prod^n e_i
    assert n % 2 == 0
    n2 = n // 2

    ks = [n2] + [1] * n2

    # find a B such that prwbc(ks, B) < eps

    B = findB(ks, eps, Bhint)

    return B, prwbc(ks, B)

def plotIndep():
    # plot prwbc for indep
    # n for different lines
    # B as the axis
    vals = {}
    for n in range(1, 20, 10):
        ks = [n] + [1] * n
        # B from 0 to 1, step 0.05
        for B in [x * 0.05 for x in range(1, 21)]:
            val = prwbc(ks, B)
            vals[(n, B)] = val
        print(n)

    import matplotlib.pyplot as plt
    for n in range(1, 20, 10):
        Bs = []
        ps = []
        for B in [x * 0.05 for x in range(1, 21)]:
            Bs.append(B)
            ps.append(vals[(n, B)])
        # print ps in log
        print(ps)
        print([log2(p) for p in ps])
        plt.plot(Bs, ps, label=f'n={n}')
    plt.yscale('log')
    plt.xlabel('B')
    plt.ylabel('prwbc')
    plt.title('prwbc for indep')
    plt.legend()
    plt.show()

#plotIndep()

def dep(n, eps, Bhint):
    # for 2n terms
    # s^n * mu^{n-1} * e
    assert n % 2 == 0
    n2 = n // 2

    ks = [n2, n2-1, 1]
    if n2 == 1:
        ks = [1, 1]

    B = findB(ks, eps, Bhint)
    return B, prwbc(ks, B)

def bounds():
    import sys
    # get exp from argv[1]
    exp = int(sys.argv[1]) if len(sys.argv) > 1 else '28'
    print("working with exp = ", exp)

    # erfc(6/sqrt(2))
    #eps = mpf(1.97318e-9)
    eps = mpf(2 ** (-exp))

    B1s = []
    B2s = []
    B3s = []

    for n2 in range(1, 32):
        #Bindep, _ = indep(n2 * 2, eps, Bhint=B1s[-1] if B1s else 1)
        Bdep, _ = dep(n2 * 2, eps, Bhint=B2s[-1] if B2s else 1)
        #B3dep, _ = totalIndep(n2 * 2, eps)
        #B1s.append(Bindep)
        B2s.append(Bdep)
        #B3s.append(B3dep)
        #print(n2, Bindep)
        print(n2, Bdep)

    #assert B1s[0] == B2s[0]
    #base = B1s[0]

    # D, D^2, D^3
    #canonBs = [B1s[0]]
    #for i in range(len(B1s) - 1):
    #    canonBs.append(canonBs[-1] * B1s[0])

    #print(B1s)
    print(B2s)
    #print(B3s)
    #print(canonBs)

    # take log and plot
    #log2Bs = [log(b, base) for b in B1s]
    #print("B1s", log2Bs)
    #log2Bs2 = [log(b, base) for b in B2s]
    #print("B2s", log2Bs2)
    #log2Bs3 = [log(b, base) for b in B3s]
    #print("B3s", log2Bs3)
    #log2canonBs = [log(b, base) for b in canonBs]
    #print("canonBs", log2canonBs)

# bounds()

def plot():
    # plot a curve for log2 Bs
    import matplotlib.pyplot as plt

    #plt.plot(range(1, len(log2Bs) + 1), log2Bs, label='indep')
    #plt.plot(range(1, len(log2Bs2) + 1), log2Bs2, label='dep')
    #plt.plot(range(1, len(log2Bs3) + 1), log2Bs3, label='total indep')
    #plt.plot(range(1, len(log2canonBs) + 1), log2canonBs, linestyle='dashed', label='canon')
    # plot the original Bs
    #plt.plot(range(1, len(B1s) + 1), B1s, label='indep')
    #plt.plot(range(1, len(B2s) + 1), B2s, label='dep')
    plt.plot(range(1, len(B3s) + 1), B3s, label='total indep')
    #plt.plot(range(1, len(canonBs) + 1), canonBs, linestyle='dashed', label='canon')
    plt.xlabel('n2')
    plt.ylabel('B')
    plt.title('B vs n2')
    plt.legend()
    plt.show()

def lambert(k, eps):
    # k = (1, 1, 1, ....)
    T = (k - 1) / 2
    c = k

    C1 = ((2 * pi) ** T) * (k ** (-1/2))

    w = - (T / c) * lambertw(- (c / T) * ((eps/C1) ** (1 / T)), -1)
    B = w ** (k/2)
    return B
    # k = (k)
    #T = 0
    #c = 1
    #w = - (T / c) * lambertw(- (c / T) * (eps ** (1 / T)), -1)
    #B = w ** (k/2)
    #return B

def lambertIndep(k, eps):
    # k = (n, 1, 1, ....)
    assert k % 2 == 0
    n2 = k // 2
    u = n2 + 1

    T = (u - 1) / 2
    c = k / (n2 ** (n2 / k))

    C1 = ((2 * pi) ** T) * (k ** (-1/2)) * (n2 ** (1/2 + (1-u) * n2 / (2 * k)))

    w = - (T / c) * lambertw(- (c / T) * ((eps/C1) ** (1 / T)), -1)
    B = w ** (k/2)
    return B

def lambderGeneral(ks, eps):
    # ks = (k1, k2, ..., kr)
    u = len(ks)
    k = sum(ks)
    T = (u - 1) / 2

    c = k / prod([ki ** (ki / k) for ki in ks])
    C1 = ((2 * pi) ** T) * (k ** (-1/2)) * prod([ki ** (1/2 + (1-u) * ki / (2 * k)) for ki in ks])

    w = - (T / c) * lambertw(- (c / T) * ((eps/C1) ** (1 / T)), -1)
    B = w ** (k/2)
    return B

def test_lambert(e):
    eps = mpf(2 ** (-e))
    Bhint = 1
    ret = []
    for n in range(1, 33):
        ks = [n] + [1] * n
        #ks = [n] + ([n-1, 1] if n > 1 else [1])
        B_lambert = lambderGeneral(ks, eps)
        ret.append(float(B_lambert))
        #B_prwbc = findB(ks, eps, Bhint)
        #Bhint = B_prwbc # update
        #print(n, B_lambert, B_prwbc)

        #prwbc_val = prwbc(ks, B_lambert)
        #prwbc_find = prwbc(ks, B_prwbc)
        #print("  prwbc at lambert:", log2(prwbc_val))
        #print("  prwbc at findB:", log2(prwbc_find))
    # print ret in c++ format, i.e. std::map<int, double> lambert_bound_eXX = {
    print(f"std::map<int, double> log2B_e{e} = {{")
    for i, v in enumerate(ret):
        print(f"    {{{i+1}, {log2(v)}}},")
    print("};")
    return ret

#test_lambert(28)
#
res = {}
for e in range(16, 40, 4):
    res[e] = test_lambert(e)
print(res)

def compare():
    eps = mpf(2 ** -28)
    Bhint = 1
    #for n in range(1, 10):
    #    B3dep, _ = totalIndep(n * 2, eps, Bhint)
    #    Bhint = B3dep # update
    #    #can, _ = totalIndep(n * 2, eps)
    #    Blambert = lambert(2 * n, eps)
    #    print(n, B3dep, Blambert)
    for n in range(1, 10):
        B3dep, _ = indep(n * 2, eps, Bhint)
        Bhint = B3dep # update
        #can, _ = totalIndep(n * 2, eps)
        Blambert = lambertIndep(2 * n, eps)
        print(n, B3dep, Blambert)

#compare()
