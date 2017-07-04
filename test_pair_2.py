from pairshuffle import PairShuffle
from Crypto.PublicKey import ElGamal
from Crypto.Random import random
from Crypto import Random
import binascii
import random


def main():
    Resoolt3 = 0
    Resoolt1 = [None] * 5
    Resoolt2 = [None] * 5
    k = 5  # Number in list
    ##key = ElGamal.generate(512, Random.new().read)
##    q = (key.p - 1) / 2
    p = 10198267722357351868598076141027380280417188309231803909918464305012113541414604537422741096561285049775792035177041672305646773132014126091142862443826263
    q = 5099133861178675934299038070513690140208594154615901954959232152506056770707302268711370548280642524887896017588520836152823386566007063045571431221913131
    g = 4
    h = random.getrandbits(512)
    H = pow(g, h, p)
    # Make client keys
    c = [None] * k
    C = [None] * k
    X = [None] * k
    Y = [None] * k
    for i in range(0, k):
        c[i] = random.getrandbits(512)
        C[i] = pow(g, c[i], p)
    ############
    for i in range(0, k):
        r = random.getrandbits(512)
        X[i] = pow(g, r, p)
        Y[i] = pow(H, r, p)
        Y[i] = (Y[i] * C[i]) % p
    pairshuffle_obj = PairShuffle(p, 5)
    Xbar, Ybar = pairshuffle_obj.go_shuffle_shuffle(p, q, g, H, X, Y)
    b_var = pairshuffle_obj.go_shuffle_verify(p, q, g, H, X, Y, Xbar, Ybar)
    print b_var


if __name__ == "__main__":
    main()
