import random
import gmpy2
import shuffle_random
from datetime import datetime
from random import randint

class SimpleShuffle(object):
    def __init__(self, modulus, k):
        self.p0Y = []
        self.p0X = []
        self.p2Theta = [None] * (2 * k)
        self.p4Zalpha = [None] * (2 * k - 1)
        self.v1Zt = 0
        self.v3Zc = 0
        self.modulus = modulus
        self.k = k

    def Prove(self, modulus, order, G, gamma, x, y):
        """Prove function in Section 3 of Neff's paper"""
        k = self.k
        if k <= 1:
            raise Exception("Can't Shuffle length 1 vector")
        if k != len(y):
            raise Exception("Mismatched vector lengths")

        # Step 0: xi = logG(Xi), Xi = G^xi, same for yi.
        # Basically creates Yi and Xi
        for i in range(k):
            self.p0X.append(pow(G, x[i], modulus))
            self.p0Y.append(pow(G, y[i], modulus))

        if None in self.p0X  or None in self.p0Y:
            raise Exception('Error, p0X or p0Y has a None')

        # Verifier Step 1: create t in Zq
        self.v1Zt = shuffle_random.shuffle_rand_int(shuffle_random.RANDOM_FUNC_CHOICE, 0, order - 1)
        t = self.v1Zt

        if self.v1Zt == None:
            raise Exception('v1Zt is None')

        # Prover step 2
        gamma_t = (gamma * t) % order  # (gamma * t) variable
        x_hat = []  # Inits
        y_hat = []

        for i in range(k):
            x_hat.append((x[i] - t) % order)  # (5)
            y_hat.append((y[i] - gamma_t) % order)  # (6)

        #(7) theta and Theta vectors: Start
        thlen = (2 * k) - 1
        theta = []
        Theta = [None] * (thlen + 1)
        for i in range((2 * k) - 1):
            theta.append(shuffle_random.shuffle_rand_int(shuffle_random.RANDOM_FUNC_CHOICE, 0, order - 1))
        Theta[0] = thenc(modulus, order, G, None, None, theta[0], y_hat[0])
        for i in range(1, k):
            Theta[i] = thenc(modulus, order, G, theta[i - 1],
                             x_hat[i], theta[i], y_hat[i])
        for i in range(k, thlen):
            Theta[i] = thenc(modulus, order, G, theta[i - 1],
                             gamma, theta[i], None)
        Theta[thlen] = thenc(
            modulus, order, G, theta[thlen - 1], gamma, None, None)
        self.p2Theta = Theta

        if None in self.p2Theta:
            raise Exception('Error, p2Theta has a None')
        #(7) theta and Theta vectors: End

        # Verifier Step 3
        self.v3Zc = shuffle_random.shuffle_rand_int(shuffle_random.RANDOM_FUNC_CHOICE, 0, order - 1)
        c = self.v3Zc

        if self.v3Zc == None:
            raise Exception('Error, v3Zc is None')

        # Prover step 4
        alpha = [None] * thlen
        runprod = c
        # (8)
        for i in range(k):
            # The one multiplication Neff was reffering to
            runprod = (runprod * x_hat[i]) % order
            # The one division Neff was referring to
            runprod = (runprod * gmpy2.invert(y_hat[i], order)) % order
            alpha[i] = (theta[i] + runprod) % order
        gammainverse = gmpy2.invert(gamma, order)
        rungamma = c
        # That is the second part of (8)
        for i in range(1, k):
            rungamma = (rungamma * gammainverse) % order
            alpha[thlen - i] = (theta[thlen - i] + rungamma) % order

        # Verifier step 5
        self.p4Zalpha = alpha

        if None in self.p4Zalpha:
            raise Exception('None in p4Zalpha')

        return 1

    def Verify(self, modulus, order, G, Gamma):
        """Verifier for Neff's SimpleShuffle"""
        X = self.p0X
        Y = self.p0Y
        Theta = self.p2Theta
        alpha = self.p4Zalpha
        # Validate vector lens
        k = len(Y)
        thlen = (2 * k) - 1
        if k <= 1 or len(Y) != k or len(Theta) != thlen + 1 or len(alpha) != thlen:
            raise Exception('Something went wrong with vector lengths')

        if None in self.p0X  or None in self.p0Y:
            raise Exception('None in p0X or p0Y')
        if self.v1Zt == None:
            raise Exception('v1Zt is None')
        if None in self.p2Theta:
            raise Exception('p2Theta has a None')
        if self.v3Zc == None:
            raise Exception('v3Zc is None')
        if None in self.p4Zalpha:
            raise Exception('p4Zalpha has a None')
        t = self.v1Zt
        c = self.v3Zc

        # Verifier step 5
        negt = (-t) % order  # find -t
        U = pow(G, negt, modulus)  # (10)
        W = pow(Gamma, negt, modulus)  # (10)
        X_hat = []  # Init
        Y_hat = []
        for i in range(k):
            X_hat.append((X[i] * U) % modulus)  # (11)
            Y_hat.append((Y[i] * W) % modulus)  # (11)
        P = 0
        Q = 0
        s = 0
        #(12) Start
        b_good = True
        b_good = b_good and thver(
            X_hat[0], Y_hat[0], Theta[0], P, Q, c, alpha[0], modulus, order)
        for i in range(1, k):
            b_good = b_good and thver(
                X_hat[i], Y_hat[i], Theta[i], P, Q, alpha[i - 1], alpha[i], modulus, order)
        for i in range(k, thlen):
            b_good = b_good and thver(
                Gamma, G, Theta[i], P, Q, alpha[i - 1], alpha[i], modulus, order)
        b_good = b_good and thver(
            Gamma, G, Theta[thlen], P, Q, alpha[thlen - 1], c, modulus, order)
        if not b_good:
            raise Exception('Incorrect poof')
        #(12) End
        return


def thenc(modulus, order, G, a, b, c, d):
    """Helper function in order to compute G^(ab-cd)"""
    ab = 0
    cd = 0
    if a == None or a == 0:
        ab = 0
    else:
        ab = (a * b) % order
    if c == None or c == 0:
        cd = 0
    else:
        if d == None or d == 0:
            cd = c
        else:
            cd = (c * d) % order
    return pow(G, (ab - cd) % order, modulus)


def thver(A, B, T, P, Q, a, b, modulus, order):
    """Helper function in order to verify Theta elements"""
    P = pow(A, a, modulus)
    Q = pow(B, ((-b) % order), modulus)
    P = (P * Q) % modulus
    if P == T:
        return True
    else:
        return False
