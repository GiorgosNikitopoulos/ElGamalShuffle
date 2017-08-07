import random
import gmpy2
import shuffle_random
from datetime import datetime
from random import randint
from simpleshuffle import SimpleShuffle

class PairShuffle(object):
    def __init__(self, modulus, k):
        self.modulus = modulus
        self.k = k
        self.p1Gamma = 0
        self.p1A = []
        self.p1C = []
        self.p1U = []
        self.p1W = []
        self.v2Zrho = []
        self.p3D = []
        self.p1Lamda1 = 0
        self.p1Lamda2 = 0
        self.v4Zlamda = -1
        self.p5Zsigma = []
        self.p5Ztau = -1
        self.pv6 = SimpleShuffle(modulus, k)

    # Alpha and beta are respectively X and Y
    def go_shuffle_prove(self, pi, modulus, order, generator, public, alpha, beta, neff_beta):
        """PairShuffle Prove function. Returns 1 if successful"""

        if len(alpha) != len(pi) or len(alpha) != len(beta):
            raise Exception("Mismatched Vector Length")
        k = self.k
        piinv = [None] * k
        for i in range(k):
            piinv[pi[i]] = i
        # Prover STEP 1
        u = []
        w = []
        a = []

        # u w a random lists
        for i in range(k):
            u.append(shuffle_random.shuffle_rand_int(shuffle_random.RANDOM_FUNC_CHOICE, 0, order - 1))
        for i in range(k):
            w.append(shuffle_random.shuffle_rand_int(shuffle_random.RANDOM_FUNC_CHOICE, 0, order - 1))
        for i in range(k):
            a.append(shuffle_random.shuffle_rand_int(shuffle_random.RANDOM_FUNC_CHOICE, 0, order - 1))
        tau0 = shuffle_random.shuffle_rand_int(shuffle_random.RANDOM_FUNC_CHOICE, 0, order - 1)
        for i in range(k):
            self.v2Zrho.append(shuffle_random.shuffle_rand_int(shuffle_random.RANDOM_FUNC_CHOICE, 0, order - 1))

        nu = shuffle_random.shuffle_rand_int(shuffle_random.RANDOM_FUNC_CHOICE, 1, order - 1)
        gamma = shuffle_random.shuffle_rand_int(shuffle_random.RANDOM_FUNC_CHOICE, 1, order - 1)

        # compute public commits
        self.p1Gamma = pow(generator, gamma, modulus)  # (21)
        wbetasum = tau0 % order
        self.p1Lamda1 = 1
        self.p1Lamda2 = 1
        for i in range(k):
            self.p1A.append(pow(generator, a[i], modulus))  # (21)
            temporary_variable = (gamma * a[pi[i]]) % order
            self.p1C.append(pow(generator, temporary_variable, modulus))  # (21)
            self.p1U.append(pow(generator, u[i], modulus))  # (21)
            temporary_variable = (gamma * w[i]) % order
            self.p1W.append(pow(generator, temporary_variable, modulus))  # (21)
            temporary_variable = (w[i] * neff_beta[pi[i]]) % order
            wbetasum = (wbetasum + temporary_variable) % order
            temporary_variable = (w[piinv[i]] - u[i]) % order
            temporary_variable_2 = pow(alpha[i], temporary_variable, modulus)
            self.p1Lamda1 = (
                self.p1Lamda1 * temporary_variable_2) % modulus  # (22)
            temporary_variable_2 = pow(beta[i], temporary_variable, modulus)
            self.p1Lamda2 = (
                self.p1Lamda2 * temporary_variable_2) % modulus  # (23)
        g_to_the_wbetasum = pow(generator, wbetasum, modulus)
        h_to_the_wbetasum = pow(public, wbetasum, modulus)
        self.p1Lamda1 = (g_to_the_wbetasum * self.p1Lamda1) % modulus  # (22)
        self.p1Lamda2 = (h_to_the_wbetasum * self.p1Lamda2) % modulus  # (23)

        # Verifier STEP 2
        if None in self.v2Zrho:
            raise Exception("Error, None in v2Zrho")

        B = [(pow(generator, self.v2Zrho[i], modulus) * gmpy2.invert(self.p1U[i], modulus)) for i in xrange(k)]

        # Prover step 3
        b = []
        for i in range(k):
            b.append((self.v2Zrho[i] - u[i]) % order)  # (25)

        d = [] # Init
        for i in range(k):
            d.append((gamma * b[pi[i]]) % order)  # (26)
            self.p3D.append(pow(generator, d[i], modulus))  # (26)
        if None in self.p3D:
            raise Exception("Error, p3D has a None")

        # Verifier step 4
        # Generate random Lamda for fourth step
        self.v4Zlamda = shuffle_random.shuffle_rand_int(shuffle_random.RANDOM_FUNC_CHOICE, 0, order - 1)

        if self.v4Zlamda == -1:
            raise Exception("Error, v4Zlamda is an inappropriate value")

        r = []
        for i in range(k):
            # lambda * betai in (27)
            temporary_variable = (self.v4Zlamda * b[i]) % order
            r.append((temporary_variable + a[i]) % order)  # (27)

        s = []
        for i in range(k):
            s.append((gamma * r[pi[i]]) % order)  # Compute sigma under (27)

        # Prover step 5
        self.p5Ztau = (-tau0) % order  # Compute -tau0 in (29)
        for i in range(k):
            self.p5Zsigma.append((w[i] + b[pi[i]]) % order)  # (28)
            self.p5Ztau = (
                self.p5Ztau + ((b[i] * neff_beta[i]) % order)) % order  # (29)

        if self.v4Zlamda == -1:
            raise Exception("Error, v4Zlamda is an inappropriate value")
        # Make the dictionary for p5
        # (30)
        return self.pv6.Prove(modulus, order, generator, gamma, r, s)

    def go_shuffle_verify(self, modulus, order, generator, public, alpha, beta, alphabar, betabar):
        """PairShuffle Verify function. Returns 1 if Verify is successful"""

        k = self.k  # length
        if len(alpha) != k or len(beta) != k or len(alphabar) != k or len(betabar) != k:
            raise Exception('Error, inappropriate vector lengths')
        # Check for error there if p1 is null
        if None in self.p1A \
            or None in self.p1U \
            or None in self.p1W \
            or None in self.p1C \
            or self.p1Gamma == 0 \
            or self.p1Gamma == None \
            or self.p1Lamda1 == 0:
            raise Exception("Error, in the p1 step")
        # Check for error there if v2 is null
        if None in self.v2Zrho:
            raise Exception("Error, None value in v2Zrho")

        # Verifier step 2
        B = [None] * k
        for i in range(k):
            P = pow(generator, self.v2Zrho[i], modulus)  # g = generator
            B = (P * gmpy2.invert(self.p1U[i], modulus)) % modulus  # (24)
        # Check for error there if p3 is null
        if None in self.p3D:
            raise Exception("Error, None value in p3D")
        # Check for error there if v4 is null
        if self.v4Zlamda == -1:
            raise Exception("Error, v4Zlamda is an inappropriate value")
        # Prover and Verifier step 6: simple k-shuffle
        self.pv6.Verify(modulus, order, generator, self.p1Gamma)

        # Check for error there if p5 is null
        if self.p5Ztau == -1:
            raise Exception("Error, None value in p5Ztau")

        # Verifier Step 7
        Phi1 = 1
        Phi2 = 1
        P = 0
        Q = 0
        for i in range(k):
            # (31)
            Phi1 = (Phi1 * pow(alphabar[i],
                               self.p5Zsigma[i], modulus)) % modulus
            Phi1 = (
                Phi1 * gmpy2.invert(pow(alpha[i], self.v2Zrho[i], modulus), modulus)) % modulus
            # (32)
            Phi2 = (Phi2 * pow(betabar[i],
                               self.p5Zsigma[i], modulus)) % modulus
            Phi2 = (
                Phi2 * gmpy2.invert(pow(beta[i], self.v2Zrho[i], modulus), modulus)) % modulus
            if pow(self.p1Gamma, self.p5Zsigma[i], modulus) != (self.p1W[i] * self.p3D[i]) % modulus:
                raise Exception("invalid PairShuffleProof")
        if (self.p1Lamda1 * pow(generator, self.p5Ztau, modulus)) % modulus != Phi1:  # (33)
            raise Exception("invalid PairShuffleProof")

        if (self.p1Lamda2 * pow(public, self.p5Ztau, modulus)) % modulus != Phi2:  # (33)
            raise Exception("invalid PairShuffleProof")
        return 1

    def go_shuffle_shuffle(self, modulus, order, generator, public, alpha, beta):
        """Start the shuffle. (Neffs Elgamal PairShuffle)"""
        k = len(alpha)  # Length of X(alpha)
        if k != len(beta):
            raise Exception("alpha,beta vectors have inconsistent length")
        pi = range(k)

        for i in range(k - 1, 0, -1):  # Permutation array
            j = shuffle_random.shuffle_rand_int(shuffle_random.RANDOM_FUNC_CHOICE, 0, i)
            if j != i:
                temporary_variable = pi[j]
                pi[j] = pi[i]
                pi[i] = temporary_variable
        neff_beta = [] # Initializing BETA
        for i in range(k):
            neff_beta.append(shuffle_random.shuffle_rand_int(shuffle_random.RANDOM_FUNC_CHOICE, 0, order - 1))
        XBar = []  # Initializing XBar
        YBar = []  # Initializing YBar

        for i in range(k):
            XBar_tmp = pow(generator, neff_beta[pi[i]], modulus)  # (17)
            XBar.append((XBar_tmp * alpha[pi[i]]) % modulus)  # (17)
            YBar_tmp = pow(public, neff_beta[pi[i]], modulus)  # (17)
            YBar.append((YBar_tmp * beta[pi[i]]) % modulus)  # (17)
        self.go_shuffle_prove(pi, modulus, order,
                              generator, public, alpha, beta, neff_beta)
        return XBar, YBar
