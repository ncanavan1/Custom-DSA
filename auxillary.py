import numpy as np
import random
import os
from Crypto.Hash import SHAKE256
from Crypto.Hash import SHAKE128


q = 8380417  # A prime number
zeta = 1753
d = 13
tau = 39
lam = 128
gamma1 = 2**17
gamma2 = (q - 1)/88
k = 4
l = 4
eta = 2    # A small positive integer
beta = tau*eta
omega = 80

zetas = [0, 4808194, 3765607, 3761513, 5178923, 5496691, 5234739, 5178987,  7778734, 3542485, 
         2682288, 2129892, 3764867, 7375178, 557458, 7159240,  5010068, 4317364, 2663378, 6705802, 
         4855975, 7946292, 676590, 7044481,  5152541, 1714295, 2453983, 1460718, 7737789, 4795319, 
         2815639, 2283733,  3602218, 3182878, 2740543, 4793971, 5269599, 2101410, 3704823, 1159875,  
         394148, 928749, 1095468, 4874037, 2071829, 4361428, 3241972, 2156050,  3415069, 1759347, 
         7562881, 4805951, 3756790, 6444618, 6663429, 4430364,  5483103, 3192354, 556856, 3870317, 
         2917338, 1853806, 3345963, 1858416,  3073009, 1277625, 5744944, 3852015, 4183372, 5157610, 
         5258977, 8106357,  2508980, 2028118, 1937570, 4564692, 2811291, 5396636, 7270901, 4158088,  
         1528066, 482649, 1148858, 5418153, 7814814, 169688, 2462444, 5046034,  4213992, 4892034, 
         1987814, 5183169, 1736313, 235407, 5130263, 3258457,  5801164, 1787943, 5989328, 6125690, 
         3482206, 4197502, 7080401, 6018354,  7062739, 2461387, 3035980, 621164, 3901472, 7153756, 
         2925816, 3374250,  1356448, 5604662, 2683270, 5601629, 4912752, 2312838, 7727142, 7921254,  
         348812, 8052569, 1011223, 6026202, 4561790, 6458164, 6143691, 1744507,  1753, 6444997, 
         5720892, 6924527, 2660408, 6600190, 8321269, 2772600,  1182243, 87208, 636927, 4415111, 
         4423672, 6084020, 5095502, 4663471,  8352605, 822541, 1009365, 5926272, 6400920, 1596822, 
         4423473, 4620952,  6695264, 4969849, 2678278, 4611469, 4829411, 635956, 8129971, 5925040,  
         4234153, 6607829, 2192938, 6653329, 2387513, 4768667, 8111961, 5199961,  3747250, 2296099, 
         1239911, 4541938, 3195676, 2642980, 1254190, 8368000,  2998219, 141835, 8291116, 2513018, 
         7025525, 613238, 7070156, 6161950,  7921677, 6458423, 4040196, 4908348, 2039144, 6500539, 
         7561656, 6201452,  6757063, 2105286, 6006015, 6346610, 586241, 7200804, 527981, 5637006,  
         6903432, 1994046, 2491325, 6987258, 507927, 7192532, 7655613, 6545891,  5346675, 8041997, 
         2647994, 3009748, 5767564, 4148469, 749577, 4357667,  3980599, 2569011, 6764887, 1723229, 
         1665318, 2028038, 1163598, 5011144,  3994671, 8368538, 7009900, 3020393, 3363542, 214880, 
         545376, 7609976,  3105558, 7277073, 508145, 7826699, 860144, 3430436, 140244, 6866265,  
         6195333, 3123762, 2358373, 6187330, 5365997, 6663603, 2926054, 7987710,  8077412, 
         3531229, 4405932, 4606686, 1900052, 7598542, 1054478, 7648983]



def IntegerToBits(x,alpha):
    """
    Computes a base-2 representation of x mod 2**alpha using little endian order.

    Input: A non-negative integer x and a positive integer alpha
    
    Output: A bit string y of length alpha
    """

    y = np.zeros(alpha, dtype=int)
    x_dash = x
    for i in range(alpha):
        y[i] = x_dash%2
        x_dash = np.floor(x_dash/2)
    return y
    

def BitsToInteger(y,alpha):
    """
    Computes the integer value expressed by a bit string using little-endian order.

    Input: A positive integer alpha and a bit string y of length alpha

    Ouput: A non-negative integer x
    """

    x = 0
    for i in range(1,alpha+1): ## +1 makes work as need to go through bits
        x = 2*x +y[alpha-i]
    return x


def IntegerToBytes(x,alpha):
    """
    Computes a base-256 representation of x mod 256**alpha using little endian order.

    Input: A non negative integer x and a positive integer alpha

    Output: A byte string y of length alpha
    """
    y = np.zeros(alpha, dtype=int)
    x_dash = x
    for i in range(alpha -1):
        y[i] = x_dash%256
        x_dash = np.floor(x_dash/256)
    return y


def BitsToBytes(y):
    """
    Converts a bit string into a byte string using little endian order
    
    Input: A but string y of length alpha

    Output: A byte string z of length ceil(alpha/8)
    """
    alpha = len(y)
    z = np.zeros(int(np.ceil(alpha/8)), dtype=int)
    for i in range(alpha): ##alpha not alpha-1
        z[int(np.floor(i/8))] = z[int(np.floor(i/8))] + y[i]*2**(i%8)
    return z    


def BytesToBits(z):
    """
    Converts a byte string into a bit string using little endian order

    Input: A byte string z of length alpha

    Output: A bit string y of length alpha*8
    """
    alpha = len(z)
    y = np.zeros(8*alpha, dtype=int)
    z_dash = z
    for i in range(alpha): ##alpha not alpha-1
        for j in range(8):
            y[(8*i) + j] = z_dash[i]%2
            z_dash[i] = np.floor(z_dash[i]/2)
    return y



def CoeffFromThreeBytes(b0, b1, b2):
    """
    Generates an element of {0,1,2,...,q-1} U {⊥}

    Input: Three bytes b0, b1, b2

    Output: An integer mod q or ⊥
    """

    b2_dash = b2
    if b2_dash > 127:
        b2_dash = b2_dash - 128
    
    z = 2**16 * b2_dash + 2**8 * b1 + b0

    if z < q:
        return z
    else:
        return "⊥"
    
def CoeffFromHalfByte(b):
    """
    Let eta in {2,4}. Generates an element of {-eta, -eta+1, ..., eta} U {⊥}

    Input: Integer b in {0,1,...,15}

    Output: An integer between - eta and eta or ⊥
    """

    if eta == 2 and b < 15:
        return 2 - b%5
    else:
        if eta == 4 and b < 9:
            return 4 - b
        else:
            return "⊥" 

def SimpleBitPack(w,b):
    """
    Encodes a polynomial w into a byte string

    Input: b in N and w in Ring such that the coefficients of w are all in [0,b]

    Output: A byte string of length 32.bitlen_b
    """
    bitlen_b = b.bit_length()
    z = []
    for i in range(256):
        z.extend(IntegerToBits(w[i], bitlen_b))
    return BitsToBytes(z)

def BitPack(w, a, b):
    """
    Encodes a polynomial w into a byte string

    Input: a,b in N and w in Ring such that the coefficients of w are all in [-a,b]

    Output: A byte string of length 32.bitlen_(ab)
    """

    bitlen_ab = (a+b).bit_length()
    z = []
    for i in range(256):
        z.extend(IntegerToBits(b-w[i], bitlen_ab))
    return BitsToBytes(z)

def SimpleBitUnpack(v,b):
    """
    Reverses the procedure of SimpleBitPack

    Input: b in N and a byte string v of length 32.bitlen_b

    Output: A polynomial w in Ring with coeffcients in [0,2**c -1], where c = bitlen_b
    """

    bitlen_b = b.bit_length()
    c = bitlen_b
    z = BytesToBits(v)
    w = np.zeros(256, dtype=int)
    for i in range(256):
        w[i] = BitsToInteger(z[i*c:(i*c + c)],c) ##remove -1 from i*c+c
    return w


def BitUnpack(v,a,b):
    """
    Reverse the procedure of BitPack

    Input: a,b in N and a byte string v of length 32.bitlen_(ab)

    Output: A polynomial w in Ring with coefficients in [b-2**c + 1, b] where c = bitlen_ab.
    When a+b+1 is a power of 2, the coefficients are in [-a,b]
    """

    bitlen_ab = (a+b).bit_length()
    c = bitlen_ab
    z = BytesToBits(v)
    w = np.zeros(256, dtype=int)
    for i in range(256):
        w[i] = b - BitsToInteger(z[i*c:(i*c + c)],c)
    return w
 

def HintBitPack(h):
    """
    Encodes a polynomial vector 𝐡 with binary coefficients into a byte string.

    Input: A polynomial vector 𝐡 ∈ 𝑅_2^k such that the polynomials 𝐡[0], 𝐡[1],...,𝐡[𝑘 − 1] have collectively at most 𝜔 nonzero coefficients.

    Output: A byte string 𝑦 of length 𝜔 + 𝑘 that encodes 𝐡 as described above
    """

    y = np.zeros(omega + k, dtype=int)
    index = 0
    for i in range(k):
        for j in range(256):
            if h[i][j] != 0:
                y[index] = j
                index = index + 1
        y[omega + i] = index
    return y

def HintBitUnpack(y):
    """
    Reverses the procedure HintBitPack

    Input: A byte string y of length ω + k that encodes h as described above

    Output: A polynomial vector 𝐡 ∈ 𝑅_2^k such or ⊥
    """

    h = np.zeros((k,256), dtype=int)
    index = 0
    for i in range(k):
        if y[omega + i] < index or y[omega + i] > omega:
            return "⊥"
        first = index
        while index < y[omega + i]:
            if index > first:
                if y[index - 1] >= y[index]:
                    return "⊥"
            h[i][y[index]] = 1
            index = index + 1
    for i in range(index, omega):
        if y[i] != 0:
            return "⊥"
    return h

def pkEncode(rho, t1):
    """ 
    Encodes a public key for ML-DSA into a byte string.

    Input: rho in B^32, t1 in Ring^k with coefficients in [0,2^bitlen(q-1)-d -1]

    Output: Public key pk in B^32+32k*bitlen(q-1)-d
    """
    bitlen_q = (q - 1).bit_length()
    pk = rho
    for i in range(k):
        pk.extend(SimpleBitPack(t1[i], 2**(bitlen_q - d) - 1))
    return pk

def pkDecode(pk):
    """
    Reverses the procedure pkEncode

    Input: Public key pk in B^32+32k*bitlen(q-1)-d

    Output: A tuple (rho, t1) where rho in B^32 and t1 in Ring^k with coefficients in [0,2^bitlen(q-1)-d -1]`
    """
    bitlen_q = (q - 1).bit_length()
    rho = pk[0:32]
    z = np.array(pk[32:]).reshape(k,32*(bitlen_q -d))
    t1 = np.zeros((k,256), dtype=int)
    for i in range(k):
        t1[i] = SimpleBitUnpack(z[i], 2**(bitlen_q - d) - 1)
    return (rho, t1)


#rho = [random.randint(0,255) for _ in range(32)]
#print(rho)
#bitlen_q = int(np.ceil(np.log2((q - 1))))
#t1 = [[random.randint(0,2**(bitlen_q - d) -1) for _ in range(256)] for _ in range(k)]
#pk = pkEncode(rho,t1)
#print(pkDecode(pk)[0])


def skEncode(rho, K, tr, s1, s2, t0):
    """
    Encodes a secrey key for ML-DSA into a byte string

    Input: 𝜌 ∈ 𝔹32, 𝐾 ∈ 𝔹32, 𝑡𝑟 ∈ 𝔹64 , 𝐬1 ∈ 𝑅ℓ with coefficients in [−𝜂, 𝜂], 𝐬2 ∈ 𝑅𝑘 with coefficients in [−𝜂, 𝜂], 𝐭0 ∈ 𝑅𝑘 with coefficients in [−2𝑑−1 + 1, 2𝑑−1].

    Output: Private key 𝑠𝑘 ∈ 𝔹32+32+64+32⋅((𝑘+ℓ)⋅bitlen (2𝜂)+𝑑𝑘)
    """
    sk = []
    sk.extend(rho)
    sk.extend(K)
    sk.extend(tr)

    for i in range(l):
        sk.extend(BitPack(s1[i], eta, eta))
    for i in range(k):
        sk.extend(BitPack(s2[i], eta, eta))
    for i in range(k):
        sk.extend(BitPack(t0[i], 2**(d-1) -1, 2**(d-1)))
    return sk

def skDecode(sk):
    """
    Reverses the procedure skEncode

    Input:  Private key 𝑠𝑘 ∈ 𝔹32+32+64+32⋅((ℓ+𝑘)⋅bitlen (2𝜂)+𝑑𝑘)

    Output:  𝜌 ∈ 𝔹32, 𝐾 ∈ 𝔹32, 𝑡𝑟 ∈ 𝔹64, 𝐬1 ∈ 𝑅ℓ , 𝐬2 ∈ 𝑅𝑘 , 𝐭0 ∈ 𝑅𝑘 with coefficients in [−2𝑑−1 + 1, 2𝑑−1].
    """
    bitlen_eta = (2*eta).bit_length()
    index_y_end = 128+int(32*bitlen_eta*l)
    index_z_end = index_y_end + int(32*bitlen_eta*k)

    rho = sk[0:32]
    K = sk[32:64]
    tr = sk[64:128]
    y = np.array(sk[128:index_y_end]).reshape((l,32*bitlen_eta))
    z = np.array(sk[index_y_end:index_z_end]).reshape((k,32*bitlen_eta))
    w = np.array(sk[index_z_end:]).reshape((k,32*d))

    for i in range(l):
        s1[i] = BitUnpack(y[i], eta, eta)
    for i in range(k):
        s2[i] = BitUnpack(z[i], eta, eta)
    for i in range(k):
        t0[i] = BitUnpack(w[i], 2**(d-1) -1, 2**(d-1))
    return (rho, K, tr, s1, s2, t0)


def generate_ring_vector(L, n, a, b):
    """
    Generates a polynomial vector of length L, with coefficents in range [a,b] with n coeffcients

    Input: a,b,n,L in N

    Output: vector in Ring^L
    """
    return [[random.randint(a,b) for _ in range(n)] for _ in range(L)]

def random_sparse_vector(k: int, n: int = 10, omega: int = 10):
    """
    Generate a vector of length k containing binary polynomials (length n).
    The total number of 1 coefficients across all polynomials is at most `omega`.
    """
    # total available positions = k * n
    total_positions = k * n
    
    # how many ones to place (≤ omega, but cannot exceed total_positions)
    num_ones = np.random.randint(0, min(omega, total_positions))
    
    # choose distinct positions for ones
    one_positions = set(random.sample(range(total_positions), num_ones))
    
    # build the vector
    vector = []
    for i in range(k):
        poly = [0] * n
        for j in range(n):
            if i * n + j in one_positions:
                poly[j] = 1
        vector.append(poly)
    
    return vector

s1 = generate_ring_vector(l, 256, -eta, eta)
s2 = generate_ring_vector(k, 256, -eta, eta)
t0 = generate_ring_vector(k, 256, -2**(d-1) + 1, 2**(d-1))
rho = [random.randint(0,256) for _ in range(32)]
K = [random.randint(0,256) for _ in range(32)]
tr = [random.randint(0,256) for _ in range(64)]

#sk = skEncode(rho, K, tr, s1, s2, t0)
sk_decode = skDecode(skEncode(rho, K, tr, s1, s2, t0))
sk_correct = [rho, K, tr, s1, s2, t0]


def sigEncode(c_tilde, z, h):
    """
    Encodes a signture with a byte string

    Input:c ̃ ∈ Bλ/4, z ∈ Rl with coefficients in [−γ1 + 1, γ1], h ∈ R2k.
    
    Output: Signature σ ∈ Bλ/4+l⋅32⋅(1+bitlen (γ1−1))+ω+k.
    """

    sigma = c_tilde.copy()
    for i in range(l):
        sigma.extend(BitPack(z[i], gamma1 -1, gamma1))

    sigma.extend(HintBitPack(h))
    return sigma

def sigDecode(sigma):
    """
    Reverses the procedure sigEncode

    Input: Signature σ ∈ Bλ/4+l⋅32⋅(1+bitlen (γ1−1))+ω+k.

    Output: c ̃ ∈ Bλ/4, z ∈ Rl with coefficients in [−γ1 + 1, γ1], h ∈ R2k, or ⊥.
    """
    c_end = int(lam/4)
    x_end = c_end + 32*(1+(gamma1 -1).bit_length())*l
    c_tilde = sigma[0:c_end]
    x = np.asarray(sigma[c_end:x_end]).reshape((l,32*(1+(gamma1 -1).bit_length())))
    y = sigma[x_end:]

    z = np.zeros((l,256), dtype=int)
    for i in range(l):
        z[i] = BitUnpack(x[i], gamma1 -1, gamma1)
    h = HintBitUnpack(y)
    return c_tilde, z, h


def w1Encode(w1):
    """
    Encodes a polynomial vector w1 into a byte string
    
    Input: w1 ∈ Rk whose polynomial coordinates have coefficients in [0, (q − 1)/(2γ2) − 1].

    Output: A byte string representing the encoded polynomial vector.A byte string representation ̃ .
    """

    w1_tilde = []
    for i in range(k):
        w1_tilde.extend(SimpleBitPack(w1[i], int((q - 1)/(2*gamma2) - 1)))
    return w1_tilde


#c_tilde = [random.randint(0,255) for _ in range(int(lam/4))]
#z = generate_ring_vector(l, 256, -gamma1 + 1, gamma1)
#h = random_sparse_vector(k, 256, omega)
#sigma = sigEncode(c_tilde, z, h)
#sig = sigDecode(sigma)


def SampleInBall(rho):
    """
    Samples a polynomial c ∈ R with coefficients from {−1, 0, 1} and Hamming weight τ ≤ 64.

    Input: A seed rho in B^lambda/4

    Output: A polynomial c in Ring
    """
    c = np.zeros(256, dtype=int)
    ctx = SHAKE256.new()    
    ctx.update(bytearray(rho))
    s = list(ctx.read(8))
    h = BytesToBits(s)
    for i in range(256 - tau, 256):
        j = list(ctx.read(1))[0]
        while j > i:
            j = list(ctx.read(1))[0]
        c[i] = c[j]
        c[j] = (-1)**h[i + tau - 256]
    return c

#rho = [random.randint(0,255) for _ in range(int(lam/4))]
#print(SampleInBall(rho))


def RejNTTPoly(rho):
    """
    Samples a polynomial ∈ Tq.

    Input: A seed ρ ∈ B34. 
    
    Output: An element â ∈ Tq.
    """
    a_tilde = np.zeros(256, dtype=int)
    j = 0
    ctx = SHAKE128.new()
    ctx.update(bytearray(rho))
    while j < 256:
        s = list(ctx.read(3))
        a_tilde[j] = CoeffFromThreeBytes(s[0], s[1], s[2])
        if a_tilde[j] != "⊥":
            j = j+1
    return a_tilde


#rho = [random.randint(0,255) for _ in range(34)]
#print(RejNTTPoly(rho))

def RejBoundedPoly(rho):
    """
    Samples an element a ∈ R with coefficients in [−η, η] computed via rejection sampling from ρ.  
    
    Input: A seed ρ ∈ B66. 
    
    Output: A polynomial a ∈ R.
    """
    a = np.zeros(256, dtype=int)
    j = 0
    ctx = SHAKE256.new()
    ctx.update(bytearray(rho))
    while j < 256:
        z = list(ctx.read(1))[0]
        z0 = CoeffFromHalfByte(z%16)
        z1 = CoeffFromHalfByte(int(np.floor(z/16)))
        if z0 != "⊥":
            a[j] = z0
            j = j+1
        if z1 != "⊥" and j < 256:
            a[j] = z1
            j = j+1
    return a

#rho = [random.randint(0,255) for _ in range(66)]
#print(RejBoundedPoly(rho)) 

def ExpandA(rho):
    """
    Samples a k × l matrix Â of elements of Tq.  

    Input: A seed ρ ∈ B32.

    Output: Matrix A ∈ (Tq)k×l.
    """
    A_hat = np.zeros((k,l,256), dtype=int) 
    for r in range(k):
        for s in range(l):
            rho_dash = rho.copy()
            rho_dash.extend(IntegerToBytes(s,1))
            rho_dash.extend(IntegerToBytes(r,1))
            A_hat[r,s] = RejNTTPoly(rho_dash)
    return A_hat

#rho = [random.randint(0,255) for _ in range(32)]
#A_hat = ExpandA(rho)
#print(A_hat)


def ExpandS(rho):
    """
    Samples vectors s1 ∈ Rl and s2 ∈ Rk, each with polynomial coordinates whose coefficients are in the interval [−η, η].  
    
    Input: A seed ρ ∈ B64. 
    
    Output: Vectors s1, s2 of polynomials in R.
    """
    s1 = np.zeros((l,256), dtype=int)
    s2 = np.zeros((k,256), dtype=int)
    for r in range(l):
        rho_1 = rho.copy()
        rho_1.extend(IntegerToBytes(r,2))
        s1[r] = RejBoundedPoly(rho_1)
    for r in range(k):
        rho_2 = rho.copy()
        rho_2.extend(IntegerToBytes(r+l,2))
        s2[r] = RejBoundedPoly(rho_2)
    return s1, s2

#rho = [random.randint(0,255) for _ in range(64)]
#s1, s2 = ExpandS(rho)
#print(s1)
#print(s2)

def ExpandMask(rho, mu):
    """
    Samples a vector y ∈ Rl such that each polynomial y[r] has coefficients between −γ1 + 1 and γ1.  
    
    Input: A seed ρ ∈ B64 and a nonnegative integer μ. 
    
    Output: Vector y ∈ Rl.
    """
    y = np.zeros((l, 256), dtype=int)
    c = 1 + (gamma1 - 1).bit_length()
    for r in range(l):
        rho_dash = []
        rho_dash.extend(rho)
        rho_dash.extend(IntegerToBytes(mu + r, 2))
        ctx = SHAKE256.new()
        ctx.update(bytearray(rho_dash))
        v = list(ctx.read(32*c))
        y[r] = BitUnpack(v, gamma1 -1, gamma1)
    return y


#rho = [random.randint(0,255) for _ in range(64)]
#mu = random.randint(0,255)
#y = ExpandMask(rho, mu)
#print(y)

def mod_pm(m: int, alpha: int) -> int:
    """
    Compute m' = m mod^{±} alpha, i.e. the signed residue of m in the interval −⌈α/2⌉ < m' ≤ ⌊α/2⌋  such that m and m' are congruent modulo alpha.
    """
    m_dash = m%alpha
    if m_dash > np.floor(alpha/2):
        m_dash = m_dash - alpha
    return m_dash


#r = 10
#print(mod_pm(r,3))

def Power2Round(r):
    """
    Decomposes r into (r1, r0) such that r ≡ r12d + r0 mod q.  
    
    Input: r ∈ Zq.  
    
    Output: Integers (r1, r0).
    """
    r_plus = r%q
    r0 = mod_pm(r_plus,2**d)
    return (r_plus - r0)/(2**d), r0

#r = random.randint(0,q-1)
#r1, r0 = Power2Round(r)
#print(r, r0, r1)

def Decompose(r):
    """
    Decomposes r into (r1, r0) such that r ≡ r1(2γ2) + r0 mod q.  
    
    Input: r ∈ Zq.  
    
    Output: Integers (r1, r0).
    """
    r_plus = r%q
    r0 = mod_pm(r_plus, int(2*gamma2))
    if r_plus - r0 == q-1:
        r1 = 0
        r0 = r0 - 1
    else:
        r1 = (r_plus - r0)/(2*gamma2)
    return r1, r0

#r = random.randint(0,q-1)
#r1, r0 = Decompose(r)
#print(r, r0, r1)

def HighBits(r):
    """
    Computes the high bits of r.  
    
    Input: r ∈ Zq.  
    
    Output: Integer r1.
    """
    r1, r0 = Decompose(r)
    return r1

def LowBits(r):
    """
    Computes the low bits of r.  
    
    Input: r ∈ Zq.  
    
    Output: Integer r0.
    """
    r1, r0 = Decompose(r)
    return r0

def MakeHint(z,r):
    """
    Computes hint bit indicating whether adding z to r alters the high bits of r.  
    
    Input: z, r ∈ Zq. 
    
    Output: Boolean.
    """
    r1 = HighBits(r)
    v1 = HighBits(r+z)
    return r1 != v1

#r = random.randint(0,q-1)
#z = random.randint(-gamma1 + 1, gamma1)
#print(MakeHint(z,r))

def UseHint(r, h):
    """
    Returns the high bits of r adjusted according to hint h.  
    
    Input: Boolean h, r ∈ Zq.  
    
    Output: r1 ∈ Z with 0 ≤ r1 ≤ q−1 .
    """
    m = (q-1)/(2*gamma2)
    r1, r0 = Decompose(r)
    if h == 1 and r0 > 0:
        return (r1 + 1)%m
    if h == 1 and r0 <= 0:
        return (r1 - 1)%m
    else:
        return r1

#r = random.randint(0,q-1)
#z = random.randint(-gamma1 + 1, gamma1)
#h = MakeHint(z,r)
#r1 = UseHint(r,h)
#print(r, z, h, r1)


def NTT(w):
    """
    Computes the NTT.  
    
    Input: Polynomial w(X) = ∑255  j=0 wjXj ∈ Rq.  
    
    Output: ŵ = (ŵ[0], ... , ŵ[255]) ∈ Tq.
    """
    w_hat = np.zeros(256, dtype=int)
    for j in range(256):
        w_hat[j] = w[j]
    m = 0
    len = 128
    while len >= 1:
        start = 0
        while start < 256:
            m = m+1
            z = zetas[m]
            for j in range(start, start + len):
                t = z*w_hat[j + len]%q
                w_hat[j+len] = (w_hat[j] - t)%q
                w_hat[j] = (w_hat[j] + t)%q
            start = start + 2*len
        len = int(np.floor(len/2))
    return w_hat

#w = [random.randint(0,q-1) for _ in range(256)]
#print(w, NTT(w))

def InvNTT(w_hat):
    """
    Computes the inverse of the NTT.  
    
    Input: ŵ = (ŵ[0], ... , ŵ[255]) ∈ Tq.  
    
    Output: Polynomial w(X) = ∑255  j=0 wjXj ∈ Rq.
    """
    w = np.zeros(256, dtype=int)
    for j in range(256):
        w[j] = w_hat[j]
    m = 256
    len = 1
    while len < 256:
        start = 0
        while start < 256:
            m = m-1
            z = -zetas[m]
            for j in range(start, start + len):
                t = w[j]
                w[j] = (t + w[j + len])%q
                w[j + len] = (t - w[j + len])%q
                w[j + len] = (z*w[j + len])%q
            start = start + 2*len
        len = 2*len
    f = 8347681
    for j in range(256):
        w[j] = f* w[j]%q
    return w

#w = [random.randint(0,q-1) for _ in range(256)]
#w_hat = NTT(w)
#w_return = InvNTT(w_hat)
#print(w == w_return)

def BitRev8(m):
    """
    Transforms a byte by reversing the order of bits in its 8-bit binary expansion.  
    
    Input: A byte m ∈ [0, 255]. 
    
    Output: A byte r ∈ [0, 255].
    """
    b = IntegerToBits(m, 8)
    b_rev = np.zeros(8, dtype=int)
    for i in range(8):
        b_rev[i] = b[7-i]
    r = BitsToInteger(b_rev, 8)
    return r

#print(BitRev8(171))  ## should be 213


def AddNTT(a_hat, b_hat):
    """
    Computes the sum â + b̂ of two elements â, b̂ ∈ Tq.
    
    Input: â, b ∈ Tq. 
    
    Output: c ̂ ∈ Tq
    """
    c_hat = np.zeros(256, dtype=int)
    for j in range(256):
        c_hat[j] = (a_hat[j] + b_hat[j])%q
    return c_hat


def MultiplyNTT(a_hat, b_hat):
    """
    Computes the product â · b̂ of two elements â, b̂ ∈ Tq.
    
    Input: â, b ∈ Tq. 
    
    Output: c ̂ ∈ Tq
    """
    c_hat = np.zeros(256, dtype=int)
    for j in range(256):
        c_hat[j] = (a_hat[j] * b_hat[j])%q
    return c_hat

def AddVectorNTT(v_hat, w_hat):
    """
    Computes the sum v̂ + ŵ of two vectors v̂, ŵ over Tq.  
    
    Input: l ∈ N, v̂ ∈ Tql, ŵ ∈ Tql.  
    
    Output: û ∈ Tql.
    """
    u_hat = np.zeros((l,256), dtype=int)
    for i in range(l):
        u_hat[i] = AddNTT(v_hat[i], w_hat[i])
    return u_hat

def ScalarVectorNTT(c_hat, v_hat):
    """
    Computes the product c ̂∘ v̂ of a scalar c ̂and a vector v̂ over Tq.  
    
    Input: c ̂ ∈ Tq, l ∈ N, v̂ ∈ Tql.  
    
    Output: ŵ ∈ Tql.
    """
    w_hat = np.zeros((l,256), dtype=int)
    for i in range(l):
        w_hat[i] = MultiplyNTT(c_hat, v_hat[i])
    return w_hat

def MatrixVectorNTT(M_hat, v_hat):
    """
    Computes the product M̂ ∘ v̂ of a matrix M̂ and a vector v̂ over Tq.  
    
    Input: k, l ∈ N, M ∈ Tqk×l ̂ q .  ̂ ,v ∈ Tl  
    
    Output: ŵ ∈ Tqk.
    """
    w_hat = np.zeros((k,256), dtype=int)
    for i in range(k):
        for j in range(l):
            w_hat[i] = AddNTT(w_hat[i], MultiplyNTT(M_hat[i,j], v_hat[j]))
    return w_hat

#print(BytesToBits(np.array([171, 171])))   
#print(BitsToBytes(BytesToBits(np.array([171, 170]))))
#print(CoeffFromThreeBytes(128, 0, 128))
#print(CoeffFromThreeBytes(255, 255, 255))
#print(CoeffFromHalfByte(1))
#w = [random.randint(0,15) for i in range(255)]
#print(SimpleBitPack(np.array(w), 15))
#w = [random.randint(-15,15) for i in range(255)]
#print(BitPack(np.array(w), 15, 15))
#w = [random.randint(0,15) for i in range(255)]
#print(SimpleBitUnpack(SimpleBitPack(np.array(w),15),15))
#print(w)
#w = np.array([random.randint(-15,15) for i in range(256)])
#print(BitUnpack(BitPack(w,15,15),15,15))
#print(w)




#h = random_sparse_vector(k=k,n=256,omega=omega)
#print(HintBitPack(h))
#h_dash = HintBitUnpack(HintBitPack(h))
#l=0
