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


rho = [random.randint(0,255) for _ in range(64)]
mu = random.randint(0,255)
y = ExpandMask(rho, mu)
print(y)
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
