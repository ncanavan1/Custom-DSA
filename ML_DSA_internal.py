import numpy as np
import random
import os
from Crypto.Hash import SHAKE256
from Crypto.Hash import SHAKE128
import auxillary as aux


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

def KeyGen_internal(zeta):
    """
    Generates a public-private key pair from a seed.  
    
    Input: Seed ξ ∈ B32  
    
    Output: Public key pk ∈ B32+32k(bitlen (q−1)−d)  and private key sk ∈ B32+32+64+32⋅((l+k)⋅bitlen (2η)+dk).
    """

    zeta.extend(aux.IntegerToBytes(k,1))
    zeta.extend(aux.IntegerToBytes(l,1))
    ctx = SHAKE256.new()
    ctx.update(bytearray(zeta))
    H_op = list(ctx.read(128))

    rho = H_op[0:32]
    rho_dash = H_op[32:96]
    K = H_op[96:128]

    A_hat = aux.ExpandA(rho)
    s1, s2 = aux.ExpandS(rho_dash)
    t = aux.ListInvNTT(aux.MatrixVectorNTT(A_hat, aux.ListNTT(s1))) + s2
    t1, t0 = aux.Power2RoundVec(t)
    pk = np.asarray(aux.pkEncode(rho, t1))
    
    ctx = SHAKE256.new()
    ctx.update(bytearray(pk))
    tr = list(ctx.read(64))

    sk = np.asarray(aux.skEncode(rho, K, tr, s1, s2, t0))

    return pk, sk

zeta = [random.randint(0,255) for _ in range(32)]
pk, sk = KeyGen_internal(zeta)
print("Public Key: ", pk.tolist())
print("Secret Key: ", sk.tolist())

