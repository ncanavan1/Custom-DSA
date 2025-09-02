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

#zeta = [random.randint(0,255) for _ in range(32)]
#pk, sk = KeyGen_internal(zeta)
#print("Public Key: ", pk.tolist())
#print("Secret Key: ", sk.tolist())

def Sign_internal(sk, M_dash, rnd): 
    """
    Deterministic algorithm to generate a signature for a formatted message M ′.  
    
    Input: Private key sk ∈ B32+32+64+32⋅((l+k)⋅bitlen (2η)+dk), formatted message M ′ ∈ {0, 1}∗, and  per message randomness or dummy variable rnd ∈ B32.  
    
    Output: Signature σ ∈ Bλ/4+l⋅32⋅(1+bitlen (γ1−1))+ω+k.
    """

    rho, K, tr, s1, s2, t0 = aux.skDecode(sk)
    s1_hat = aux.ListNTT(s1)
    s2_hat = aux.ListNTT(s2)
    t0_hat = aux.ListNTT(t0)
    A_hat = aux.ExpandA(rho)
    
    ctx = SHAKE256.new()
    inter = list(aux.BytesToBits(tr))
    inter.extend(M_dash)
    ctx.update(bytearray(inter))
    mu = list(ctx.read(64))

    ctx = SHAKE256.new()
    inter = K.copy()
    inter.extend(rnd)
    inter.extend(mu)
    ctx.update(bytearray(inter))
    rho_dash_dash = list(ctx.read(64))

    kappa = 0
    z = "⊥"
    h = "⊥"

    while z == "⊥" or h == "⊥": ##maybe and??
        y = aux.ExpandMask(rho_dash_dash, kappa)
        w = aux.ListInvNTT(aux.MatrixVectorNTT(A_hat, aux.ListNTT(y)))
        w1 = aux.HighBitsVec(w) ##need componet wise version
        
        ctx = SHAKE256.new()
        inter = mu.copy()
        inter.extend(aux.w1Encode(w1))
        ctx.update(bytearray(inter))

        c_tilde = list(ctx.read(int(lam/4)))
        c = aux.SampleInBall(c_tilde)
        c_hat = aux.NTT(c)

        cs1 = aux.ListInvNTT([aux.MultiplyNTT(c_hat, s1_hat_vec) for s1_hat_vec in s1_hat])
        cs2 = aux.ListInvNTT([aux.MultiplyNTT(c_hat, s2_hat_vec) for s2_hat_vec in s2_hat])

        z = y + cs1

        r0 = aux.LowBitsVec(w - cs2)  ###need component wise version
        
        if np.max(np.abs(z)) >= gamma1 - beta or np.max(np.abs(r0)) >= gamma2 - beta:
            print(np.max(np.abs(z)), np.max(np.abs(r0)), kappa)
            z = "⊥"
            h = "⊥"
        else:
            ct0 = aux.InvNTT(aux.MultiplyNTT(c_hat, t0_hat))
            h = aux.MakeHint(-ct0, w - cs2 + ct0)
            if np.max(np.abs(ct0)) >= gamma2 or np.sum(h) > omega: ##must be all ones in sum
                z = "⊥"
                h = "⊥"
        kappa = kappa + l
    sigma  = aux.sigEncode(c_tilde, aux.mod_pm(z,q),h)
    return sigma



sk_len = 32 + 32 + 64 + 32 * ((l + k) * (2 * eta).bit_length() + d*k)
sk = [random.randint(0,255) for _ in range(sk_len)]

rnd = [random.randint(0,255) for _ in range(32)]

M_dash = [random.randint(0,1) for _ in range(500000)]

sigma = Sign_internal(sk, M_dash, rnd)

