import numpy as np
import os
import ML_DSA_internal as internal


def KeyGen():
    """
    Generates a public-private key pair

    Output: Public key, pk and private key, sk
    """
    zeta = os.urandom(32)
    if zeta == None:
        raise ValueError("Failed to generate random bytes")
    return internal.KeyGen_internal(zeta)


def Sign(sk, M, ctx):
    """
    Generates an ML-DSA signature

    Inputs private key, sk, message, M and context string, ctx

    Outputs signature sigma
    """
    if np.abs(ctx) > 255:
        raise ValueError("Context string is too long")

    rnd = os.urandom(32)
    if rnd == None:
        raise ValueError("Failed to generate random bytes")
    
    M_dash = internal.BytesToBits(internal.IntegerToByes(0,1) + internal.IntegerToBytes(np.abs(ctx),1) + ctx) + M
    sigma = internal.Sign_internal(sk, M_dash, rnd)
    return sigma