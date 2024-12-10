from lab4.src.const import ip_permutations, f_perm, pc1_permutations, e_bit_permutation, pc2_permutations, p_permutations, \
    s_blocks




def ip_transmutation(data: str):
    if len(data) != 64:
        raise ValueError("Invalid input length for IP transmutation. Expected 64 bits.")
    result = "".join(data[i - 1] for i in ip_permutations)
    return result

def f_transmutation(data: str):
    if len(data) != 64:
        raise ValueError("Invalid input length for IP transmutation. Expected 64 bits.")
    result = "".join(data[i - 1] for i in f_perm)
    return result

def pc1_transmutation(data: str):
    if len(data) != 64:
        raise ValueError("Invalid input length for PC1 transmutation. Expected 64 bits.")
    result = "".join(data[i - 1] for i in pc1_permutations)
    return result

def e_bit_expansion(data: str):
    if len(data) != 32:
        raise ValueError("Invalid input length for E-bit expansion. Expected 32 bits.")
    result = "".join(data[i - 1] for i in e_bit_permutation)
    return result

def bit_addition(first: str, second: str):
    if len(first) != len(second) or len(first) != 32:
        raise ValueError("Invalid input length for bit addition. Expected 32-bit strings.")
    int1 = int(first, 2)
    int2 = int(second, 2)
    result = (int1 + int2) % (1 << 32)
    return bin(result)[2:].zfill(32)

def pc2_transmutation(data: str):
    if len(data) != 56:
        return ValueError("Invalid input length for PC2 transmutation. Expected 56 bits.")
    result = ""
    for i in pc2_permutations:
        result += data[i - 1]

    return result

def left_shift(data: str):
    if len(data) != 28:
        return "not valid length"
    first_char = data[0]
    result = data[1:] + first_char
    return result

def xor(first: str, second: str):
        if len(first) != len(second) :
            raise ValueError("Invalid input length for XOR operation. Expected 48 bits.")
        return ''.join('1' if bit1 != bit2 else '0' for bit1, bit2 in zip(first, second))

def func_rk(r: str, k: str):

    def p_transmutation(data: str):
        if len(data) != 32:
            raise ValueError("Invalid input length for P transmutation. Expected 32 bits.")
        return ''.join(data[i - 1] for i in p_permutations)

    expanded_r = e_bit_expansion(r)
    b = xor(expanded_r, k)

    blocks = [b[i:i + 6] for i in range(0, len(b), 6)]
    if len(blocks) != 8:
        raise ValueError("Invalid number of blocks. Expected 8 blocks of 6 bits each.")

    s_block = ""
    for index, block in enumerate(blocks):
        if len(block) != 6:
            raise ValueError(f"Block at index {index} has invalid length: {len(block)}.")
        i = int(block[0] + block[5], 2)  # Row bits
        j = int(block[1:5], 2)          # Column bits
        s = s_blocks[index][i][j]
        s_block += bin(s)[2:].zfill(4)

    return p_transmutation(s_block)

def des_encrypt(m:str, k:str)->str:
    ip = ip_transmutation(m)
    k_plus = pc1_transmutation(k)

    Cn = k_plus[:28]
    Dn = k_plus[28:]
    Ln = ip[:32]
    Rn = ip[32:]

    for i in range(16):
        if i not in [0, 1, 8, 15]:
            Cn = left_shift(Cn)
            Dn = left_shift(Dn)

        ln = Ln
        rn = Rn
        Cn = left_shift(Cn)
        Dn = left_shift(Dn)

        CnDn = Cn + Dn
        Kn = pc2_transmutation(CnDn)

        Ln = rn
        rk = func_rk(rn, Kn)
        Rn = xor(ln, rk)

        # print(f"C{i}=" + Cn + f"\tD{i}=" + Dn)
        # print(f"L{i+1}=" + Ln + f"\tR{i+1}=" + Rn)
        # print(f"K{i+1}=" + Kn + "\n")

    crp = f_transmutation(Rn + Ln)
    return crp

