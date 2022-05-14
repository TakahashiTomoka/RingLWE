f = 9
bound_q = 500

for p in range(100):
    if not is_prime(p): continue

    for q in range(2*p+1, bound_q, p*2):#qも奇数だから、p*2づつ増加
        if not is_prime(q): continue
        if q%f != 1: continue
        print("p: {}, q: {}".format(p,q))

    



