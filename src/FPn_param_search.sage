p = 13
f = 9
bound_q = 100
bound_d = 1000

print("p: {}".format(p))

for q in range(2*p+1, bound_q, p*2):#qも奇数だから、p*2づつ増加
    if not is_prime(q): continue
#    if q%f != 1: continue
    print("q: {}".format(q))

    
    Root = []
    for i in range(q):
        root = (i^f)%q
        if root not in Root:
            Root.append(root)#f乗根が存在するもの
    
    for d in range(bound_d):
        if d%q == 0: continue
        d_mod = d%q
        if d_mod not in Root:
            if d^(1/f) not in ZZ:
                print("q: {0}, d: {1}".format(q, d))
                #with open('Fp3_param.txt', mode ='a+',encoding = 'utf-8') as f:
                #    f.write("p: {}, q: {}, d: {}\n".format(p, q, d))




