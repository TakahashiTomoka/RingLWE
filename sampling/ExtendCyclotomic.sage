###
# I will try to sample from the field bQ(zeta_p, \sqrt{d})
# and perform the degree two attack on it.
##
load('sampling/mega.sage')

def new_t_matrix(n, prec = 100):
    """
    return the unitary matrix 1/sqrt{2} (1,1,-i, i) \otimes I.
    """
    C = ComplexField(prec)
    eyer = Matrix.identity(C, n//2)

    return C(1/sqrt(2))*block_matrix([[eyer, eyer], [-C(I)*eyer, C(I)*eyer]])



class ExtendCyclotomic:
    """
    The sampler for the supposedly very vulnerable field Q(zeta_p, sqrt(d))
    """
    def __init__(self, p, d, f, sigma0 = 1, prec = 100):
        # checks
        if not ZZ(p).is_prime():
            raise ValueError('p must be prime.')
        #elif d <= 0 or not ZZ(d).is_squarefree or Mod(d,4) == 1 or Mod(d,p) == 0:
        #    raise ValueError('d does not meet criteria')
        self.p = p
        self.d = d
        self.f = f
        self.n = f*(p-1)
        self.q = None

        # computing matrices
        self.sigma0 = sigma0
        self.C = ComplexField(prec)
        self.R = RealField(prec)
        self.Bc = self.compute_Bc()
        self.T = new_t_matrix(self.n)
        #print ('T computed')
        B0 = _real_part(self.T*self.Bc)
        #print(' B0 computed')
        self.B = B0
        #self.B = block_matrix([[B0, B0*self.R(d**(0.5))], [B0,-B0*self.R(d**(0.5))]])

        # discriminant
        self.dK = self.compute_dK()
        self.D = MyLatticeSampler(self.B, sigma = self.sigma0, method = None)

    def __repr__(self):
        return 'RLWE error sampler on the field Q(zeta_%s, %s^(1/%s))'%(self.p, self.d, self.f)

    def sigma0(self):
        return self.sigma0

    def setq(self,q):
        p, d = self.p, self.d
        if not ZZ(q).is_prime() or Mod(q,p) != 1 or kronecker(d,q) != -1:
            raise ValueError('q does not meet all requirements')
        self.q = q
        self.k = GF(q)
        PolyRing.<x> = self.k[];
        self.F = self.k.extension(x^self.f-self.d,'a')

        # a p-th root in Fq
        self.alpha = cyclotomic_polynomial(p).change_ring(self.k).roots(multiplicities = False)[0]



    def getq(self):
        return self.q

    def compute_Bc(self):
        # A_{i,j} = A_{i,(p-1)/2 + j}^{-1}
        p = self.p
        f = self.f
        d = self.d
        C = self.C
        zetap = C.zeta(p)
        zeta2f = C.zeta(2*f)
        drootf = C(d^(1/f))
        #tt = Matrix(C, [[zetap^(i*j) for j in range(p-1)] for i in range((p-1)//2)])
        #bb = Matrix(C, [[zetap^((-i)*j) for j in range(p-1)] for i in range((p-1)//2)])
        tt = Matrix(C, [[zetap^(i*j) for j in range(1, p)] for i in range(1, (p-1)//2 + 1)])
        bb = Matrix(C, [[zetap^((-i)*j) for j in range(1, p)] for i in range(1, (p-1)//2 + 1)])

        top = [[tt*drootf^j*zeta2f^((2*i+1)*j) for j in range(f)] for i in range(f)]
        bot = [[bb*drootf^j*zeta2f^((-2*i-1)*j) for j in range(f)] for i in range(f)]
        return block_matrix(C, top+bot)

        #Index = _index_matrix(p)
        #return  Matrix(C, [[zetap**Index[i][j] for j in range(p-1)] for i in range(p-1)])

    def compute_dK(self):
        """
        compute the discriminant. For a sanity check.
        """
        return (self.p)**(2*(self.p-2))*(4*self.d)**(self.p-1)

    def __call__(self, reduced = True):
        if reduced:
            return self.reduce_modq(self.D()[1])
        else:
            return self.D()[1]

    def reduce_modq(self, vector_ZZ):
        if self.q is None:
            raise ValueError('q is unset.')
        q = self.q
        k = self.k
        F = self.F
        sqrtd = F.gen(0)
        alpha = self.alpha
        v = vector_ZZ
        v1,  v2 = v[:len(v)//2], v[len(v)//2:]
        assert len(v1) == len(v2)
        return self.F(sum([v1[i]*alpha**i for i in range(len(v1))])) + self.F(sum([v2[i]*alpha**i for i in range(len(v2))]))*sqrtd



