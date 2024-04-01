class HypProblem:
    def __init__(self, T, S, C, B, F, G, X0, Y0) -> None:
        self.set_T(T)
        self.set_S(S)
        self.set_C(C)
        self.set_B(B)
        self.set_F(F)
        self.set_G(G)
        self.set_X0(X0)
        self.set_Y0(Y0)

    def set_T(self, T):
        self.T0 = T[0]
        self.T1 = T[1]


    def set_S(self, S):
        self.S0 = S[0]
        self.S1 = S[1]

    def set_C(self, C):
        self.C1 = C[0]
        self.C2 = C[1]

    def set_B(self, B):
        self.B11 = B[0][0]
        self.B12 = B[0][1]
        self.B21 = B[1][0]
        self.B22 = B[1][1]

    def set_F(self, F):
        self.F1 = F[0]
        self.F2 = F[1]

    def set_G(self, G):
        self.G11 = G[0][0]
        self.G12 = G[0][1]
        self.G21 = G[1][0]
        self.G22 = G[1][1]

    def set_X0(self, X0):
        self.x0 = X0

    def set_Y0(self, Y0):
        self.y0 = Y0

class ConjHypProblem(HypProblem):
    def __init__(self, hyp_problem: HypProblem):
        self.B11 = lambda s, t: -hyp_problem.B11(hyp_problem.S1-s, hyp_problem.T1-t)
        self.B12 = lambda s, t: -hyp_problem.B12(hyp_problem.S1-s, hyp_problem.T1-t)
        self.B21 = lambda s, t: -hyp_problem.B21(hyp_problem.S1-s, hyp_problem.T1-t)
        self.B22 = lambda s, t: -hyp_problem.B22(hyp_problem.S1-s, hyp_problem.T1-t)

        self.G11 = lambda t: -hyp_problem.G11(hyp_problem.T1-t)
        self.G12 = lambda t: -hyp_problem.G12(hyp_problem.T1-t)
        self.G21 = lambda t: -hyp_problem.G21(hyp_problem.T1-t)
        self.G22 = lambda t: -hyp_problem.G22(hyp_problem.T1-t)

        self.F1 = lambda s, t: -hyp_problem.F1(hyp_problem.S1-s, hyp_problem.T1-t)
        self.F2 = lambda s, t: -hyp_problem.F2(hyp_problem.S1-s, hyp_problem.T1-t)

        self.x0 = hyp_problem.x0
        self.y0 = hyp_problem.y0

        self.set_T([hyp_problem.T0,hyp_problem.T1])
        self.set_S([hyp_problem.S0,hyp_problem.S1])
        self.set_C([hyp_problem.C1,hyp_problem.C2])