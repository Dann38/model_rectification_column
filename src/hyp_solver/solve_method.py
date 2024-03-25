from .problem import HypProblem
from typing import Tuple
import numpy as np

class SolveMethod:
    def __init__(self) -> None:
        pass

    def start_solve(self, problem: HypProblem, s:float, t:float,
                    sl:float, tl:float, xl:float, yl:float, hl:float, 
                    sr:float, tr:float, xr:float, yr:float, hr:float) -> Tuple[float, float]:

        B11, B12, B21, B22 = problem.B11, problem.B12, problem.B21, problem.B22,
        F1, F2 = problem.F1, problem.F2

        A = [[1 - hr/2*B11(s, t), -hr/2*B12(s, t)],
             [-hl/2*B21(s, t), 1-hl/2*B22(s, t)]]
        b = [xr + hr/2*(B11(sr, tr)*xr+B12(sr, tr)*yr + F1(s, t) + F1(sr, tr)),
             yl + hl/2*(B21(sl, tl)*xl+B22(sl, tl)*yl + F2(s, t) + F2(sl, tl))]
        return np.linalg.solve(A, b)


    def left_solve(self, problem: HypProblem, s:float, t:float,
                   sl:float, tl:float, xl:float, yl:float, hl:float, 
                   sr:float, tr:float, xr:float, yr:float, hr:float) -> Tuple[float, float]:
        
        B11 = problem.B11
        B12 = problem.B12
        F1 = problem.F1
        G21 = problem.G21
        G22 = problem.G22

        A = [[1 - hr / 2 * B11(s, t), -hr / 2 * B12(s, t)],
                [-hl / 2 * G21(t), 1 - hl / 2 * G22(t)]]
        b = [xr + hr / 2 * (B11(sr, tr) * xr + B12(sr, tr) * yr + F1(s, t) + F1(sr, tr)),
                yl + hl / 2 * (G21(tl) * xl + G22(tl) * yl)]
        return np.linalg.solve(A, b)
    
    
    def right_solve(self, problem: HypProblem, s:float, t:float,
                   sl:float, tl:float, xl:float, yl:float, hl:float, 
                   sr:float, tr:float, xr:float, yr:float, hr:float) -> Tuple[float, float]:
        
        B21 = problem.B21
        B22 = problem.B22
        F2 = problem.F2
        G12 = problem.G12
        G11 = problem.G11

        A = [[1 - hr / 2 * G11(t), -hr / 2 * G12(t)],
                [-hl / 2 * B21(s, t), 1 - hl / 2 * B22(s, t)]]
        b = [xr + hr / 2 * (G11(tr) * xr + G12(tr) * yr),
                yl + hl / 2 * (B21(sl, tl) * xl + B22(sl, tl) * yl + F2(s, t) + F2(sl, tl))]
        return np.linalg.solve(A, b)
    
    def finish_solve(self, problem: HypProblem, s:float, t:float,
                   sl:float, tl:float, xl:float, yl:float, hl:float, 
                   sr:float, tr:float, xr:float, yr:float, hr:float) -> Tuple[float, float]:
        return self.start_solve(problem, s, t, sl, tl, xl, yl, hl, sr, tr, xr, yr, hr)

    def center_solve(self, problem: HypProblem, s:float, t:float,
                   sl:float, tl:float, xl:float, yl:float, hl:float, 
                   sr:float, tr:float, xr:float, yr:float, hr:float) -> Tuple[float, float]:
        return self.start_solve(problem, s, t, sl, tl, xl, yl, hl, sr, tr, xr, yr, hr)