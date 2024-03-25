from .problem import HypProblem
from .mesh import Mesh
from .solve_method import SolveMethod

class HypSolver:
    def __init__(self, mesh: Mesh, solve_method: SolveMethod) -> None:
        self.solve_method = solve_method
        self.mesh = mesh
        

    def solve(self, problem: HypProblem) -> dict:
        self.mesh.solve(problem, self.solve_method)
        return {}