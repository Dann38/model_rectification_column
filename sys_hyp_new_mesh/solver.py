from mesh import Mesh
if __name__ == '__main__':
    mesh = Mesh()
    mesh.create_mesh(revers_time=True)
    # mesh.solve()
    # for node in mesh.result[-3][-2]:
    #     print(node)
    # mesh.plot_border_right()
    mesh.plot_mesh()
    # mesh.plot_final()