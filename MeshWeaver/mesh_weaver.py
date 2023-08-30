

class MeshWeaver:
    def __init__(self, s, t, c, count_node_s):
        """
        Строит сетку из узлов
        :param s: диапазон абсцисс
        :param t: диапазон ординат
        :param c: углы наклона (первый отрицательный, записывается без знака)
        :param count_node_s: (кол-во узлов по s)
        """
        self.s = s
        self.t = t
        self.c = c
        self.count_node_s = count_node_s

