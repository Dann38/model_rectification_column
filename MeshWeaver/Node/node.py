
class Node:
    def __init__(self, s, t, bottom_left_node=None, bottom_right_node=None, top_left_node=None, top_right_node=None):
        self.s = s
        self.t = t
        self.bottom_left_node = bottom_left_node
        self.bottom_right_node = bottom_right_node
        self.top_left_node = top_left_node
        self.top_right_node = top_right_node

    def character_in_t(self, c, t):
        b0 = self.s - c*self.t
        s = c*t + b0
        return s

