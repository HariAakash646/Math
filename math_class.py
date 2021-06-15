import turtle

class Math:

    def __init__(self):
        self.log = {"lcm": [], "pair_of_lin_eq_2_vars": [],
                    "hcf": [], "pythagoras_theorem": [],
                    "quadratic_eq": []}
    
    def lcm(self, x, y):
        """This method takes two
        integers and returns the L.C.M."""

        if x > y:
            greater = x
        else:
            greater = y

        while(True):
            if((greater % x == 0) and (greater % y == 0)):
                lcm = greater
                break
            greater += 1
        
        self.log["lcm"].append({"input1": x, "input2": y, "output": lcm})

        return lcm

    def hcf(self, x, y):
        """This method takes two
        integers and returns the H.C.F"""

        lower = min(x, y)

        while True:
            if (x % lower == 0) and (y % lower == 0):
                hcf = lower
                break
            lower -= 1

        self.log["hcf"].append({"input1": x, "input2": y, "output": hcf})
        return hcf


    def dec_to_int(self, a, b, c):
        a_int_ratio = a.as_integer_ratio()
        b_int_ration = b.as_integer_ratio() 
        common_val = self.lcm(a_int_ratio[1], b_int_ration[1])
        return int(a * common_val), int(b * common_val), c * common_val 
    
    def pair_of_linear_eq(self, a1=1, b1=1, c1=1, a2=1, b2=1, c2=1):
        """This method solves a pair of 
        linear equations given the co-efficients"""

        if type(a1) is float or type(b1) is float:
            a1, b1, c1 = self.dec_to_int(a1, b1, c1)
        if type(a2) is float or type(b2) is float:
            a2, b2, c2 = self.dec_to_int(a2, b2, c2)
        a_lcm = self.lcm(abs(a1), abs(a2))
        up_b1 = (a_lcm/a1) * b1
        up_b2 = (a_lcm/a2) * b2
        up_c1 = (a_lcm/a1) * c1
        up_c2 = (a_lcm/a2) * c2
        if a1 and a2 > 0 or a1 and a2 < 0:
            b = up_b1 - up_b2
            c = up_c1 - up_c2
        else:
            b = up_b1 + up_b2
            c = up_c1 + up_c2

        try:
            y = c/b
            x = (c1 - b1 * y) / a1
        except ZeroDivisionError:
            if a1/a2 == b1/b2 == c1/c2:
                print("There are infinite solutions")
            elif a1/a2 == b1/b2 != c1/c2:
                print("There are no solutions")
            return None, None
        return x, y

    @staticmethod
    def coeff_of_var(eq, var):
        """Given an eq and var returns
        the coeffecient of the variable"""

        out_str = ""
        out_sign = 1
        var_idx = eq.index(var)
        for i in range(var_idx-1, -1, -1):
            if eq[i] == '-':
                out_sign = -1
                break
            elif eq[i] == '+':
                break
            elif eq[i] == ' ':
                next
            else:
                out_str += eq[i]
        if '.' in out_str:
            out = float(out_str[::-1]) * out_sign
        elif '/' in out_str:
            temp = out_str[::-1]
            val1, val2 = temp.split('/')
            out = (int(val1)/int(val2)) * out_sign
        else:
            try:
                out = int(out_str[::-1]) * out_sign
            except ValueError:
                out = out_sign

        return out

    

    def str_lin_eq(self, eq):
        """This method returns the coefficients
        of a linear equation given the equation 
        as a string"""

        a = self.coeff_of_var(eq, 'x')
        b = self.coeff_of_var(eq, 'y')
        c_str = ""
        c_sign = 1
        for i in eq[::-1]:
            if i == "-":
                c_sign = -1
                break
            elif i == "+" or i == " " or i == "=":
                break
            else:
                c_str += i
        
        if '.' in c_str:
            c = float(c_str[::-1]) * c_sign
        elif '/' in c_str:
            temp = c_str[::-1]
            val1, val2 = temp.split('/')
            c = (int(val1)/int(val2)) * c_sign
        else:
            c = int(c_str[::-1]) * c_sign
        return a, b, c

    def solve_pair_lin_eq(self, eq1, eq2):
        """This method solves a pair of lin eq
        returning the x and y values
        input format: ax + by = c"""

        a1, b1, c1 = self.str_lin_eq(eq1)
        a2, b2, c2 = self.str_lin_eq(eq2)
        out = self.pair_of_linear_eq(a1, b1, c1, a2, b2, c2)

        self.log["pair_of_lin_eq_2_vars"].append({"eq1": eq1, "eq2": eq2, "x": out[0], "y": out[1]})

        return out

    def quadratic_eq(self, a, b, c):
        """This method returns the possible
        values of x given the coefficients
        of x and x ** 2"""

        x1 = (-b + (b ** 2 - 4*a*c)**(1/2))/(2*a)
        x2 = (-b - (b ** 2 - 4*a*c)**(1/2))/(2*a)

        return x1, x2

    def str_quadratic_eq(self, eq):
        """This method returns the coefficients
        of a quadratic equation given the equation
        as a string"""

        a = self.coeff_of_var(eq, "x^2")
        eq = eq.replace("x^2", '')
        if "x" in eq:
            b = self.coeff_of_var(eq, "x")
        else:
            b = 0
        c = self.coeff_of_var(eq, "=")

        return a, b, c

    def solve_quadratic_eq(self, eq):
        """This method returns the possible
        values of x given a quadratic equation
        input format: ax^2 + bx + c = 0"""
        
        a, b, c = self.str_quadratic_eq(eq)
        x1, x2 = self.quadratic_eq(a, b, c)

        self.log["quadratic_eq"].append({"eq": eq, "possible_x": (x1, x2)})
        return x1, x2

    def pythagoras_theorem(self, a=None, b=None, c=None):
        """Given two any two sides of a 
        right-angled triangle returns
        the third side"""

        if a is None and b is None or b is None and c is None or a is None and c is None:
            print("Give the values of atleast two sides of the triangle")
            raise ValueError
        if a is None:
            a = (c ** 2 - b ** 2) ** (1/2)
        elif b is None:
            b = (c ** 2 - a ** 2) ** (1/2)
        elif c is None:
            c = (a ** 2 + b ** 2) ** (1/2)
        
        self.log["pythagoras_theorem"].append({"a": float(a), "b": float(b), "c": float(c)})
        return float(c), float(b), float(a)

    def __str__(self):
        return str(self.log)


class Trignometry(Math):
    
    def __init__(self):
        super().__init__()
        self.sine = {0: 0, 30: 1/2,
                    45: 1/(2**(1/2)),
                    60: (3 ** (1/2))/2,
                    90:1}
        
        self.cosine = {0: 1,
                        30: (3 ** (1/2))/2,
                        45: 1/(2**(1/2)),
                        60: 1/2, 90:0}
        
        self.tan = {0: 0,
                    30: 1/(3 ** (1/2)),
                    45: 1,
                    60: 3**(1/2),
                    90: None}

        self.sine1 = 0.01745240643728351281941897851632
        self.cosine1 = 0.99984769515639123915701155881391

    @staticmethod
    def draw_triangle(angle1, angle2, angle3, hypotenuse, adjacent, opposite):
        """This method draws the the specified
        triangle using the turtle"""

        board = turtle.Turtle()
 
        board.forward(adjacent) # draw base
        
        board.left(angle2)
        board.forward(opposite)
        
        board.left(90 + angle1)
        board.forward(hypotenuse)

        board.hideturtle()
 
        turtle.done()

    def find_sine_cos_tan(self, angle):
        """This method returnd the cos, sin
        and tan of any given angle"""

        if angle in self.sine.keys():
            return self.sine[angle], self.cosine[angle], self.tan[angle]
        
        poss_angles = [i for i in self.sine.keys() if i < angle]
        max_poss = max(*poss_angles, 0)
        prev_sine = self.sine[max_poss]
        prev_cos = self.cosine[max_poss]
        while max_poss < angle:
            curr_sine = prev_sine * self.cosine1 + prev_cos * self.sine1
            curr_cos = self.cosine1 * prev_cos - prev_sine * self.sine1
            prev_sine = curr_sine
            prev_cos = curr_cos
            max_poss += 1

        return curr_sine, curr_cos, curr_sine/curr_cos

    def trignometric_sides_finder(self, angle, hypotenuse=None, opposite=None, adjacent=None):
        """This method returns all three sides 
        of a triangle given an angle and either
        of the hypotenuse, opposite or adjacent"""

        sine, cos, tan = self.find_sine_cos_tan(angle)
        
        if hypotenuse:
            opposite_out = sine * hypotenuse
            adjacent_out = cos * hypotenuse
            hypotenuse_out = hypotenuse * 1.
        elif opposite:
            hypotenuse_out = opposite / sine
            adjacent_out = opposite / tan
            opposite_out = opposite * 1.
        elif adjacent:
            hypotenuse_out = adjacent / cos
            opposite_out = tan * adjacent
            adjacent_out = adjacent * 1.
        else:
            print("Please enter the value of either one of the sides")
            raise ValueError

        return hypotenuse_out, opposite_out, adjacent_out

    def trignometric_angle_finder(self, hypotenuse=None, opposite=None, adjacent=None):
        """This method returns the angle
        given any two of the hypotenuse,
        opposite or adjacent side lengths"""

        curr_diff = 10000.
        if hypotenuse and opposite:
            val = opposite /hypotenuse
            for i in range(0, 91):
                sin, _, _ = self.find_sine_cos_tan(i)
                if abs(sin - val) < curr_diff:
                    curr_diff = abs(sin-val)
                    best_angle = i
        elif hypotenuse and adjacent:
            val = adjacent / hypotenuse
            for i in range(0, 91):
                _, cos, _ = self.find_sine_cos_tan(i)
                if abs(cos - val) < curr_diff:
                    curr_diff = abs(cos - val)
                    best_angle = i
        elif adjacent and opposite:
            val = opposite / adjacent
            for i in range(0, 90):
                _, _, tan = self.find_sine_cos_tan(i)
                if abs(tan - val) < curr_diff:
                    curr_diff = abs(tan - val)
                    best_angle = i
        
        return best_angle
    
    def trignometric_solver(self, angle=None, hypotenuse=None, opposite=None, adjacent=None, print_out=False, draw=False, drawing_magnific=1):
        """This method returns all the dimensions
        and angles of a right angled triangle given 
        a pair of sides or a side and an angle.
        Can also print out the output.
        Can also draw the triangle."""

        if angle is None:
            angle = self.trignometric_angle_finder(hypotenuse, opposite, adjacent)
            hypotenuse, opposite, adjacent = self.pythagoras_theorem(adjacent, opposite, hypotenuse)
        else:
            hypotenuse, opposite, adjacent = self.trignometric_sides_finder(angle, hypotenuse, opposite, adjacent)
        third_angle = 180 - (90 + angle)
        if print_out:
            print("Triangle ABC:")
            print()
            print(f"Angle A: {angle}")
            print("Angle B: 90")
            print(f"Angle C: {third_angle}")
            print(f"AB: {adjacent}")
            print(f"BC: {opposite}")
            print(f"AC: {hypotenuse}")

        if draw:
            self.draw_triangle(angle, 90, third_angle, hypotenuse*drawing_magnific, adjacent*drawing_magnific, opposite*drawing_magnific)

        return angle, third_angle, hypotenuse, opposite, adjacent


# solver = Math()
# output = solver.solve_pair_lin_eq("x + y = 9", "88x - 11y = 0")
# print(output)
# solver.hcf(20, 8)
# solver.pythagoras_theorem(a=8, b=9)
# solver.solve_quadratic_eq("-6x^2+ 6 = 0")

trig_solver = Trignometry()
# print(trig_solver.trignometric_sides_finder(34, hypotenuse=10))
# print(trig_solver.trignometric_angle_finder(opposite=3, adjacent=4))
# print(trig_solver.find_sine_cos_tan(32))
print("*************************************************************")
trig_solver.trignometric_solver(angle=1, hypotenuse=0, opposite=4, adjacent=None, print_out=True, draw=True, drawing_magnific=5)