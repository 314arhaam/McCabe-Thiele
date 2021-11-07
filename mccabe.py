import os
import numpy as np
from scipy.integrate import trapz
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

os.system("cls")
print(r"""
ooo        ooooo             .oooooo.              .o8                 
`88.       .888'            d8P'  `Y8b            "888                 
 888b     d'888   .ooooo.  888           .oooo.    888oooo.   .ooooo.  
 8 Y88. .P  888  d88' `"Y8 888          `P  )88b   d88' `88b d88' `88b 
 8  `888'   888  888       888           .oP"888   888   888 888ooo888 
 8    Y     888  888   .o8 `88b    ooo  d8(  888   888   888 888    .o 
o8o        o888o `Y8bod8P'  `Y8bood8P'  `Y888""8o  `Y8bod8P' `Y8bod8P' 
                                                                       
                                                                       
                                                                       
""")

class distillColumn:
    def __init__(self, feed, xb, xf, xd, q, fxy, r, name=""):
        self.system_name = name
        #
        self.feed = feed
        self.q = q
        self.f = fxy #lambda x: alpha * x / (1 + (alpha - 1) * x)
        self.x_B, self.x_F, self.x_D = xb, xf, xd
        self.R = r
        #
        self.D = (xf - xb) / (xd - xb) * feed
        self.B = feed - self.D
        self.x_mid = q!=1 and (xd/(r+1)+xf/(q-1))/((q/(q-1))-(r/(r+1))) or xf
        self.y_mid = self.upper_line(self.x_mid)
        #
        self.L_upper = self.D * self.R
        self.V_upper = (1 + self.R) * self.D
        self._a = (self.y_mid - self.x_B) / (self.x_mid - self.x_B)
        self._b = (1 - self._a) * self.x_B
        self.L_lower = self.B / (self._a - 1)
        self.V_lower = self.L_lower - self.B
        #
        self.n_trays = 0
        self.Rmin = 0
        #
        self._rmin()
        self._azeocheck()
        print(self)
    def __repr__(self):
        r = f"""
        {self.system_name}
        ---------------------------
        feed        {self.feed:.3f}
        buttom      {self.B:.3f}
        top         {self.D:.3f}
        xB          {self.x_B:.3f}
        xD          {self.x_D:.3f}
        xF          {self.x_F:.3f}

        NO. trays   {self.n_trays}
        min RR      {self.Rmin:.3f}
        min. Tray   {self.fensk():.3f}

        ave. alpha  {self.average_alpha:.3}
        azeotrope   {self.x_azeo:.3f}
        """
        return r
    def _azeocheck(self):
        self.x_azeo = -1
        fun = lambda x: self.f(x) - x
        x = fsolve(fun, 0.5)[0]
        if (0 < x < 1) and (x > 1e-2 or x < 1-1e-2):
            self.x_azeo = x
        return 0
    def upper_line(self, x):
        return self.R / (self.R+1) * x + self.x_D / (self.R+1)
    def lower_line(self, x):
        return self._a * x + self._b
    def plot(self):
        N = 50
        x_lower = np.linspace(self.x_B, self.x_mid, N)
        x_upper = np.linspace(self.x_mid, self.x_D, N)
        x = np.linspace(0, 1, N)
        fig = plt.figure("McCabe-Thiele", figsize=(7, 7))
        plt.title(self.system_name)
        plt.plot(x_lower, self.lower_line(x_lower), "r")
        plt.plot(x_upper, self.upper_line(x_upper), "m")
        plt.plot([self.x_F, self.x_mid], [self.x_F, self.upper_line(self.x_mid)], "k")
        plt.legend(["Stripping", "Rectifying", "q-line"], fontsize=8)
        plt.plot(x, self.f(x), "k", lw=1)
        plt.plot(x, x, "k--", lw=1)
        plt.plot([self.x_B, self.x_B], [0, self.x_B], "b--")
        plt.plot([self.x_D, self.x_D], [0, self.x_D], "b--")
        plt.plot([self.x_F, self.x_F], [0, self.x_F], "b--")
        plt.xticks(np.linspace(0, 1, 11))
        plt.yticks(np.linspace(0, 1, 11))
        plt.xlabel("liquid mole-fraction")
        plt.ylabel("vapour mole-fraction")
        plt.axis("square")
        plt.xlim(0, 1)
        plt.ylim(0, 1)
        return 0
    def run(self):
        linv = lambda fun: lambda y: (y - fun(0))/(fun(1) - fun(0))
        self.plot()
        X = self.x_B
        F = self.lower_line
        n_trays = -1
        while X < self.x_D and n_trays < 100:
            Y = self.f(X)
            plt.plot([X, X], [F(X), Y], "k", lw=0.5)
            if linv(F)(Y) > self.x_mid:
                F = self.upper_line
            X_ = linv(F)(Y)
            if linv(F)(Y) > self.x_D:
                plt.plot([X, Y], [Y, Y], "k", lw=0.5)
            else:
                plt.plot([X, X_], [Y, Y], "k", lw=0.5)
            n_trays += 1
            #plt.annotate(f"{n_trays}", (0.01 + X, (F(X)+Y)/2))
            X = X_
        self.n_trays = n_trays
        if n_trays < 100:
            return 0
        else:
            return -1
    def fensk(self, N=100):
        alphaFun = lambda x: (self.f(x)/x)/((1-self.f(x))/(1-x))
        x = np.linspace(1e-6, 1, N, endpoint=False)
        y = alphaFun(x)
        a = trapz(y, x)
        self.average_alpha = a
        n = np.log(self.x_D/(1 - self.x_D) * (1 - self.x_B)/self.x_B) / np.log(a)
        return n
    def _rmin(self):
        if self.q == 1:
            x = self.x_F
        else:
            self.eq = lambda x: self.f(x)-self.q/(self.q-1)*x+self.x_F/(self.q-1)
            x = fsolve(self.eq, self.x_F)[0]
        y = self.f(x)
        k = (y - self.x_D) / (x - self.x_D)
        self.Rmin = k / (1 - k)
        return 0

if __name__ == "__main__":
    # default example
    b = distillColumn(1000, 0.15, 0.65, 0.9, 0.5, lambda x, alpha=2.8: alpha * x / (1 + (alpha - 1) * x), 1)
    b.run()
    plt.show()

