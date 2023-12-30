import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import trapz
from scipy.optimize import fsolve
from typing import List, Tuple, Callable



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
    def __init__(self, feed: float, xb: float, xf: float, xd: float, q: float,
                 fxy: Callable, r: float, name: str = ""):
        """__init__ _summary_

        Args:
            feed (float): flowrate of the feed stream
            xb (float): molar fraction of the liquid at the buttoms
            xf (float): molar fraction of the liquid in the feed stream
            xd (float): molar fraction of the liquid at the distilate
            q (float): thermodynamical state of the feed
            fxy (Callable): xy equilibrium
            r (float): reflux ratio
            name (str, optional): name of the object/tower. Defaults to "".
        """
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
        self.n_trays = 0    # initializing the tray numbers
        self.Rmin = 0       # initializing the reflux ratio
        # almost private methods
        self._rmin()        # compute the minimum reflux ratio
        self._azeocheck()   # check if there are any azeotropes
    
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
        """_azeocheck check if there are any azeotropes in the solution.
        distillation process stops at the azeotrope point, if exists.
        """
        self.x_azeo = -1    # initial value
        fun = lambda x: self.f(x) - x
        x = fsolve(fun, 0.5)[0] # solve: f(x) = x where f(x) is the equilibrim
        if (0 < x < 1) and (x > 1e-2 or x < 1-1e-2):
            self.x_azeo = x
        return None
    
    def upper_line(self, x: float):
        """upper_line the function for operating line at the upper section of
        the column.
        R/(R+1)*x + x_d/(R+1)
        latex: y = \frac{R}{R+1}x + \frac{x_{d}}{R+1}

        Args:
            x (float): input parameter, liquid mole fraction

        Returns:
            _type_: vapour mole fraction (y)
        """
        return self.R / (self.R+1) * x + self.x_D / (self.R+1)
    
    def lower_line(self, x: float):
        """lower_line the function for operating line at the lower section of
        the column.

        Args:
            x (float): input parameter, liquid mole fraction

        Returns:
            _type_: vapour mole fraction (y)
        """
        return self._a * x + self._b
    
    def plot(self):
        N = 50 # controls the resolution
        x_lower = np.linspace(self.x_B, self.x_mid, N)
        x_upper = np.linspace(self.x_mid, self.x_D, N)
        x = np.linspace(0, 1, N)
        _ = plt.figure("McCabe-Thiele", figsize = (7, 7))
        # plot title
        plt.title(self.system_name)
        # plotting
        plt.plot(x_lower, self.lower_line(x_lower), "r")    # stripping section
        plt.plot(x_upper, self.upper_line(x_upper), "m")    # rectifying section
        plt.plot([self.x_F, self.x_mid], 
                 [self.x_F, self.upper_line(self.x_mid)], 
                 "k")                                       # q-line
        plt.legend(["Stripping", "Rectifying", "q-line"], fontsize = 8)
        plt.plot(x, self.f(x), "k", lw = 1)                 # equilibrium
        plt.plot(x, x, "k--", lw = 1)                       # y = x
        plt.plot([self.x_B, self.x_B], [0, self.x_B], "b--")
        plt.plot([self.x_D, self.x_D], [0, self.x_D], "b--")
        plt.plot([self.x_F, self.x_F], [0, self.x_F], "b--")
        # axis ticks and labels
        plt.xticks(np.linspace(0, 1, 11))
        plt.yticks(np.linspace(0, 1, 11))
        plt.xlabel("liquid mole-fraction")
        plt.ylabel("vapour mole-fraction")
        # axis type
        plt.axis("square")
        # range of the x, y parameters; between 0 and 1
        plt.xlim(0, 1)
        plt.ylim(0, 1)
        return 0
    
    def run(self):
        # compute the inverse of the linear function
        # y = ax + b ; x = (y - b)/a
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
            plt.annotate(f"{n_trays}", (0.01 + X, (F(X)+Y)/2))
            X = X_
        self.n_trays = n_trays
        # control if there are infinite number of trays calculated which means
        # an azeotrope exists.
        if n_trays < 100:
            return 0
        else:
            return -1
        
    def fensk(self, N: int = 100):
        """fensk Fensk equation: an estimate of the minimum required trays for
        a given system; A rule of thumb.

        Args:
            N (int, optional): number of points in the calculation. 
            Defaults to 100.

        Returns:
            _type_: minimum number of the trays.
        """
        alphaFun = lambda x: (self.f(x)/x)/((1-self.f(x))/(1-x))
        x = np.linspace(1e-6, 1, N, endpoint = False)
        y = alphaFun(x)
        a = trapz(y, x)
        self.average_alpha = a  # average volatility of the solution
        n = np.log(self.x_D/(1 - self.x_D) * (1 - self.x_B)/self.x_B) / np.log(a)
        return n
    
    def _rmin(self):
        """_rmin compute the minimum reflux ratio
        """
        if self.q == 1:
            x = self.x_F
        else:
            self.eq = lambda x: self.f(x)-self.q/(self.q-1)*x+self.x_F/(self.q-1)
            x = fsolve(self.eq, self.x_F)[0]
        y = self.f(x)
        k = (y - self.x_D) / (x - self.x_D)
        self.Rmin = k / (1 - k)
        return None

if __name__ == "__main__":
    # default example
    b = distillColumn(1000, 0.15, 0.65, 0.9, 0.5, lambda x, alpha=2.8: alpha * x / (1 + (alpha - 1) * x), 1)
    b.run()
    plt.show()

