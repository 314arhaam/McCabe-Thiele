# $\textbf{McCabe-Thiele}$
Functional & Object-Oriented programs to rule all distillation columns.  
  
![The most beautiful mccabe code!](https://github.com/314arhaam/McCabe-Thiele/blob/main/img/haskell-mccabe.png)
## Introduction
The [McCabe-Thiele method](https://en.wikipedia.org/wiki/McCabe%E2%80%93Thiele_method) is a technique that is commonly employed in the field of chemical engineering to model the separation of two substances by a distillation column. It uses the fact that the composition at each theoretical tray is completely determined by the mole fraction of one of the two components. This method is based on the assumptions that the distillation column is isobaric - i.e the pressure remains constant - and that the flow rates of liquid and vapor do not change throughout the column (i.e., constant molar overflow). The assumption of constant molar overflow requires that: 
1. The heat needed to vaporize a certain amount of liquid of the feed components are equal,
2. For every mole of liquid vaporized, a mole of vapor is condensed, and
3. Heat effects such as heat needed to dissolve the substance(s) are negligible.  

The method was first published by Warren L. McCabe and Ernest Thiele in 1925, both of whom were working at the Massachusetts Institute of Technology (MIT) at the time.

## Example
### Script
```python
import mccabe

# simple equilibrium function with constant volatility of alpha
# custom function with real-wrold data could be used
equil = lambda x, alpha=2.8: alpha * x / (1 + (alpha - 1) * x)

kwargs = {"feed": 1000, # feed stream molar flow-rate
          "xb": 0.15, # bottoms mole-fraction
          "xf": 0.65, # feed mole-fraction
          "xd": 0.9, # distillate mole-fraction
          "q": 0.5, # q-line slope
          "fxy": equil # equilibrim function}

col = distillColumn(**kwargs)

col.run()
plt.show()
```
### Result
```
"""

ooo        ooooo             .oooooo.              .o8
`88.       .888'            d8P'  `Y8b            "888
 888b     d'888   .ooooo.  888           .oooo.    888oooo.   .ooooo.
 8 Y88. .P  888  d88' `"Y8 888          `P  )88b   d88' `88b d88' `88b
 8  `888'   888  888       888           .oP"888   888   888 888ooo888
 8    Y     888  888   .o8 `88b    ooo  d8(  888   888   888 888    .o
o8o        o888o `Y8bod8P'  `Y8bood8P'  `Y888""8o  `Y8bod8P' `Y8bod8P'






        ---------------------------
        feed        1000.000
        buttom      333.333
        top         666.667
        xB          0.150
        xD          0.900
        xF          0.650

        NO. trays   6
        min RR      0.597
        min. Tray   3.856

        ave. alpha  2.77
        azeotrope   -1.000

"""
```

<p align="center">
  <img src="https://github.com/314arhaam/macaboo/blob/main/img/McCabe-Thiele.png" width="500" title="plot">
</p>
