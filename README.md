# macaboo
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/8eac61e7edde444ba1855ffef2ecbc50)](https://www.codacy.com/gh/314arhaam/macaboo/dashboard?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=314arhaam/macaboo&amp;utm_campaign=Badge_Grade)  
An Object-Oriented program to rule all distillation columns!

## Introduction
[McCabe-Thiele method](https://en.wikipedia.org/wiki/McCabe%E2%80%93Thiele_method), is a simple method for designing distillation columns. It is in the syllabus of Unit Operations I course, in Chemical Engineering.

## Example
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

        NO. trays   7
        min RR      0.597
        min. Tray   3.856

        ave. alpha  2.77
        azeotrope   -1.000
"""
```

<p align="center">
  <img src="https://github.com/314arhaam/macaboo/blob/main/McCabe-Thiele.png" width="500" title="plot">
</p>
