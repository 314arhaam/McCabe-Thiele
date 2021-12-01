-- ghc --make distill.hs
-- distill.exe (windows)

-- McCabe-Thiele method q = 1

{- 
operating line as a function of
    xb  buttoms liquid mole fraction
    xd  distillate liquid mole fraction
    zf  feed concentration
    r   reflux ratio
-}
operationLine :: Float -> Float -> Float -> Float -> (Float -> Float)
operationLine xb xd zf r = \x -> case () of
                               _| x >  xm -> a1 * x + b1 -- stripping
                                | x <= xm -> a2 * x + b2 -- rectifying
                                where
                                    xm = zf
                                    a1 = r / (r + 1)
                                    b1 = xd / (r + 1)
                                    a2 = (a1 * xm + b1 - xb) / (xm - xb)
                                    b2 = xb - a2 * xb

{-
solves distillation column for
    x           recursive parameter, mole fraction
    oLine       combined operating lines as a bijective function
    invEquil    inversed of the equilibrium function
-}
distil :: Float -> (Float -> Float) -> (Float -> Float) -> [Float]
distil x oLine invEquil | x_ >  oLine x_ = [x]
                        | x_ <= oLine x_ = x : distil (oLine x_) oLine invEquil
                        where
                            x_ = invEquil x


main :: IO ()
main = do
    let xb  = 0.1 -- buttoms molar concentration
        xd  = 0.9 -- tops molar concentration
        zf  = 0.5 -- feed concentration
        rr  = 0.6 -- reflux ratio
        al  = 6.0 -- relative volatility (alpha)
    let inverseEquilibrium = \y -> y / (al - (al - 1) * y) -- invEq
    let operatingLinesFun = operationLine xb xd zf rr -- apply McCabe
    let trays = distil xd operatingLinesFun inverseEquilibrium
    print trays
    print("Number of Trays", length trays - 1)
