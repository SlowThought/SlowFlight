#lang racket/base
(require "zfoil.rkt"; should provide racket/base, /vector, /math, and chosen zfoil functions
         "airfoils.rkt"
         plot)

(let[(xs(vector-map(λ(z)(real-part z))naca-4412))
     (ys(vector-map(λ(z)(imag-part z))naca-4412))
     (cps(vector-map(λ(v)(-(Cp(magnitude v))))(V (/ pi 18)(make-geometry naca-4412))))
     (vectorize(λ(A B)(map(λ(a b)(vector a b))(vector->list A)(vector->list B))))]
  (plot(list(axes)
            (lines(vectorize xs ys))
            (points(vectorize xs cps)))
       #:x-min -0.1
       #:y-label"-Cp"
       #:title"NACA 4412 @ 10 deg alpha"))